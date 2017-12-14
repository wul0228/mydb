#!/usr/bin/env python
# --coding:utf-8--
# date: 2017/12/08
# author:wuling
# emai:ling.wu@myhealthgene.com

#this model set  to download,extract,standard insert and select gene data from clinVar_variant

import sys
sys.path.append('../')
sys.setdefaultencoding = ('utf-8')
from share import *
from config import *  

__all__ = ['downloadData','extractData','standarData','insertData','updateData','selectData']

version  = 1.0

model_name = psplit(os.path.abspath(__file__))[1]

(clinVar_variant_model,clinVar_variant_raw,clinVar_variant_store,clinVar_variant_db,clinVar_variant_map) = buildSubDir('clinVar_variant')

log_path = pjoin(clinVar_variant_model,'clinVar_variant.log')

# main code
def downloadData(redownload = False):
    '''
    this function is to download the raw data from go gene FTP WebSite
    args:
    redownload-- default False, check to see if exists an old edition before download
                       -- if set to true, download directly with no check
    rawdir-- the directory to save download file
    '''
    if  not redownload:

        (choice,existgoFile) = lookforExisted(clinVar_variant_raw,'variant')

        if choice != 'y':
            return

    if redownload or not existgoFile or  choice == 'y':

        ftp = connectFTP(**clinVar_varient_ftp_infos)

        filename = clinVar_varient_filename

        mt =  ftp.sendcmd('MDTM {}'.format(filename)).replace(' ','')

        savefilename = '{}_{}_{}.txt.gz'.format(filename.rsplit('.txt',1)[0].strip(),mt,today)

        remoteabsfilepath = pjoin(clinVar_varient_ftp_infos['logdir'],'{}'.format(filename))

        save_file_path = ftpDownload(ftp,filename,savefilename,clinVar_variant_raw,remoteabsfilepath)

        # gunzip file
        gunzip = 'gunzip {}'.format(save_file_path)

        os.popen(gunzip)

    # create log file
    if not os.path.exists(log_path):

        with open(log_path,'w') as wf:

            json.dump({'clinVar_variant':[(mt,today,model_name)]},wf,indent=8)

    print  'datadowload completed !'

    filepath = save_file_path.split('.gz')[0].strip()

    return (filepath,today)

def extractData(filepath,version):

    process = clinVar_parser(version)

    process.sum(filepath)

    print 'extract and insert completed'

    return (filepath,version)

def updateData(insert=False,_mongodb='../_mongodb/'):

    clinVar_variant_log = json.load(open(log_path))

    ftp = connectFTP(**clinVar_varient_ftp_infos)

    filename = clinVar_varient_filename

    mt =  ftp.sendcmd('MDTM {}'.format(filename)).replace(' ','')

    if mt != clinVar_variant_log['clinVar_variant'][-1][0]:

        filepath,version = downloadData(redownload=True)

        extractData(filepath,version)

        clinVar_variant_log['clinVar_variant'].append((mt,today,model_name))

        # create new log
        with open(log_path,'w') as wf:

            json.dump(clinVar_variant_log,wf,indent=2)

        print  '{} \'s new edition is {} '.format('clinVar_variant',mt)
        
        bakeupCol('clinvar_variant_{}'.format(version),'clinvar_variant',_mongodb)
        
    else:

        print  '{} {} is the latest !'.format('clinVar_variant',mt)

def selectData(querykey = 'GeneID',value='1'):
    '''
    this function is set to select data from mongodb
    args:
    querykey -- a specified field in database
    queryvalue -- a specified value for a specified field in database
    '''
    conn = MongoClient('127.0.0.1', 27017 )

    db = conn.mydb

    colnamehead = 'clinvar_variant_'

    dataFromDB(db,colnamehead,querykey,queryvalue=None)

class dbMap(object):

    #class introduction

    def __init__(self,version):

        self.version = version

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        colname = 'clinvar_variant_{}'.format(self.version)

        col = db.get_collection(colname)

        self.col = col

        self.docs = col.find({})

        self.colname = colname

    def mapGene2AlleleID(self):

        geneid2gensym = dict()

        geneid2alleleid = dict()

        genesym2alleleid = dict()

        for doc in self.docs:

            gene_id = doc.get('GeneID')
            gene_sym = doc.get('GeneSymbol')
            AlleleID = doc.get('AlleleID')

            if gene_id and gene_id != '-1':

                if gene_id not in geneid2gensym :

                    geneid2gensym[gene_id] = list()

                if gene_id not in geneid2alleleid :

                    geneid2alleleid[gene_id] = list()

                geneid2alleleid[gene_id].append(AlleleID)

                if gene_sym:

                    geneid2gensym[gene_id] += gene_sym

                    for sym in gene_sym:

                        if  sym not in genesym2alleleid:

                            genesym2alleleid[sym] = list()

                        genesym2alleleid[sym].append(AlleleID)

        genesym2geneid = value2key(geneid2gensym)

        map_dir = pjoin(clinVar_variant_map,self.colname)

        createDir(map_dir)

        save = {'geneid2gensym':geneid2gensym,'genesym2geneid':genesym2geneid,
                      'geneid2alleleid':geneid2alleleid, 'genesym2alleleid':genesym2alleleid }

        for name,dic in save.items():

            dedupdic = dict()

            for key,val in dic.items():

                dedupdic[key] = list(set(val))

            with open(pjoin(map_dir,'{}.json'.format(name)),'w') as wf:
                json.dump(dedupdic,wf,indent=2)

    def mapping(self):

        self.mapGene2AlleleID()

        # n = 0

        # for doc in self.docs:

        #     AlleleID = doc.get('AlleleID')

        #     gene_sym = doc.get('GeneSymbol')

        #     gene_sym = [i.strip() for i in gene_sym.split(';') if i ]

        #     self.col.update(
        #         {'AlleleID':AlleleID},
        #         {'$set':{'GeneSymbol':gene_sym}}
        #         )

        #     n += 1

        #     print n,AlleleID

class clinVar_parser(object):

    """docstring for clinVar_parser"""
    def __init__(self, version):

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        col = db.get_collection('clinvar_variant_{}'.format(version))

        self.col = col

        self.version = version

    def sum(self,filepath):

        tsvfile = open(filepath)

        n = 0

        for line in tsvfile:

            if line.startswith('#'):

                keys = line.replace('#','',1).strip().split('\t')

                print keys

            else:

                data = [i.strip() for i in line.strip().split('\t')]

                dic = dict([(key,val) for key,val in zip(keys,data)])
            
                PhenotypeIDS = dic.pop('PhenotypeIDS')

                PhenotypeList = dic.pop('PhenotypeList')

                GeneSymbol = dic.pop('GeneSymbol')

                PhenotypeList_list = [ i.strip() for i in PhenotypeList.split(',') if i]
                PhenotypeIDS_list = [ i.strip() for i in PhenotypeIDS.split(',') if i]
                # GeneSymbol_list = [ i.strip() for i in GeneSymbol.split(',') if i]
                GeneSymbol_list = [ i.strip() for i in GeneSymbol.split(';') if i]

                dic.update({
                    'PhenotypeIDS':PhenotypeIDS_list,
                    'PhenotypeList':PhenotypeList_list,
                    'GeneSymbol':GeneSymbol_list
                    })

                location_keys =  ['Assembly','ChromosomeAccession','Chromosome','Start','Stop', 'Cytogenetic','ReferenceAllele','AlternateAllele']
                
                location = dict()

                for key in location_keys:

                    val = dic.pop(key)

                    location.update({key:val})

                AlleleID = dic.pop('AlleleID')

                self.col.update(
                    {'AlleleID':AlleleID},
                    {'$set':dic,
                     '$push':{'location':location}
                     },
                     True
                     )

                print 'line',n,AlleleID

            n += 1

def main():

    modelhelp = model_help.replace('&'*6,'CLINVAR_VARIANT').replace('#'*6,'clinVar_variant')

    funcs = (downloadData,extractData,updateData,selectData,dbMap,clinVar_variant_store)

    getOpts(modelhelp,funcs=funcs)
        
if __name__ == '__main__':
    
    main()
    # filepath,version = downloadData(redownload = True)

    # filepath = '/home/user/project/dbproject/mydb_v1/clinVar_variant/dataraw/variant_summary_21320171204103823_171208183634.txt'
    # version = '171208183634'
    # extractData(filepath,version)
    man = dbMap('171208183634')
    man.mapping()