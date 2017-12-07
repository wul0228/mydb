#!/usr/bin/env python
# --coding:utf-8--
# date: xxxxxx
# author:xxxxxx
# emai:xxxxxx

#this model set  to xxxxxx

import sys
sys.path.append('../')
sys.setdefaultencoding = ('utf-8')
from share import *
from config import *  

__all__ = ['downloadData','extractData','standarData','insertData','updateData','selectData']

version  = 1.0

model_name = psplit(os.path.abspath(__file__))[1]

(clinVar_variant_model,clinVar_variant_raw,clinVar_variant_store,clinVar_variant_db,clinVar_variant_map) = buildSubDir('clinVar_variant')

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
    log_path = pjoin(clinVar_variant_model,'clinVar_variant.log')

    if not os.path.exists(pjoin(clinVar_variant_model,'clinVar_variant.log')):

        with open('./clinVar_variant.log','w') as wf:

            json.dump({'clinVar_variant':[(mt,today,model_name)]},wf,indent=8)

    print  'datadowload completed !'

    return (save_file_path,mt)

def extractData(filepath,version):

    process = clinVar_parser(version)

    process.sum(filepath)

    print 'extract and insert completed'

def updateData():

    clinVar_variant_log = json.load(open('./clinVar_variant.log'))

    rawdir = pjoin(clinVar_variant_raw,'pathway_update_{}'.format(today))

    ftp = connectFTP(**clinVar_varient_ftp_infos)

    filename = clinVar_varient_filename

    mt =  ftp.sendcmd('MDTM {}'.format(filename)).replace(' ','')

    if mt != clinVar_variant_log['clinVar_variant'][-1][0]:

        createDir(rawdir)

        filepath,version = downloadData(redownload=True)

        extractData(filepath,version)

        clinVar_variant_log['clinVar_variant'].append((mt,today,model_name))

        print  '{} \'s new edition is {} '.format('clinVar_variant',mt)
        # create new log
        with open('./clinVar_variant.log','w') as wf:

            json.dump(clinVar_variant_log,wf,indent=2)

    else:

        print  '{} is the latest !'.format('clinVar_variant')


def selectData(querykey = 'GeneID',value='1'):
    '''
    this function is set to select data from mongodb
    args:
    querykey -- a specified field in database
    queryvalue -- a specified value for a specified field in database
    '''
    conn = MongoClient('127.0.0.1', 27017 )

    db = conn.mygene

    colnamehead = 'clinvar_variant_'

    dataFromDB(db,colnamehead,querykey,queryvalue=None)

class dbMap(object):

    #class introduction

    def __init__(self):

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mygene')

        col = db.get_collection('clinvar_variant_{}'.format(version))

        self.col = col

        self.version = version

    def mapGene2AlleleID(self):

        docs = self.col.find({})

        geneid2alleleid = dict()

        genesym2alleleid = dict()

        for doc in docs:

            gene_id = doc.get('GeneID')

            gene_sym = doc.get('GeneSymbol')

            AlleleID = doc.get('AlleleID')

            if gene_id and gene_id not in geneid2alleleid:

                geneid2alleleid[gene_id] = list()

            geneid2alleleid[gene_id].append(AlleleID)

            if gene_sym and gene_sym not in genesym2alleleid:

                genesym2alleleid[gene_sym] = list()

            genesym2alleleid[gene_sym].append(AlleleID)


        with open(pjoin(clinVar_varient_map,'geneid2alleleid.json'),'w') as wf:
            json.dump(geneid2alleleid,wf,indent=2)

        with open(pjoin(clinVar_varient_map,'genesym2alleleid.json'),'w') as wf:
            json.dump(genesym2alleleid,wf,indent=2)


    def mapping(self):

        self.mapXX2XX()

class clinVar_parser(object):

    """docstring for clinVar_parser"""
    def __init__(self, version):

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mygene')

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

                PhenotypeList_list = [ i.strip() for i in PhenotypeList.split(',') if i]
                PhenotypeIDS_list = [ i.strip() for i in PhenotypeIDS.split(',') if i]

                if PhenotypeIDS:

                    dic.update({'PhenotypeIDS':PhenotypeIDS_list})

                if PhenotypeList:

                    dic.update({'PhenotypeList':PhenotypeList_list})

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
    # downloadData(redownload = True)

    # filepath = '/home/user/project/dbproject/mygene_v1/clinVar_variant/dataraw/variant_summary_21320171204103823_171206113223.txt'
    # version = '171206113223'
    # extractData(filepath,version)
    # print 'clinVar_variant'.upper()