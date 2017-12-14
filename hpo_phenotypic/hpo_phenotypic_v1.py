#!/usr/bin/env python
# --coding:utf-8--
# date: 2017/12/12
# author:wuling
# emai:ling.wu@myhealthgene.com

#this model set  to xxxxxx

import sys
sys.path.append('../')
sys.setdefaultencoding = ('utf-8')
from share import *
from config import *  

__all__ = ['downloadData','extractData','standarData','insertData','updateData','selectData']

version  = 1.0

model_name = psplit(os.path.abspath(__file__))[1]

(hpo_phenotypic_model,hpo_phenotypic_raw,hpo_phenotypic_store,hpo_phenotypic_db,hpo_phenotypic_map) = buildSubDir('hpo_phenotypic')

log_path = pjoin(hpo_phenotypic_model,'hpo_phenotypic.log')

# main code
def downloadData( redownload=False ):
    '''
    this function is to download the raw data from go gene FTP WebSite
    args:
    redownload-- default False, check to see if exists an old edition before download
                       -- if set to true, download directly with no check
    rawdir-- the directory to save download file
    '''
    if  not redownload:

        (choice,existhpoFile) = lookforExisted(phenotypic_raw,'phenotypic')

        if choice != 'y':
            return

    if redownload or not existhpoFile or  choice == 'y':

        rawdir = pjoin(hpo_phenotypic_raw,'phenotypic_{}'.format(today))

        createDir(rawdir)

        process = hpo_parser(today)
    
        mt = process.getMt()   

        # for  url in hpo_download_urls:

        #     process.wget(url,mt,rawdir)

        func = lambda x:process.wget(x,mt,rawdir)

        multiProcess(func,hpo_download_urls,size=3)

    if not os.path.exists(log_path):

        with open(log_path,'w') as wf:
            json.dump({'hpo_phenotypic':[(mt,today,model_name)]},wf,indent=8)

    print  'datadowload completed !'

    filapths = [pjoin(rawdir,filename) for filename in listdir(rawdir)]

    return (filapths,today)

def extractData(filepaths,version):
    
    process = hpo_parser(version)

    for filepath in filepaths:

        if filepath.endswith('obo'):

            process.obo(filepath)

    for filepath in filepaths:

        if filepath.endswith('txt'):

            process.tsv(filepath)

        elif filepath.endswith('.tab'):

            process.tab(filepath)

    print 'extract and insert complete '

    return (filepaths,version)

def updateData(insert=False,_mongodb='../_mongodb/'):

    hpo_phenotypic_log = json.load(open(log_path))

    process = hpo_parser(today)

    mt = process.getMt()

    if mt != hpo_phenotypic_log['hpo_phenotypic'][-1][0]:

        filepath,version = downloadData(redownload=True)

        extractData(filepath,version)

        hpo_phenotypic_log['hpo_phenotypic'].append((mt,today,model_name))

        # create new log
        with open(log_path,'w') as wf:

            json.dump(hpo_phenotypic_log,wf,indent=2)

        print  '{} \'s new edition is {} '.format('hpo_phenotypic',mt)

        bakeupCol('hpo_phenotypic_{}'.format(version),'hpo_phenotypic',_mongodb)

    else:
        print  '{} {} is the latest !'.format('hpo_phenotypic',mt)

def selectData(querykey = 'id',value='HP:0011220'):
    '''
    this function is set to select data from mongodb
    args:
    querykey -- a specified field in database
    queryvalue -- a specified value for a specified field in database
    '''
    conn = MongoClient('127.0.0.1', 27017 )

    db = conn.mydb

    colnamehead = 'hpo_phenotypic_'

    dataFromDB(db,colnamehead,querykey,queryvalue=None)

class dbMap(object):

    #class introduction

    def __init__(self,version):

        self.version = version

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        colname = 'hpo_phenotypic_{}'.format(version)

        col = db.get_collection(colname)

        self.col = col

        self.docs = col.find({})

        self.colname = colname

    def maphpoid2geneid(self):
            
        hpoid2geneid =dict()
        hpoid2diseaseid = dict()
        diseaseid2hpoid = {'OMIM':{},'ORPHA':{},'DECIPHER':{}}
        diseaseid2geneid= {'OMIM':{},'ORPHA':{}}
        geneid2disease= dict()

        for doc in self.docs:

            hpoid = doc.get('id')

            ORPHA = doc.get('ORPHA')

            OMIM = doc.get('OMIM')

            DECIPHER = doc.get('DECIPHER')

            if  OMIM:

                diseaseid = OMIM.keys()

#----------------diseaseid2hpoid and hpoid2diseaseid-------------------

                for disid in diseaseid:

                    if disid not in diseaseid2hpoid['OMIM']:

                        diseaseid2hpoid['OMIM'][disid] = list()

                    diseaseid2hpoid['OMIM'][disid].append(hpoid)

                if hpoid not in hpoid2diseaseid:

                    hpoid2diseaseid[hpoid] = dict()

                hpoid2diseaseid[hpoid]['OMIM'] = diseaseid

#----------------diseaseid2geneid and hpoid2geneid and geneid2disease-------------------

                for omimid,val in OMIM.items():

                    gene_list = val.get('gene')

                    if not gene_list:
                        continue

                    for gene in gene_list:

                        geneid = gene.get('gene-id(entrez)')
                        genesym = gene.get('gene-symbol')

                        if omimid not in  diseaseid2geneid['OMIM']:
                            diseaseid2geneid['OMIM'][omimid] = list()

                        diseaseid2geneid['OMIM'][omimid].append(geneid)

                        if hpoid not in hpoid2geneid:

                            hpoid2geneid[hpoid] = list()

                        hpoid2geneid[hpoid].append(geneid)

                        if geneid not in geneid2disease:

                            geneid2disease[geneid] = dict()

                        if 'OMIM' not  in geneid2disease[geneid]:

                            geneid2disease[geneid]['OMIM'] = list()

                        geneid2disease[geneid]['OMIM'].append(omimid)


            if  ORPHA:

                diseaseid = ORPHA.keys()

#----------------diseaseid2hpoid and hpoid2diseaseid-------------------

                for disid in diseaseid:

                    if disid not in diseaseid2hpoid['ORPHA']:

                        diseaseid2hpoid['ORPHA'][disid] = list()

                    diseaseid2hpoid['ORPHA'][disid].append(hpoid)

                if hpoid not in hpoid2diseaseid:

                    hpoid2diseaseid[hpoid] = dict()

                hpoid2diseaseid[hpoid]['ORPHA'] = diseaseid

#----------------diseaseid2geneid and hpoid2geneid and geneid2disease-------------------

                for orphaid,val in ORPHA.items():

                    gene_list = val.get('gene')

                    if not gene_list:
                        continue

                    for gene in gene_list:

                        geneid = gene.get('gene-id(entrez)')
                        genesym = gene.get('gene-symbol')

                        if orphaid not in  diseaseid2geneid['ORPHA']:
                            diseaseid2geneid['ORPHA'][orphaid] = list()

                        diseaseid2geneid['ORPHA'][orphaid].append(geneid)

                        if hpoid not in hpoid2geneid:

                            hpoid2geneid[hpoid] = list()

                        hpoid2geneid[hpoid].append(geneid)

                        if geneid not in geneid2disease:

                            geneid2disease[geneid] = dict()

                        if 'ORPHA' not  in geneid2disease[geneid]:

                            geneid2disease[geneid]['ORPHA'] = list()

                        geneid2disease[geneid]['ORPHA'].append(orphaid)

            if  DECIPHER:

                diseaseid = DECIPHER.keys()

                for disid in diseaseid:

                    if disid not in diseaseid2hpoid['DECIPHER']:

                        diseaseid2hpoid['DECIPHER'][disid] = list()

                    diseaseid2hpoid['DECIPHER'][disid].append(hpoid)

                if hpoid not in hpoid2diseaseid:

                    hpoid2diseaseid[hpoid] = dict()

                hpoid2diseaseid[hpoid]['DECIPHER'] = diseaseid
        
        geneid2hpoid = value2key(hpoid2geneid)

        map_dir = pjoin(hpo_phenotypic_map,self.colname)

        createDir(map_dir)

        save = {'hpoid2diseaseid':hpoid2diseaseid,'diseaseid2hpoid':diseaseid2hpoid,
                        'diseaseid2geneid':diseaseid2geneid,'geneid2disease':geneid2disease,
                        'hpoid2geneid':hpoid2geneid,'geneid2hpoid':geneid2hpoid}

        for name,dic in save.items():

            if name in ['diseaseid2geneid','diseaseid2hpoid','geneid2disease']:
                dic = dedupDicVal(dic)

            elif name in ['geneid2hpoid']:
                for key,val in dic.items():
                    dic[key] = list(set(val))
                    
            with open(pjoin(map_dir,'{}.json'.format(name)),'w') as wf:
                json.dump(dic,wf,indent=2)

    def mapping(self):

        self.maphpoid2geneid()

class hpo_parser(object):
    """docstring for hpo_parser"""
    def __init__(self, version):

        super(hpo_parser, self).__init__()

        self.version = version

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        col = db.get_collection('hpo_phenotypic_{}'.format(version))

        self.col = col

    def getMt(self):

        headers = {'User-Agent':'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Ubuntu Chromium/37.0.2062.120 Chrome/37.0.2062.120 Safari/537.36'}

        web = requests.get(hpo_download_web,headers = headers,verify=False)

        soup = bs(web.content,'lxml')

        h1 = soup.find(name='h1',attrs={'class':'build-caption page-headline'})

        mt = h1.text.strip().split('\n')[1].strip().split('AM')[0].split('(')[1].strip().replace(':','').replace(' ','').replace(',','')

        return mt  
    
    def wget(self,url,mt,rawdir):

        filename = url.rsplit('/',1)[1].strip()

        if filename.endswith('obo'):

            savename = '{}_{}_{}.obo'.format(filename.split('.obo',1)[0].strip(),mt.replace('.','*'),today)

        elif filename.endswith('.tab'):

            savename = '{}_{}_{}.tab'.format(filename.split('.tab',1)[0].strip(),mt.replace('.','*'),today)

        else:
            savename = '{}_{}_{}.txt'.format(filename.split('.txt',1)[0].strip(),mt.replace('.','*'),today)

        storefilepath = pjoin(rawdir,savename)

        command = 'wget -O {} {}'.format(storefilepath,url)

        os.popen(command)

        return storefilepath

    def obo(self,filepath):

        obofile = open(filepath)

        # skip  term head 
        n = 1 

        for line in obofile:

            if line.count('[Term]'):

                break

            n += 1
        
        aset = dict()

        for line in obofile:

            if line.count('[Term]'):

                if aset:

                    self.col.insert(aset)
                     
                    print 'hpo_phenotypic line','obo',n

                aset = dict()

            else:

                line = line.strip()

                if  bool(line):

                    (key,val) = tuple(line.strip().split(':',1))

                    key = key.strip()
                    val = val.strip()

                    if key in ['name','id']:

                        aset[key] = val

                    else:

                        if key not in aset:

                            aset[key] = list()

                        aset[key].append(val)

            n += 1

        if aset:

            self.col.insert(aset)
                     
            print 'hpo_phenotypic line','obo',n

        print 'hpo_phenotypic completed! '

    def tsv(self,filepath):

        tsvfile = open(filepath)

        keys = ['diseaseId','gene-symbol','gene-id(entrez)','HPO-ID','HPO-term-name']

        n = 0 

        for line in tsvfile:

            if line.startswith('#'):
                continue

            data = line.strip().split('\t')

            dic = dict([(key,val) for key,val in zip(keys,data)])

            hpo_id = dic.pop('HPO-ID')

            diseaseId = dic.pop('diseaseId')

            dic.pop('HPO-term-name')

            db = diseaseId.split(':')[0].strip()
            _id = diseaseId.split(':')[1].strip()

            self.col.update(
                {'id':hpo_id},
                {'$push':{'{}.{}.gene'.format(db,_id):dic}},
                )

            print 'hpo_phenotypic txt line',db,diseaseId,n,hpo_id

            n += 1

    def tab(self,filepath):

        tabfile = open(filepath)

        keys = ['db','diseaseId','DB_Name','Qualifier','HPO-ID','DB:Reference','Evidence_code','Onset_modifier','Frequency_modifier','with','Aspect','Synonym','Date','Assigned_by']

        n = 0 

        for line in tabfile:

            if line.startswith('#'):
                continue

            data = line.strip().split('\t')

            dic = dict([(key,val) for key,val in zip(keys,data)])

            db = dic.pop('db')

            hpo_id = dic.pop('HPO-ID')

            diseaseId = dic.pop('diseaseId')

            Assigned_by = dic.pop('Assigned_by')

            dic.pop('with')

            dic['Assigned_by'] = Assigned_by.strip().split(';')

            self.col.update(
                {'id':hpo_id},
                {'$push':{'{}.{}.annotation'.format(db,diseaseId):dic}},
                )

            print 'hpo_phenotypic tab line',db,diseaseId,n,hpo_id

            n += 1

def main():

    modelhelp = model_help.replace('&'*6,'HPO_Phenotypic').replace('#'*6,'hpo_phenotypic')

    funcs = (downloadData,extractData,updateData,selectData,dbMap,hpo_phenotypic_store)

    getOpts(modelhelp,funcs=funcs)
        
if __name__ == '__main__':
    main()
    # downloadData( redownload=True)
    rawdir = '/home/user/project/dbproject/mydb_v1/hpo_phenotypic/dataraw/phenotypic_171212143159'
    filepaths = [pjoin(rawdir,filename) for filename in listdir(rawdir)]
    extractData(filepaths,'171212143159')

    man = dbMap('171212143159')
    man.mapping()
    # 