#!/usr/bin/env python
# --coding:utf-8--
# date: 2017/12/07
# author:wuling
# emai:ling.wu@myhealthgene.com

#this model set  to download,extract,standard insert and select gene data from disgenet

import sys
sys.path.append('../')
sys.setdefaultencoding = ('utf-8')
from share import *
from config import *  

__all__ = ['downloadData','extractData','standarData','insertData','updateData','selectData']

version  = 1.0

model_name = psplit(os.path.abspath(__file__))[1]

(disgenet_disease_model,disgenet_disease_raw,disgenet_disease_store,disgenet_disease_db,disgenet_disease_map) = buildSubDir('disgenet_disease')

log_path = pjoin(disgenet_disease_model,'disgenet_disease.log')

# main code
def downloadData(redownload=False):
    '''
    this function is to download the raw data from go gene FTP WebSite
    args:
    redownload-- default False, check to see if exists an old edition before download
                       -- if set to true, download directly with no check
    rawdir-- the directory to save download file
    '''
    if  not redownload:

        (choice,existgoFile) = lookforExisted(disgenet_disease_raw,'all_gene_disease')

        if choice != 'y':
            return

    if redownload or not existgoFile or  choice == 'y':

        process = disgenet_parser(today)
    
        mt = process.getMt()   

        filepath = process.wget(disgenet_download_url,mt,disgenet_disease_raw)

    if not os.path.exists(log_path):

        with open(log_path,'w') as wf:
            json.dump({'disgenet':[(mt,today,model_name)]},wf,indent=8)

    print  'datadowload completed !'

    return (filepath,today)

def extractData(filepath,version):
    
    process = disgenet_parser(version)

    process.tsv(filepath)
    
    print 'extract and insert completed'
    
    return (filepath,version)

def updateData():

    disgenet_disease_log = json.load(open(log_path))

    process = disgenet_parser(today)

    mt = process.getMt()

    if mt != disgenet_disease_log['disgenet_disease'][-1][0]:

        filepath,version = downloadData(redownload=True)

        extractData(filepath,version)

        disgenet_disease_log['disgenet_disease'].append((mt,today,model_name))

        # create new log
        with open(log_path,'w') as wf:

            json.dump(disgenet_disease_log,wf,indent=2)

        print  '{} \'s new edition is {} '.format('disgenet_disease',mt)

        bakeupCol('disgenet_disease_{}'.format(version),'disgenet_disease')

    else:
        print  '{} {} is the latest !'.format('disgenet_disease',mt)

def selectData(querykey = 'geneId',value='1'):
    '''
    this function is set to select data from mongodb
    args:
    querykey -- a specified field in database
    queryvalue -- a specified value for a specified field in database
    '''
    conn = MongoClient('127.0.0.1', 27017 )

    db = conn.mydb

    colnamehead = 'disgenet_disease_'

    dataFromDB(db,colnamehead,querykey,queryvalue=None)


class dbMap(object):

    #class introduction

    def __init__(self,version):

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        col = db.get_collection('disgenet_disease_{}'.format(version))

        self.col = col

        self.version = version

        self.docs = self.col.find({})

    def mapGeneID2Dis(self):

        (id2dis,sym2dis,dis2id,dis2sym) = ({},{},{},{})

        for doc in  self.docs:

            gene_id = doc.pop('geneId')

            gene_sym = doc.pop('geneSymbol')

            doc.pop('_id')

            disease = doc.keys()

            if gene_id not in id2dis:

                id2dis[gene_id] = list()

            id2dis[gene_id] +=  disease

            if gene_sym not in sym2dis:

                sym2dis[gene_sym] = list()

            sym2dis[gene_sym] +=  disease

        dis2id = value2key(id2dis)

        dis2sym = value2key(sym2dis)

        save = {'id2dis':id2dis,'sym2dis':sym2dis,'dis2id':dis2id,'dis2sym':dis2sym}

        for filename in save.keys():

            print filename,len(save[filename])

            with open(pjoin(disgenet_disease_map,'{}_{}.json'.format(filename,self.version)),'w') as wf:
                json.dump(save[filename],wf,indent=8)

    def mapping(self):

        self.mapGeneID2Dis()

class disgenet_parser(object):

    def __init__(self, version):

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        col = db.get_collection('disgenet_disease_{}'.format(version))

        self.col = col

        self.version = version

    def getMt(self):

        web = requests.get(disgenet_download_web).content

        soup = bs(web,'lxml')

        version = soup.select('#content > div > div > div > div.span8 > div > ul > li:nth-of-type(1)')[0].text.split('version')[1].split(')')[0].strip()

        return version


    def wget(self,url,mt,rawdir):

        filename = url.rsplit('/',1)[1].strip().replace('.tsv.gz','')

        savename = '{}_{}_{}.tsv.gz'.format(filename,mt.replace('.','*'),today)

        storefilepath = pjoin(rawdir,savename)

        command = 'wget -O {} {}'.format(storefilepath,url)

        os.popen(command)

        # gunzip file
        gunzip = 'gunzip {}'.format(storefilepath)

        os.popen(gunzip)

        return storefilepath.rsplit('.gz',1)[0].strip()

    def tsv(self,filepath):

        tsvfile = open(filepath)

        n = 0

        for line in tsvfile:

            if n == 0:

                keys =[key.strip() for key in  line.split('\t') if key]

            else:

                data =[i for i in line.split('\t') if i]

                dic = dict([(key,val ) for key,val in zip(keys,data)])

                originalSource = dic.get('originalSource').strip()

                if originalSource not in ['CTD_human','HPO','ORPHANET','PSYGENET','UNIPROT']:

                    n += 1
                    continue

                geneId = dic.pop('geneId').strip()

                geneSymbol = dic.pop('geneSymbol').strip()

                diseaseName = dic.pop('diseaseName').strip()

                self.col.update(
                    {'geneId':geneId},
                    {'$set':{'geneSymbol':geneSymbol},
                      '$push':{'{}'.format(diseaseName):dic}
                      },
                      True)

                print 'disgenet line',n,geneId

            n += 1

def main():

    modelhelp = model_help.replace('&'*6,'DisGeNET').replace('#'*6,'disgenet')

    funcs = (downloadData,extractData,updateData,selectData,dbMap,disgenet_disease_store)

    getOpts(modelhelp,funcs=funcs)
        
if __name__ == '__main__':
    main()
   
