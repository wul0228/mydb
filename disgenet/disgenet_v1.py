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

(disgenet_model,disgenet_raw,disgenet_store,disgenet_db,disgenet_map) = buildSubDir('disgenet')

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

        (choice,existgoFile) = lookforExisted(disgenet_raw,'all_gene_disease')

        if choice != 'y':
            return

    if redownload or not existgoFile or  choice == 'y':

        process = disgenet_parser(today)
    
        mt = process.getMt()   

        filepath = process.wget(disgenet_download_url,mt,disgenet_raw)
    
    log_path = pjoin(disgenet_model,'disgenet.log')

    if not os.path.exists(pjoin(disgenet_model,'disgenet.log')):

        with open('./disgenet.log','w') as wf:
            json.dump({'disgenet':[(mt,today,model_name)]},wf,indent=8)

    print  'datadowload completed !'

    return (filepath,today)

def extractData(filepath,version):
    
    process = disgenet_parser(version)

    process.tsv(filepath)
    
    print 'extract and insert completed'

def updateData():

    disgenet_log = json.load(open('./disgenet.log'))

    process = disgenet_parser(today)

    mt = process.getMt()

    if mt != disgenet_log['disgenet'][-1][0]:

        filepath,version = downloadData(redownload=True)

        extractData(filepath,version)

        disgenet_log['disgenet'].append((mt,today,model_name))

        print  '{} \'s new edition is {} '.format('disgenet',mt)
        # create new log
        with open('./disgenet.log','w') as wf:

            json.dump(disgenet_log,wf,indent=2)

    else:
        print  '{} is the latest !'.format('disgenet')

def selectData():

    #function introduction
    #args:
    
    return

class dbMap(object):

    #class introduction

    def __init__(self):
        pass


    def mapXX2XX(self):
        pass

    def mapping(self):

        self.mapXX2XX()

class disgenet_parser(object):

    def __init__(self, version):

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mygene')

        col = db.get_collection('disgenet_{}'.format(version))

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

    funcs = (downloadData,extractData,updateData,selectData,dbMap,disgenet_store)

    getOpts(modelhelp,funcs=funcs)
        
if __name__ == '__main__':
    main()
  
