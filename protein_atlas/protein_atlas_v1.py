#!/usr/bin/env python
# --coding:utf-8--
# date: 2017/12/08
# author:wuling
# emai:ling.wu@myhealthgene.com

#this model set  set  to download,extract,standard insert and select gene data from proteinAtlas

import sys
sys.path.append('../')
sys.setdefaultencoding = ('utf-8')
from share import *
from config import *  

__all__ = ['downloadData','extractData','standarData','insertData','updateData','selectData']

version  = 1.0

model_name = psplit(os.path.abspath(__file__))[1]

(protein_atlas_model,protein_atlas_raw,protein_atlas_store,protein_atlas_db,protein_atlas_map) = buildSubDir('protein_atlas')

log_path = pjoin(protein_atlas_model,'protein_atlas.log')

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

        (choice,existgoFile) = lookforExisted(protein_atlas_raw,'all_gene_disease')

        if choice != 'y':
            return

    if redownload or not existgoFile or  choice == 'y':

        process = protein_parser(today)
    
        download_url,mt = process.getMt()   

        filepath = process.wget(download_url,mt,protein_atlas_raw)

    if not os.path.exists(pjoin(protein_atlas_model,'protein_atlas.log')):

        with open(log_path,'w') as wf:
            json.dump({'protein_atlas':[(mt,today,model_name)]},wf,indent=8)

    print  'datadowload completed !'

    return (filepath,today)

def extractData(filepath,version):
    
    process = protein_parser(version)

    process.tsv(filepath)
    
    print 'extract and insert completed'

    return (filepath,version)
    
def updateData():

    protein_atlas_log = json.load(open(log_path))

    process = protein_parser(today)

    (download_url,mt) = process.getMt()

    if mt != protein_atlas_log['protein_atlas'][-1][0]:

        filepath,version = downloadData(redownload=True)

        extractData(filepath,version)

        protein_atlas_log['protein_atlas'].append((mt,today,model_name))

        print  '{} \'s new edition is {} '.format('protein_atlas',mt)
        # create new log
        with open('./protein_atlas.log','w') as wf:

            json.dump(protein_atlas_log,wf,indent=2)

        bakeupCol('protein_atlas_{}'.format(version),'protein_atlas')

    else:
        print  '{} is the latest !'.format('protein_atlas')


def selectData(querykey = 'Ensembl',value='ENSG00000000003'):
    '''
    this function is set to select data from mongodb
    args:
    querykey -- a specified field in database
    queryvalue -- a specified value for a specified field in database
    '''
    conn = MongoClient('127.0.0.1', 27017 )

    db = conn.mydb

    colnamehead = 'protein_atlas_'

    dataFromDB(db,colnamehead,querykey,queryvalue=None)


class dbMap(object):

    #class introduction

    def __init__(self,version):

        self.version = version

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        col = db.get_collection('protein_atlas_{}'.format(self.version))

        self.col = col

    def mapXX2XX(self):
        pass

    def mapping(self):

        self.mapXX2XX()

class protein_parser(object):

    def __init__(self, version):

        self.version = version

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        col = db.get_collection('protein_atlas_{}'.format(self.version))

        self.col = col


    def getMt(self):

        headers = {'User-Agent':'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Ubuntu Chromium/37.0.2062.120 Chrome/37.0.2062.120 Safari/537.36',}

        web = requests.get(protein_atlas_download_web,headers = headers,verify=False)

        soup = bs(web.content,'lxml')

        down = soup.find(text = 'proteinatlas.tsv.zip')

        a = down.findParent('a')

        href = a.attrs.get('href')

        tr = a.findParent('tr')

        mt = tr.text.split('version')[1].strip().split(' ')[0].strip()

        download_url = 'https://www.proteinatlas.org/' + href

        return (download_url,mt)

    def wget(self,url,mt,rawdir):

        filename = url.rsplit('/',1)[1].strip().replace('.tsv.zip','')

        savename = '{}_{}_{}.tsv.zip'.format(filename,mt.replace('.','*'),today)

        storefilepath = pjoin(rawdir,savename)

        command = 'wget -O {} {}'.format(storefilepath,url)

        os.popen(command)

        # unzip file
        unzip = 'unzip -d {}  "{}" '.format(rawdir,storefilepath)

        os.popen(unzip)

        savefilepath = storefilepath.replace('.tsv.zip','.tsv',1)

        os.rename(pjoin(rawdir,'proteinatlas.tsv'),savefilepath)

        os.remove(storefilepath)

        return savefilepath

    def tsv(self,filepath):
        
        tsvfile = open(filepath)

        n = 0

        for line in tsvfile:

            if  n == 0:

                keys = [i.strip().replace(' ','&').replace('.','*') for i in line.strip().split('\t')]

            else:

                data =[i.strip() for i in  line.split('\t')]

                dic = dict([(key,val) for key,val in zip(keys,data)])


                synonym =  [i.strip() for i in dic.pop('Gene&synonym').split(',')]

                protein_class = [i.strip() for i in dic.pop('Protein&class').split(',')]

                Antibody = [i.strip() for i in dic.pop('Antibody').split(',')]

                Prognostic_p_value = [i.strip() for i in dic.pop('Prognostic&p-value').split(',')]

                RNA_TS_TPM  =  [i.strip() for i in dic.pop('RNA&TS&TPM').split(',')]

                RNA_CS_TPM =  [i.strip() for i in dic.pop('RNA&CS&TPM').split(';')]

                dic.update({
                    'Gene&synonym':synonym,
                    'Protein&class':protein_class,
                    'Antibody':Antibody,
                    'Prognostic&p-value':Prognostic_p_value,
                    'RNA&TS&TPM':RNA_TS_TPM,
                    'RNA&CS&TPM':RNA_CS_TPM
                    })

                self.col.insert(dic)

                print 'protein atlas line',n ,dic.get('Gene')

            n += 1

def main():

    modelhelp = model_help.replace('&'*6,'PROTEIN ATLAS').replace('#'*6,'protein atlas')

    funcs = (downloadData,extractData,updateData,selectData,dbMap,protein_atlas_store)

    getOpts(modelhelp,funcs=funcs)
        
if __name__ == '__main__':

    main()
