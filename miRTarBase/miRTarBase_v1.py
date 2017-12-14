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

(miRTarBase_model,miRTarBase_raw,miRTarBase_store,miRTarBase_db,miRTarBase_map) = buildSubDir('miRTarBase')

log_path = pjoin(miRTarBase_model,'miRTarBase.log')

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

        (choice,existgoFile) = lookforExisted(miRTarBase_raw,'miRTarBase')

        if choice != 'y':
            return

    if redownload or not existgoFile or  choice == 'y':

        process = mirtarbase_parser(today)

        download_url,mt = process.getMt()

        filepath = process.wget(download_url,mt,miRTarBase_raw)

    if not os.path.exists(log_path):

        with open(log_path,'w') as wf:
            json.dump({'miRTarBase':[(mt,today,model_name)]},wf,indent=8)

    print  'datadowload completed !'

    return (filepath,today)

def extractData(filepath,version):

    process = mirtarbase_parser(version)

    process.xlsx(filepath)

    print 'extract and insert completed'

    return (filepath,version)

def updateData(insert=False,_mongodb='../_mongodb/'):

    miRTarBase_log = json.load(open(log_path))

    process = mirtarbase_parser(today)

    (download_url,mt) = process.getMt()

    if mt != miRTarBase_log['miRTarBase'][-1][0]:

        filepath,version = downloadData(redownload=True)

        extractData(filepath,version)

        miRTarBase_log['miRTarBase'].append((mt,today,model_name))

        print  '{} \'s new edition is {} '.format('miRTarBase',mt)
        # create new log
        with open('./miRTarBase.log','w') as wf:

            json.dump(miRTarBase_log,wf,indent=2)

        bakeupCol('miRTarBase_{}'.format(version),'miRTarBase',_mongodb)

    else:
        print  '{} is the latest !'.format('miRTarBase')

def selectData(querykey = 'miRTarBase&ID',value='MIRT042081'):
    '''
    this function is set to select data from mongodb
    args:
    querykey -- a specified field in database
    queryvalue -- a specified value for a specified field in database
    '''
    conn = MongoClient('127.0.0.1', 27017 )

    db = conn.mydb

    colnamehead = 'miRTarBase'

    dataFromDB(db,colnamehead,querykey,queryvalue=None)

class dbMap(object):

    #class introduction

    def __init__(self,version):

        self.version = version

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        colname = 'miRTarBase_{}'.format(self.version)

        col = db.get_collection(colname)

        self.col = col

        self.colname = colname

    def mapgenId2mirID(self):
        
        docs = self.col.find({})

        mirID_geneID = dict()

        for doc in docs:

            mirID = doc.get('miRTarBase&ID')

            geneID = doc.get('gene').keys() 

            if mirID not in mirID_geneID:

                mirID_geneID[mirID] = list()

            mirID_geneID[mirID] += geneID

        geneID_mirID = value2key(mirID_geneID)

        map_dir = pjoin(miRTarBase_map,self.colname)

        createDir(map_dir)

        with open(pjoin(map_dir,'geneID2mirID.json'),'w') as wf:
            json.dump(geneID_mirID,wf,indent=8)

        with open(pjoin(map_dir,'mirID2geneID.json'),'w') as wf:
            json.dump(mirID_geneID,wf,indent=8)

    def mapping(self):

        self.mapgenId2mirID()

class mirtarbase_parser(object):

    def __init__(self, version):

        # super(mirtarbase_parser, self).__init__()

        self.version = version

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        col = db.get_collection('miRTarBase_{}'.format(self.version))

        self.col = col

    def getMt(self):

        headers = {'User-Agent':'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Ubuntu Chromium/37.0.2062.120 Chrome/37.0.2062.120 Safari/537.36',}

        web = requests.get(miRTarbase_download_web,headers = headers,verify=False)

        soup = bs(web.content,'lxml')

        down = soup.find(text = 'miRTarBase Download')

        a =  down.find_next(name='a')

        href  =a.attrs.get('href')

        mt  = href.split('download/')[1].split('/')[0].strip()

        download_url = miRTarbase_homepage + href.split('/',1)[1].strip()

        return (download_url,mt) 

    def wget(self,url,mt,rawdir):

        filename = url.rsplit('/',1)[1].strip().replace('.xlsx','')

        savename = '{}_{}_{}.xlsx'.format(filename,mt.replace('.','*'),today)

        storefilepath = pjoin(rawdir,savename)

        command = 'wget -O {} {}'.format(storefilepath,url)

        os.popen(command)

        return storefilepath

    def xlsx(self,filepath):

        excel = xlrd.open_workbook(filepath)

        for booksheet in excel.sheets():

            keys = list()

            for row in xrange(booksheet.nrows):

                # create keys list
                
                data = list()

                for col in xrange(booksheet.ncols):

                    cel = booksheet.cell(row, col)

                    val = cel.value

                    if row == 0:

                        keys.append(val)

                    else:

                        data.append(val)

                    # print row,col,val

                if row != 0:

                    dic = dict([(key,val) for key,val in zip(
                        keys,data)])

                    mirtarbase_id = dic.pop('miRTarBase ID')

                    mirtarbase = dic.pop('miRNA')

                    target_gene_id = str(int(dic.pop('Target Gene (Entrez ID)')))

                    target_gene = dic.pop('Target Gene')

                    for key in ['Species (miRNA)','Species (Target Gene)']:
                        dic.pop(key)

                    self.col.update(
                        {'miRTarBase&ID':mirtarbase_id},
                        {'$set':{'miRNAs':mirtarbase,'gene.{}.Target&Gene'.format(target_gene_id):target_gene},
                         '$push':{'gene.{}.exprement&refs'.format(target_gene_id):dic}
                        },
                        True)

                    print 'miRTarBase line',row,'target_gene_id',target_gene_id,len(dic)
        
def main():

    modelhelp = model_help.replace('&'*6,'miRTarBase').replace('#'*6,'miRTarBase')

    funcs = (downloadData,extractData,updateData,selectData,dbMap,miRTarBase_store)

    getOpts(modelhelp,funcs=funcs)
        
if __name__ == '__main__':
    main()
     # downloadData(redownload = True)
     # filepath = '/home/user/project/dbproject/mydb_v1/miRTarBase/dataraw/miRTarBase_MTI_7*0_171211094532.xlsx'
     # version = '171211094532'
     # extractData(filepath,version)

    # man = dbMap('171211094532')
    # man.mapping()