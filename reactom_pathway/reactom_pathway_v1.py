#!/usr/bin/env python
# --coding:utf-8--
# date: 2017/12/04
# author:wuling
# emai:ling.wu@myhealthgene.com

#this model set  to download,extract,standard insert and select pathway data from reactom pathway

import sys
sys.path.append('../')
sys.setdefaultencoding = ('utf-8')
from share import *
from config import *  

__all__ = ['downloadData','extractData','standarData','insertData','updateData','selectData']

version  = 1.0

model_name = psplit(os.path.abspath(__file__))[1]

(reactom_pathway_model,reactom_pathway_raw,reactom_pathway_store,reactom_pathway_db,reactom_pathway_map) = buildSubDir('reactom_pathway')

# main code

def downloadData(redownload = False,rawdir = None):

    '''
    this function is to download the raw data from reactom web
    args:
    redownload-- default False, check to see if exists an old edition before download
                       -- if set to true, download directly with no check
    rawdir-- the directory to save download file
    '''
    if  not redownload:

        (choice,existreactomFile) = lookforExisted(reactom_pathway_raw,'pathway')

        if choice != 'y':
            return

    if redownload or not existreactomFile or  choice == 'y':

        if not rawdir:

            rawdir = pjoin(reactom_pathway_raw,'pathway_{}'.format(today))

            createDir(rawdir)

        process = reactom_parser(today)

        # get mt for every file,store in reactomfile_mt
        newest_mt = process.getMt()

        pathurl_mt,summary_mt = process.getUrl()

        # download graph file
        func  = lambda x:process.wget(x,pathurl_mt[x],rawdir)

        multiProcess(func,pathurl_mt.keys(),size=50)

        # download summary file
        for url,mt in summary_mt.items():

            process.wget(url,mt,rawdir)

    # create log file
    log_path = pjoin(reactom_pathway_model,'reactom_pathway.log')

    if not os.path.exists(pjoin(reactom_pathway_model,'reactom_pathway.log')):

        with open('./reactom_pathway.log','w') as wf:
            json.dump({'reactom_pathway':[(newest_mt,today,model_name)]},wf,indent=8)

    print  'datadowload completed !'

    filepaths = [pjoin(rawdir,filename) for filename in rawdir]

    return (filepaths,today)

def extractData(filepaths,version):

    rawdirname = psplit(psplit(filepaths[0])[0])[1].strip()

    storedir = pjoin(reactom_pathway_store,rawdirname)

    createDir(storedir)

    process = reactom_parser(version)

    jsonpaths = [file for file in filepaths if file.endswith('.json')]

    sumpaths = [file for file in filepaths if file.endswith('.txt')]

    # deal json file 
    func = lambda x : process.graph(x,storedir)

    multiProcess(func,jsonpaths,size=20)

    # deal summary file
    process.sum(sumpaths[0])

    print 'extract and insert completed'

def updateData():

    reactom_pathway_log = json.load(open('./reactom_pathway.log'))

    rawdir = pjoin(reactom_pathway_raw,'pathway_update_{}'.format(today))

    process = reactom_parser(today)

    newest_mt = process.getMt()

    new = False

    if newest_mt != reactom_pathway_log['reactom_pathway'][-1][0]:

        createDir(rawdir)

        new = True 

        filepaths,version = downloadData(redownload = True)

        extractData(filepaths,version)

        reactom_pathway_log['reactom_pathway'].append((newest_mt,today,model_name))

        print  '{} \'s new edition is {} '.format('reactom_pathway',newest_mt)

    else:

        print  '{} is the latest !'.format('reactom_pathway')

    if new:

        # create new log
        with open('./reactom_pathway.log','w') as wf:

            json.dump(reactom_pathway_log,wf,indent=2)

def selectData(querykey = 'stId',value='R-HSA-76071'):
    '''
    this function is set to select data from mongodb
    args:
    querykey -- a specified field in database
    queryvalue -- a specified value for a specified field in database
    '''
    conn = MongoClient('127.0.0.1', 27017 )

    db = conn.mygene

    colnamehead = 'reactom_pathway'

    dataFromDB(db,colnamehead,querykey,queryvalue=None)

class dbMap(object):

    #class introduction

    def __init__(self,version):

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mygene')

        col = db.get_collection('reactom_pathway_{}'.format(version))

        self.col = col

        docs = col.find({})

        self.docs = docs

        self.version = version

    def mapID2Name(self):

        id_name = dict()

        name_id = dict()

        for doc in self.docs:

            path_id = doc.get('stId')

            path_name = doc.get('name')

            if path_name:

                id_name.update({path_id:path_name})

                if path_name not in name_id:

                    name_id[path_name] = list()

                name_id[path_name].append(path_id)

            else:
                print path_id,'no name'

        with open(pjoin(reactom_pathway_map,'id2name_{}.json'.format(self.version)),'w') as wf:
            json.dump(id_name,wf,indent=2)

        with open(pjoin(reactom_pathway_map,'name2id_{}.json'.format(self.version)),'w') as wf:
            json.dump(name_id,wf,indent=2)

    def mapping(self):

        self.mapID2Name()

class reactom_parser(object):

    def __init__(self,version):

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mygene')

        col = db.get_collection('reactom_pathway_{}'.format(version))

        self.col = col
        
        self.version = version

    def Mt(text,filename):

        mt = text.split(filename)[1].strip().rsplit(' ',1)[0]

        for sym in [' ',':','-','.json','.txt']:

            mt = mt.replace(sym,'')

        mt = '213' + mt

        return mt

    def getMt(self):

            web = requests.get(reactome_download_web1)

            soup = bs(web.content,'lxml')

            trs = soup.select('body > table > tr ')

            for tr in trs:

                text = tr.text

                if text.count('diagram/'):

                     mt = Mt(text,'diagram/')

            return mt

    def getUrl(self):

        '''
        this function is to get all files update time from http://download web page
        '''
        # diagram file
        # for url in [reactome_download_web1,reactome_download_web2]:

        pathurl_mt = dict()

        summary_mt = dict()

        for url in [reactome_download_web1,reactome_download_web2]:
            
            web = requests.get(url)

            soup = bs(web.content,'lxml')

            trs = soup.select('body > table > tr ')

            for tr in trs:

                text = tr.text

                filename = text.split('.json')[0].strip()

                if filename.startswith('R-HSA') and filename.endswith('.graph'):

                    mt = Mt(text,filename)

                    key = url + filename + '.json'

                    pathurl_mt[key] = mt

                elif filename.split('.txt')[0].strip() == 'pathway2summation':

                    mt = Mt(text,'pathway2summation.txt')

                    key = url + 'pathway2summation.txt'

                    summary_mt[key] = mt

        return (pathurl_mt,summary_mt)

    def wget(self,url,mt,rawdir):


        filename = url.rsplit('/',1)[1].strip().replace('.graph.json','').replace('.txt','')

        if filename == 'pathway2summation':

            savename = '{}_{}.txt'.format(filename,mt)

        else:

            savename = '{}_{}.graph.json'.format(filename,mt)

        storefilepath = pjoin(rawdir,savename)

        command = 'wget -O {} {}'.format(storefilepath,url)

        os.popen(command)

    def graph(self,filepath,storedir):

        filename = psplit(filepath)[1].strip()

        dic = json.load(open(filepath))

        nodes = dic.pop('nodes')

        # delete origin nodes list
        
        # create new nodes dict
        nodeinfo = dict()

        for node in nodes:

            schemaClass = node.pop('schemaClass')

            if not schemaClass:

                schemaClass = 'Other'


            dbId = node.pop('dbId')

            if schemaClass not in nodeinfo:

                nodeinfo[schemaClass] = dict()

            if dbId not in nodeinfo[schemaClass]:

                nodeinfo[schemaClass][str(dbId)] = node

        dic.update({'nodes':nodeinfo})

        with open(pjoin(storedir,filename) ,'w') as wf:
            json.dump(dic,wf,indent=2)

        self.col.insert(dic)
        # pprint.pprint(dic)
    
    def sum(self,filepath):

        filename = psplit(filepath)[1].strip()

        tsvfile = open(filepath)

        keys = ['stId','name','summation']

        for line in tsvfile:

            if line.startswith('#'):
                continue

            data =[i.strip() for i in  line.strip().split('\t') if i]

            dic = dict([(key,val) for key,val in zip(keys,data)])

            stId = dic.pop('stId')

            summation = dic['summation'].decode('unicode_escape').encode('utf8')
            dic['summation'] = summation

            self.col.update(
                {'stId':stId},
                {'$set':dic})

def main():

    modelhelp = model_help.replace('&'*6,'Reactom_Pathway').replace('#'*6,'reactom_pathway')

    funcs = (downloadData,extractData,updateData,selectData,dbMap,reactom_pathway_store)

    getOpts(modelhelp,funcs=funcs)
        
if __name__ == '__main__':
    main()

    # filepaths,version = downloadData(redownload = True)

    # extractData(filepaths,version)
