#!/usr/bin/env python
# --coding:utf-8--
# date: 20171129
# author:wuling
# emai:ling.wu@myhealthgene.com

#this model set  to download,extract,standard insert and select gene data from go

import sys
sys.path.append('../')
sys.setdefaultencoding = ('utf-8')
from share import *
from config import *  

__all__ = ['downloadData','extractData','standarData','insertData','updateData','selectData']

version  = 1.0

model_name = psplit(os.path.abspath(__file__))[1]

(go_gene_model,go_gene_raw,go_gene_store,go_gene_db,go_gene_map) = buildSubDir('go_gene')


log_path = pjoin(go_gene_model,'go_gene.log')

# main code

def downloadOne(go_gene_ftp_infos,filename,rawdir):
    '''
    this function is to download  one file under  a given remote dir 
    args:
    ftp -- a ftp cursor for a specified
    filename --  the name of file need download
    rawdir -- the directory to save download file
    '''
    while  True:

        try:

            ftp = connectFTP(**go_gene_ftp_infos)

            mt =  ftp.sendcmd('MDTM {}'.format(filename)).replace(' ','')

            if not filename.endswith('.gz'):

                savefilename = '{}_{}_{}'.format(filename,mt,today)

            else:
                savefilename = '{}_{}_{}.gz'.format(filename.rsplit('.',1)[0].strip(),mt,today)

            remoteabsfilepath = pjoin(go_gene_ftp_infos['logdir'],'{}'.format(filename))

            save_file_path = ftpDownload(ftp,filename,savefilename,rawdir,remoteabsfilepath)

            print filename,'done'

            return (save_file_path,mt)

        except:

            ftp = connectFTP(**go_gene_ftp_infos)

def downloadData(redownload = False):
    '''
    this function is to download the raw data from go gene FTP WebSite
    args:
    redownload-- default False, check to see if exists an old edition before download
                       -- if set to true, download directly with no check
    rawdir-- the directory to save download file
    '''
    if  not redownload:

        (choice,existgoFile) = lookforExisted(go_gene_raw,'gene')

        if choice != 'y':
            return

    if redownload or not existgoFile or  choice == 'y':

        rawdir = pjoin(go_gene_raw,'gene_{}'.format(today))

        createDir(rawdir)

        func = lambda x:downloadOne(go_gene_ftp_infos,x,rawdir)

        # dodownload gpa and gpi file
        multiProcess(func,go_gene_filenames,size=16)

        # download obo file
        func = lambda x:downloadOne(go_obo_ftp_infos,x,rawdir)

        multiProcess(func,go_obo_filenames,size=16)

        # wget = 'wget -P {}  {}'.format(rawdir,go_gene_obo_link)
        # os.popen(wget)
        # os.rename(pjoin(rawdir,'go.obo'),pjoin(rawdir,'go_human.obo_21320{}_{}'.format(today,today)))

    #gunzip all file
    gunzip = 'gunzip {}/*'.format(rawdir)

    os.popen(gunzip)

    if not os.path.exists(log_path):

        initLogFile('go_gene',model_name,go_gene_model,rawdir=rawdir)

    update_file_heads =dict()

    for filename in listdir(rawdir):

        head = filename.split('_213')[0].strip()

        update_file_heads[head] = filename


    with open(pjoin(go_gene_db,'gene_{}.files'.format(today)),'w') as wf:
        json.dump(update_file_heads,wf,indent=2)

    print  'datadowload completed !'

    filepaths = [pjoin(rawdir,filename) for filename in rawdir]

    return (filepaths,today)

def extractData(filepaths,version):

    for filepath in filepaths:

        filename = psplit(filepath)[1].strip()

        process = gene_parser(filepath,version)

        if filename.count('gpa'):

            process.gpa_main()

        elif filename.count('obo'):

            process.obo()
            
        else:
            pass

        print filename,'done'

    print 'extract and insert complete '

    return (filepaths,version)

def updateData(insert=False):

    go_gene_log = json.load(open(log_path))

    rawdir = pjoin(go_gene_raw,'gene_update_{}'.format(today))

    new = False

    filenames = go_gene_filenames + go_obo_filenames

    for filename in filenames:

        if filename.count('gpa'):

            ftp = connectFTP(**go_gene_ftp_infos)
            ftp_infos = go_gene_ftp_infos

        elif filename.count('obo'):

            ftp = connectFTP(**go_obo_ftp_infos)
            ftp_infos = go_obo_ftp_infos

        mt =  ftp.sendcmd('MDTM {}'.format(filename))

        if mt != go_gene_log.get(filename)[-1][0]:

            new = True

            createDir(rawdir)

            downloadOne(ftp_infos,filename,rawdir)

            go_gene_log[filename].append((mt,today,model_name))

            print  '{} \'s new edition is {} '.format(filename,mt)

        else:
            print  '{} {} is the latest !'.format(filename,mt)

    if new:

        with open(log_path,'w') as wf:

            json.dump(go_gene_log,wf,indent=2)

        (latest_file,version) = createNewVersion(go_gene_raw,go_gene_db,rawdir,'gene_',today)

        if insert:

            insertUpdatedData(go_gene_raw,latest_file,'gene_',version,extractData)

            bakeupCol('go_gene_{}'.format(version),'go_gene')

def selectData(querykey = 'id',queryvalue='GO:0005765'):
    '''
    this function is set to select data from mongodb
    args:
    querykey -- a specified field in database
    queryvalue -- a specified value for a specified field in database
    '''
    conn = MongoClient('127.0.0.1', 27017 )

    db = conn.mydb

    colnamehead = 'go'

    col_names =[col_name for col_name in  db.collection_names() if col_name.startswith(colnamehead)]

    col_names.sort(key = lambda x:x.split('_')[1].strip())
    
    print  '*'*80
    print 'existed collections','\n'

    for index,ver in enumerate(col_names):

        print  'index {}  edition {} '.format(index,ver)
        print ''

    edition = raw_input('chose edition index or enter to latest : ')

    if edition == '':
        col_name = col_names[-1]
    else:
        # col_name = col_names[-1]
        col_name = col_names[int(edition)]

    col = db.get_collection(col_name)

    print '*'*80

    while True:

        queryvalue = str(raw_input('input %s  (q to quit) : ' %  querykey))
        
        if queryvalue == 'q' or queryvalue =='Q':

            break

        else:
            # select with go id
            if querykey == 'id':

                # get go_id basic info
                basic_docs = col.find({querykey:queryvalue})
               
               # gene annotation info
                annotation_docs = col.find({'GO.{}'.format(queryvalue):{'$exists':'true'}})
                
                print 'annotation_info:'
                # out put anotions info
                for doc in annotation_docs:

                    print doc.get('DB_Object_ID')
                    
                    annos = doc['GO'][queryvalue]

                    for a in annos:

                        print a ,'\n'

                    print '-'*50

               # out put basic info
                print '~'*100

                print 'basic_info:','\n'

                for doc in basic_docs:

                    for key,val in doc.items():

                        if key in ['_id',]:
                            continue

                        print key,':',val,'\n'

            elif querykey == 'gene_id':

                annotation_docs = col.find({'DB_Object_ID':queryvalue})

                for doc in annotation_docs:

                    gos = doc.get('GO')

                    for go_id,annos in gos.items():

                        print '~'*100

                        go_basic = col.find_one({'id':go_id})

                        # output go_basic
                        print 'basic_info:','\n'

                        for key,val in go_basic.items():
                            if key in ['_id',]:
                                continue
                            print key,' : ',val
                            print 

                        print '-'*50

                        print 'annotation_info:','\n'

                        # output annotations
                        for a in annos:
                            print a
                            print 

class dbMap(object):

    #class introduction

    def __init__(self,version):

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        col = db.get_collection('go_gene_{}'.format(version))

        self.col = col

        def mapXX2XX():
            pass

    def mapping(self):

        self.mapXX2XX()

class gene_parser(object):
    """docstring for gene_parser"""
    def __init__(self,filepath,version):

        self.version = version

        file = open(filepath)

        self.file = file

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        col = db.get_collection('go_gene_{}'.format(version))

        self.col = col

    def gpa(self):

            n = 0  #440715

            keys = [
            'DB','DB_Object_ID','Qualifier','GO ID','DB:Reference',
            'ECO evidence code','With/From','Interacting taxon ID','Date',
            'Assigned_by','Annotation Extension','Annotation Properties']

            for line in self.file:

                if line.startswith('!'):
                    continue

                line = line.strip().split('\t')

                data = dict([key,val] for key,val in zip(keys,line))

                go_id = data.get('GO ID')

                DB_Object_ID = data.get('DB_Object_ID')

                # for key in ['GO ID','DB_Object_ID']:
                #     data.pop(key)
                
                for key,val in data.items():
                    if key in ['GO ID','DB_Object_ID'] or not val:
                        data.pop(key)

                self.col.update(
                    {'id':go_id},
                    {'$push':{'DB_Object_ID.{}'.format(DB_Object_ID):data}},
                    )

                n += 1

                print 'go.gpa line',n,go_id

                '''
                go.gpa line 315309
                pymongo.errors.WriteError: Resulting document after update is larger than 16777216

                '''
    def gpa_main(self):

            n = 0  #440715

            keys = [
            'DB','DB_Object_ID','Qualifier','GO ID','DB:Reference',
            'ECO evidence code','With/From','Interacting taxon ID','Date',
            'Assigned_by','Annotation Extension','Annotation Properties']

            for line in self.file:

                if line.startswith('!'):
                    continue

                line = line.strip().split('\t')

                data = dict([key,val] for key,val in zip(keys,line))

                go_id = data.get('GO ID')

                DB = data.get('DB')

                DB_Object_ID = data.get('DB_Object_ID')

                for key in ['GO ID','DB_Object_ID','DB']:
                    data.pop(key)

                self.col.update(
                    {'DB':DB,'DB_Object_ID':DB_Object_ID},
                    {'$push':{'GO.{}'.format(go_id):data}},
                    True
                    )

                n += 1

                print 'go.gpa_main line',n,go_id

    def obo(self):

        # skip  term head 
        n = 1 #633722

        for line in self.file:

            if line.count('[Term]'):

                break

            n += 1
        
        aset = dict()

        for line in self.file:

            if line.count('[Term]') or line.count('[Typedef]'):

                if aset:

                    # go_id = aset.get('id')

                    # aset.pop('id')

                    # self.col.update(
                    #     {'go_id':go_id},
                    #     {'$set':aset},
                    #     True
                    #     )
                    self.col.insert(aset)
                     
                    print 'go.obo line',n

                aset = dict()

                if line.count('[Typedef]'):
                    break

            else:

                line = line.strip()

                if  bool(line):

                    (key,val) = tuple(line.strip().split(':',1))

                    key = key.strip()
                    val = val.strip()

                    if key in ['name','namespace','def','id']:

                        aset[key] = val

                    else:

                        if key not in aset:

                            aset[key] = list()

                        aset[key].append(val)

            n += 1
            
        print 'go_obo completed! '

def main():

    modelhelp = model_help.replace('&'*6,'GO_GENE').replace('#'*6,'go_gene')

    funcs = (downloadData,extractData,updateData,selectData,dbMap,go_gene_store)

    getOpts(modelhelp,funcs=funcs)
        
if __name__ == '__main__':
    main()
    # rawdir = '/home/user/project/dbproject/mydb_v1/go_gene/dataraw/gene_171129140122'

    # filepaths = [pjoin(rawdir,filename) for filename in listdir(rawdir)]
    # extractData(filepaths,'171129140122')
