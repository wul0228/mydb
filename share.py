#!/usr/bin/env python
# --coding:utf-8--
# date:20171023
# author:wuling
# emai:ling.wu@myhealthgene.com

#+++++++++++++++++++++++++ packages ++++++++++++++++++++++++++++++++++++++#

from config import *

reload(sys)
sys.path.append('..')
sys.setdefaultencoding = ('utf-8')


#+++++++++++++++++++++++++ main code ++++++++++++++++++++++++++++++++++++++#

def createDir(dirpath):
    '''
    this function is to create directory if not exist
    '''
    if not os.path.exists(dirpath):

        os.mkdir(dirpath)

    return dirpath

def buildSubDir(name):
    '''
    this function is to build all sub-directory for specified
    '''
    db = pjoin(mydb_path,name)

    # _load = pjoin(mymol_path,name,'dataload')

    _raw = pjoin(mydb_path,name,'dataraw')

    _store =  pjoin(mydb_path,name,'datastore')

    _db = pjoin(mydb_path,name,'database')

    _map = pjoin(mydb_path,name,'datamap')


    for _dir in [ db,_raw,_store,_db,_map]:

        createDir(_dir)      

    return (db,_raw,_store,_db,_map)

def initLogFile(parser,modelname,storedir,mt=None,rawdir=None):

    with open(pjoin(storedir,'{}.log'.format(parser)),'w') as wf: 

        log_dict = dict()

        if any([parser.startswith(start) for start in [ 'ncbi','go','reactom','wiki'] ]):

            for filename in listdir(rawdir):

                if not filename.endswith('.gz'):
                    
                    name = filename.split('_213')[0].strip()

                else:
                    name = filename.split('_213')[0].strip() + '.gz'

                mt = '213 ' + filename.split('_213')[1].strip().split('_',1)[0].strip()

                if name not in log_dict:

                 log_dict[name] = list()

                log_dict[name].append((mt,today,modelname))

            json.dump(log_dict,wf,indent=2)

        elif parser.startswith('ensembl'):

            for filename in listdir(rawdir):

                mt = '213 ' + filename.split('_213')[1].strip().split('_',1)[0].strip()

                ver = ['GRCh37','GRCh38']

                for v in ver:

                    if filename.count(v):

                        for key,mark in ensembl_file_mark.items():

                            if key.count(v):

                                if filename.count(mark):

                                    log_dict[key] = list()

                                    log_dict[key].append((mt,today,modelname))

            json.dump(log_dict,wf,indent=2)

        elif parser.startswith('kegg'):

            for filename in listdir(rawdir):

                mt =filename.split('_',1)[1].strip().split('_',1)[0].strip()

                name =filename.split('_',1)[0].strip() + '.json'

                if name not in log_dict:

                    log_dict[name] = list()

                log_dict[name].append((mt,today,modelname))

            json.dump(log_dict,wf,indent=2)
    
def lookforExisted(datadir,dirnamehead):

        existFile = filter(lambda x:x.startswith(dirnamehead),listdir(datadir))

        tips = '''
        there have been stored  editions of  below: \n
        {} \n
        if you still want to download again? 
        chose  y/n :'''.format(existFile)
     
        choice  = str(raw_input(tips))

        if choice != 'y':
            return  (None,None)

        else:
            return (choice,existFile)

def choseExisted(existFile):

    tips =  '''
    there are {}  edition below, please chose one of them to continue ?
    {}
    input a index like 0,1,2... (input 'q' to quit):'''.format(len(existFile), \
    [ "{} {} ".format(str(index),name) for index,name in enumerate(existFile) ])

    while True:

        index = raw_input(tips)

        if str(index) == 'q':
            return 'q'
        try:
            edition = int(index)
        except  Exception,e:
            print e

        if edition in range(len(existPubChemFile)):
            return edition

        else:
            print '\n !!! index out of range.please check again \n'

def connectFTP(host,user=None,passwd=None,logdir=None):
    '''
    this function is to connect  ftp site 
    '''
    ftp = FTP(host)

    if user or passwd:

        ftp.login(user,passwd)

    else:

        ftp.login()

    ftp.cwd(logdir)

    return ftp

def  ftpDownload(ftp,filename,savefilename,rawdir,remoteabsfilepath):
    '''
    this function is to download specified file from ftp site
    '''
    bufsize=1024 

    save_file_path =pjoin(rawdir,savefilename)

    file_handle=open(save_file_path,'wb')

    ftp.retrbinary('RETR {}'.format(remoteabsfilepath) ,file_handle.write ,bufsize) 
   
    ftp.set_debuglevel(0) 

    return save_file_path


def connect2DB(server = 'mongodb',host='localhots',port=27017,dbname='ChEBI'):
    '''
    this function is set to connect  database mymol in localhost mysql server
    '''
    if server == 'mongodb':

        db = MongoClient('mongodb://{}:{}/{}'.format(host,port,dbname))

        return db

    elif server == 'mysql':

        connection = MySQLdb.connect(host=host,port=port,user='root',passwd='281625',db=dbname)

        cursor = connection.cursor()

        return cursor

    else:

        print  'no server input'

def atomPair(): 

    #构建需要replace的带电原子类型与其对应的中性原子的pair对
    patts= (
    # Imidazoles
    ('[n+;H]','n'),
    # Amines
    ('[N+;!H0]','N'),
    # Carboxylic acids and alcohols ('[$([O-]);!$([O-][#7])]','O'), # Thiols
    ('[S-;X1]','S'),
    # Sulfonamides 
    ('[$([N-;X2]S(=O)=O)]','N'), 
    # Enamines 
    ('[$([N-;X2][C,N]=C)]','N'), 
    # Tetrazoles 
    ('[n-]','[nH]'),
    # Sulfoxides 
    ('[$([S-]=O)]','S'),
    # Amides
    ('[$([N-]C=O)]','N'),
    #
    ('[O-;X1]',"O"),
    #
    ('[O+;X3]',"O"),
    #
    ('[$([O-]=C)]','O'),
    #20170308
    ('[C-;X3]','C'),
    #20170308
    # ('[Se-;H2]','Se'),
    ('[c-;X3]','c')
    )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

def neutrCharge(smiles):
    """
    this function is to transform the charged smiles to electroneutral
    """
    try:
        mol = mfsmi(smiles)
    except:
        return smiles

    if not mol:
        return smiles
 
    atomPairs= atomPair()

    #circulate all pairs to find if there is a substructure in pre-defined pairs
    for i,(reactant, product) in enumerate( atomPairs):

        while mol.HasSubstructMatch(reactant): 

            #ReplaceSubstructures可选择Replacement = True（默认为False）一步替换所有
            rms = AllChem.ReplaceSubstructs(mol, reactant, product) 

            #rms是一个tuple,内含多个重复的mol，原因不明
            mol = rms[0] 

    # ****** C[O+](C)C
    # ****** [HH].CC1CCC2C(C(OC3C24C1CCC(O3)(O4)C)OC5C(C6CCC(C7C68C(O5)(OC(O8)(CC7)C)C)C)C)C
    try:
        return mtsmi(mol)
    except:
        print '******',smiles
        return smiles


def dataFromDB(database,colnamehead,querykey,queryvalue=None):

    '''
    this function is set to select data from mongodb
    '''
    # get all edition collection
    col_names =[col_name for col_name in  database.collection_names() if col_name.startswith(colnamehead)]

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

    col = database.get_collection(col_name)

    print '*'*80

    while True:

        queryvalue = str(raw_input('input %s  (q to quit) : ' %  querykey))
        
        if queryvalue == 'q' or queryvalue =='Q':

            break

        else:
            
            docs = col.find({querykey:queryvalue})
           
            n = 0

            if docs:

                print '\n','Result: ','\n'

                docnum = 0

                for doc in docs:
                    docnum += 1

                    pprint.pprint(doc,indent=4,width=80,depth=None)
                    
                    doc.pop('_id')
                    
                    with open('./out_{}.json'.format(docnum),'w') as wf:
                        json.dump(doc,wf,indent=8)

                    print '~'*50
                   
                    n += 1
                
                print 'allfind:',n
       
            else:
                print 'No record'

            print '-'*80

def bakeupCol(colname,colhead):

    conn = MongoClient('localhost',27017)

    db = conn.get_database('mydb')

    # bakeup new version to _mongodb directory
    bakeup =  '/usr/local/mongodb/bin/mongodump -d mydb -c {}  -o  ../_mongodb/'.format(colname)

    os.popen(bakeup)

    print '{} bakeup completed'.format(colname)

    # remove old version
    col_names =[col_name for col_name in db.collection_names() if col_name.startswith(colhead)]

    for name in col_names:

        if name != colname:

            col = db.get_collection(name)

            col.drop()
            

def value2key(dic):
    '''
    this function is to overture a value as the key and the key as vaule, just like {a:[1,2],b:[3,2]} return {1:a,2:[a,b],3:b}
    note that the value of key must be list and the element of list must be hashable
    '''
    value_key = dict()

    for key,val in dic.items():

        if isinstance(val,unicode):

            if val not in value_key:
                value_key[val] = list()
            value_key[val].append(key)

        elif isinstance(val,list):

            for v in val :
                if v not in value_key:
                    value_key[v] = list()
                value_key[v].append(key)

    return value_key

def deBlankDict(dic):
    '''
    this function is to 
    a. delete key from dic if val is None
    b. if val is list but only contain  a element  so list transfer to this element
    c. if val is list  and have multi elements ,first dedup ,an then  if dict included , iterate to deblank
    d. if val is dict but only have one key , so ,delete the key from parent-dict and update with sub-dict
    '''
    
    for key,val in dic.items():

        if not val:
            
            dic.pop(key)

        elif isinstance(val,list):

            if len(val) == 1:

                dic[key] =val[0]

            else:
                # dedup
                _val = list()

                for v in val:

                    if v not in _val:

                        if isinstance(v,dict):

                            v = deBlankDict(v)

                        _val.append(v)

                dic[key] = _val

        elif isinstance(val,dict):

            if len(val.keys()) ==1:

                val_key = val.keys()[0]

                val_val = val[val_key]

                dic.pop(key)

                dic[val_key] = val_val

            else:
                deBlankDict(val)

    return dic


def strAndList(content):
    '''
    this function is to combine the unicode and list containing multi unicode like ['a',['b','c']] ,return a ['a','b','c']
    '''
    all_unicode = list()
    
    for i in content:

        if isinstance(i,unicode):

            all_unicode.append(i)

        elif isinstance(i,list):

            all_unicode += i

    return list(set(all_unicode))

def strAndDict(content,key):
    '''
    this function is to combine the unicode and dict value containing multi unicode like ['a',{'b':'c'}] ,return a and the valule of specified key like 'b', ['a','c']
    '''

    if isinstance(content,dict):

        return [content[key],]

    elif isinstance(content,list):

        tmpfunc = lambda x:x[key] if isinstance(x,dict) else x

        return [tmpfunc(i) for i in content]

def retuOneDictValue(content,key):

    if isinstance(content,dict):

        return content.get(key)

    elif isinstance(content,list):

       for i in content:
            if isinstance(i,dict):
                 return i.get(key)

def multiProcess(func,args,size=16):
    '''
    this function is to concurrent processing
    size -- run  processes all at once
    func --- the  function to run
    args  -- argument for function
    '''

    pool = ThreadPool(size)

    results = pool.map(func,args)

    pool.close

    pool.join

def createNewVersion(rawdir,dbdir,updated_storedir,latest_file_head,version):
    '''
    this function is to create a new .files in ****/database/  to record the newest  version
    args:
    updated_storedir -- the directory  store  updated data
    '''
    update_file_heads = dict()

    for filename in listdir(updated_storedir):

        head = filename.split('_213')[0].strip()

        update_file_heads[head] = filename

        # get the latest .files file that contain the latest files name in mongodb 
    filenames = [name for name in listdir(dbdir) if name.endswith('.files')]

    print 'filenames',filenames

    latest = sorted(filenames,key=lambda x:x.split(latest_file_head)[1].strip())[-1]
        
    latest_file = json.load(open(pjoin(dbdir,latest)))

        # update the latest  files   with  updated  file
    for head ,name in update_file_heads.items():

        latest_file[head] = name

    with open(pjoin(dbdir,'{}{}.files').format(latest_file_head,version),'w') as wf:

        json.dump(latest_file,wf,indent=2)

    return (latest_file,version)

def insertUpdatedData(rawdir,latest_file,latest_file_head,version,extractData): 
    '''
    this function is to generate a  new collection in mongodb PubChem  with updated date and the old
    args:
    latest_file -- the file record the latest file names in newest version
    '''
    # update mongodb ,create new edition

        # all file name in new vrsion
    insertFiles = latest_file.values()

    _rawdirs =[dir_name for dir_name in listdir(rawdir) if dir_name.startswith(latest_file_head)]

    num = 0

    filepaths = list()
    # circulate to insert all files in insertFiles
    for _dir  in _rawdirs:

        dir_path = pjoin(rawdir,_dir)

        for filename in listdir(dir_path):

            if filename in insertFiles:

                filepath = pjoin(dir_path,filename)

                filepaths.append(filepath)

    print filepaths

    extractData(filepaths,version)

    print 'insertUpdated completed'

def getOpts(modelhelp,funcs,insert=True):
    
    (downloadData,extractData,updateData,selectData,dbMap,model_store) = funcs
    try:

        (opts,args) = getopt.getopt(sys.argv[1:],"haumf:",['--help',"--all","--update","--map","--field=="])

        querykey,queryvalue=("","")

        for op,value in opts:

            if op in ("-h","--help"):

                print modelhelp

            elif op in ('-a','--all'):

                save,version = downloadData(redownload=True)
                store,version = extractData(save,version)

                _map = dbMap(version)
                _map.mapping()

            elif op in ('-u','--update'):
                updateData()

            elif op in ('-f','--field'):
                selectData(value)
               
    except getopt.GetoptError:

        sys.exit()
      
class DateEncoder(json.JSONEncoder):  

    def default(self, obj):  

        if isinstance(obj, datetime):  

            return obj.strftime('%Y-%m-%d %H:%M:%S')  

        # elif isinstance(obj, date):  

        #     return obj.strftime("%Y-%m-%d")  

        else:  

            return json.JSONEncoder.default(self, obj) 

def main():

    print neutrCharge('c1ccccc1CC([NH3+])C(=O)[O-]')

if __name__ == '__main__' :

    main()