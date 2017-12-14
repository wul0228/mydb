#!/usr/bin/env python
# --coding:utf-8--
# date:20171023
# author:wuling
# emai:ling.wu@myhealthgene.com

#+++++++++++++++++++++++++ packages ++++++++++++++++++++++++++++++++++++++#
import sys
reload(sys)
sys.path.append('..')
sys.setdefaultencoding = ('utf-8')
from config import *

#+++++++++++++++++++++++++ main code ++++++++++++++++++++++++++++++++++++++
def initDB():

    conn = MongoClient('localhost',27017)

    db = conn.get_database('mydb')

    cols = db.collection_names()

    db_cols = dict()

    db_colnames = dict()

    for colname in cols:

        modelname = colname.rsplit('_',1)[0].strip()

        db_cols[modelname] = db.get_collection(colname)

        db_colnames[modelname] = colname

    return (db,db_cols,db_colnames)

def createDir(dirpath):
    '''
    this function is to create directory if not exist
    '''
    if not os.path.exists(dirpath):

        os.mkdir(dirpath)

    return dirpath

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

                    pprint(doc,indent=4,width=80,depth=None)
                    
                    doc.pop('_id')
                    
                    # with open('./out_{}.json'.format(docnum),'w') as wf:
                    #     json.dump(doc,wf,indent=8)

                    print '~'*50
                   
                    n += 1
                
                print 'allfind:',n
       
            else:
                print 'No record'

            print '-'*80

def getOpts(modelhelp,funcs):
    
    (update,select) = funcs
    
    try:

        (opts,args) = getopt.getopt(sys.argv[1:],"hqu:",['--help',"--update","--query"])

        querykey,queryvalue=("","")

        for op,value in opts:

            if op in ("-h","--help"):

                print modelhelp

            elif op in ('-u','--update'):
                update()

            elif op in ('-q','--query'):
                select()
               
    except getopt.GetoptError:

        sys.exit()

def main():

    pass
if __name__ == '__main__' :

    main()