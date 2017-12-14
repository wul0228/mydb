#!/usr/bin/env python
# ---coding:utf-8---
# date:20171123
# author:wuling
# emai:ling.wu@myhealthgene.com

# this model is set to
__all__ = ['manage_help','common_help','topic_help','py_template','topic_intros']

manage_help = '''

Usage: python manage.py  [OPTION]...[MODELNAME]...

Manage topic and sub-topic 

options:

-h, --help                :give this help

-i, --init    [topicname] :init a new topic

-u,--update  [topicname] :update a topic in current directory,if topicname=all,update all

-t, --topicbase [topicname] : query data record from this topic

-v, --value               : look for all database with filed = value

-o,-- output              : the  directory  path to store query result 
    eg:
        python manage.py  -t  gene  -v  TP53  -o  _result/
'''

common_help = '''
'''
topic_help ='''
Usage: python manage.py  [OPTION]...[NAME]...

-h, --help               :give this help
-u, --update             :update map and build a new edition topic base of ######
-q, --query              :look for database with default field
'''
py_template = '''

#!/usr/bin/env python
# --coding:utf-8--
# date:2017/12/13
# author:wuling
# emai:ling.wu@myhealthgene.com

# this topic is set  to construct a topic base

import sys
reload(sys)
sys.path.append('../')
sys.setdefaultencoding = ('utf-8')
from config import *
from share import *

current_path = psplit(os.path.abspath(__file__))[0]

class dbMap(object):

    def __init__(self, ):

        super(dbMap, self).__init__()
    
        (db,db_cols,db_colnames) = initDB()

        self.db = db

        self.db_cols = db_cols

        self.db_colnames = db_colnames

        map_path = pjoin(current_path,'map_{}'.format(today))

        if not os.path.exists(map_path):
            os.mkdir(map_path)

        self.map_path = map_path

    def xxxx2xxxx(self):
        pass

    def mapping(self):
        xxxx2xxxx = self.xxxx2xxxx()
        return (savefilepath,savemappath)

class dbCreate(object):

    def __init__(self,filepath,maplogpath):

        super(dbCreate, self).__init__()

        version = psplit(psplit(filepath)[0])[1].split('map_')[1]

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mytopic')

        col = db.get_collection('gene_{}'.format(version))

        self.col = col

        self.conn = conn

        self.version = version

        self.filepath = filepath

        self.maplogpath = maplogpath

    def build(self):
        pass
def update():

    man = dbMap()
    (filepath,maplogpath) = man.mapping()

    man = dbCreate(filepath,maplogpath)
    man.build()

def select(queryvalue='SNRPD3'):
    
    #this function is set to select data from mongodb
    #args:
    #querykey -- a specified field in database
    #queryvalue -- a specified value for a specified field in database
    
    conn = MongoClient('127.0.0.1', 27017 )

    db = conn.mytopic

    colnamehead = 'gene_'

    dataFromDB(db,colnamehead,querykey='symbol',queryvalue=None)

def main():

    topichelp = topic_help.replace('#'*6,'gene')

    funcs = (update,select)
    
    getOpts(topichelp,funcs=funcs)

if __name__  == '__main__':
    main()
'''

topic_intros = '''
++++++ mygene_v1 Documentation ++++++

edited@xxxxxx

please direct all questions to author@xxxxxx.com

1. brief introduction of sub-files

2. description about 


Functions
 
       
Design  

Usage: python ######.py  [OPTION]

-h, --help               :give this help
-u, --update             :update map and build a new edition topic base of ######
-q, --query              :look for database with default field
   

++++++ mygene_v1  Documentation ++++++
'''