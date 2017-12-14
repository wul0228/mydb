#!/usr/bin/env python
# ---coding:utf-8---
# date:20171123
# author:wuling
# emai:ling.wu@myhealthgene.com

# this model is set to
import sys
sys.path.append('../..')
sys.setdefaultencoding = ('utf-8')
from share import *
from config import *
from template import *

version = '1.0'

model_name = psplit(os.path.abspath(__file__))[1]

current_path = psplit(os.path.abspath(__file__))[0]

topics =  [name for name in listdir('./') if  not any([name.endswith(x) for x in ['.py','.pyc','.readme','.git']]) and not name.startswith('_')]

for topic in topics:

    import_topic = "from {} import  {}".format(topic,topic)

    exec(import_topic)

class manager(object):

    #his class is to mange all topics  under this directory 

    def __init__(self,topicname):

        self.topicname = topicname

    def helper(self):

        print  manage_help

    def initTopic(self):
    
        #this function is to init a new topic with specified topic name
        #topicname --- the specified topic's name
    
        # check to see if topicname  existed
        print '-'*50

        if self.topicname in topics:

            tips = 'the topic {} existed ,do  you still want to  create it ? (y/n) : '.format(self.topicname)

            choice = raw_input(tips)

            if choice == 'n':
                return
        # create major dir
        _topic = pjoin('./',self.topicname)
        createDir(_topic)

        # create moldename_v1.py
        pyload = open(pjoin(_topic,'{}.py'.format(self.topicname)),'w')
        pyload.write(py_template.replace('*'*6,self.topicname).strip() + '\n')
        pyload.close()

        initload = open(pjoin(_topic,'__init__.py'),'w')
        initload.close()

        introload = open(pjoin(_topic,'{}.readme'.format(self.topicname)),'w')
        introload.write(topic_intros.replace('*'*6,self.topicname) + '\n')
        introload.close()

        print 'topic %s  created successfully !' % self.topicname

    def importTopic(self,allupdate=False):
        '''
        this function is to return a update function of all topic under current directory
        topicname --- the specified topic's name
        allupdate -- default False,if set to true, update all topic one by one
        '''
        updates = {
        'gene':gene.update,
        }

        return updates if allupdate else updates.get(self.topicname)

    def updateTopic(self):
        '''
        this function is to update the specified mode 
        topicname ---the specified topic's name,if == 'all',all topic would be updated
        '''
        if self.topicname != 'all':

            if self.topicname not in topics:

                print 'No topic named {} '.format(self.topicname)
                sys.exit()

            else:
                try:
                    update_fun = self.importTopic(allupdate=False)

                    print '-'*50
                    update_fun()

                except Exception,e:
                    print e
        else:

            update_funs =self.importTopic(allupdate=True)

            n = 1
            for topic,fun in update_funs.items():

                try:
                    print '*'*80
                    print n,topic,'\n'
                    fun()
                except Exception,e:
                    print 'x'*50
                    print e
                    print 'x'*50
                n += 1

class query(object):

    """docstring for query"""
    
    def __init__(self):

        conn = MongoClient('localhost',27017)

        db  = conn.get_database('mytopic')

        self.db = db

    def selectFromTopic(self,topicname,value,outputdir):

        topic_field = {'gene':'symbol'}

        col_names =[col_name for col_name in  self.db.collection_names() if col_name.startswith(topicname)]

        col_names.sort(key = lambda x:x.split('_')[1].strip())

        # default selcet from the newest edition
        col_name = col_names[-1]

        col = self.db.get_collection(col_name)

        docs = col.find({topic_field[topicname]:value})

        if docs:

            createDir(outputdir)

            n = 0 

            for doc in docs:

                doc.pop('_id')

                with open(pjoin(outputdir,'{}_{}_{}.json'.format(topicname,value,n)),'w') as wf:

                    json.dump(doc,wf,indent=8)

                    n += 1

            print 'allfind:',n

        else:

            print 'no record'

def main():
    
    try:

        (opts,args) = getopt.getopt(sys.argv[1:],"hi:u:t:v:o:",['--help','--init=','--update=','--topicbase=','--field=','--value=','--output='])

        (base,val,out) = ("","","")

        for op,value in opts:

            man = manager(value)

            if op in ("-h","--help"):
                man.helper()

            elif op in ('-i','--init'):
                man.initTopic()

            elif op in ('-u','--update'):
                man.updateTopic()

            elif op in ('-t','--topicbase'):
                base = value

            elif op in ('-v','--value'):
                val = value

            elif op in ('-o','--output'):
                out = value

            if base  and val and out:
                man = query()
                man.selectFromTopic(base,val,out)

    except getopt.GetoptError:

        sys.exit()

if __name__ == '__main__':

    main()
