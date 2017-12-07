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

models =  [name for name in listdir('./') if  not any([name.endswith(x) for x in ['.py','.pyc','.readme']])]

class manager(object):

   
    #his class is to mange all models  under this directory 

    def __init__(self,modelname):

        self.modelname = modelname

    def helper(self):

        print  manage_help

    def initModel(self):
    
        #this function is to init a new model with specified model name
        #modelname --- the specified model's name
    
        # check to see if modelname  existed
        print '-'*50

        if self.modelname in models:

            tips = 'the model {} existed ,do  you still want to  create it ? (y/n) : '.format(self.modelname)

            choice = raw_input(tips)

            if choice == 'n':
                return
        # create major dir
        createDir(pjoin('./',self.modelname))

        # create dataload,dataraw,datastore and database  
        (_model,_raw,_store,_db,_map) = buildSubDir(self.modelname)
        createDir(pjoin(_db,'docs'))

        # create moldename_v1.py
        pyload = open(pjoin(_model,'{}_v1.py'.format(self.modelname)),'w')
        pyload.write(py_template.replace('*'*6,self.modelname).strip() + '\n')
        pyload.close()

        initload = open(pjoin(_model,'__init__.py'),'w')
        initload.close()

        introload = open(pjoin(_model,'{}.readme'.format(self.modelname)),'w')
        introload.write(model_intros.replace('*'*6,self.modelname) + '\n')
        introload.close()

        print 'model %s  created successfully !' % self.modelname

def main():
    
    try:

        (opts,args) = getopt.getopt(sys.argv[1:],"hi:",['--help',"--init="])

        for op,value in opts:

            man = manager(value)

            if op in ("-h","--help"):
                man.helper()

            elif op in ('-i','--init'):
                man.initModel()

    except getopt.GetoptError:

        sys.exit()

if __name__ == '__main__':

    main()
