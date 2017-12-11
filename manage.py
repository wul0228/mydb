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
from ncbi_gene import ncbi_gene_v1
from ensembl_gene import ensembl_gene_v1
from go_gene import  go_gene_v1
from kegg_pathway import kegg_pathway_v1
from reactom_pathway import reactom_pathway_v1
from wiki_pathway import wiki_pathway_v1
from disgenet_disease import disgenet_disease_v1
from clinVar_variant import clinVar_variant_v1
from protein_atlas import protein_atlas_v1
from miRTarBase import miRTarBase_v1
from hpo_phenotypic import hpo_phenotypic_v1

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

    def importModel(self,allupdate=False):
        '''
        this function is to return a update function of all model under current directory
        modelname --- the specified model's name
        allupdate -- default False,if set to true, update all model one by one
        '''
        updates = {
        'ncbi_gene':ncbi_gene_v1.updateData,
        'ensembl_gene':ensembl_gene_v1.updateData,
        'kegg_pathway':kegg_pathway_v1.updateData,
        'reactom_pathway':reactom_pathway_v1.updateData,
        'wiki_pathway':wiki_pathway_v1.updateData,
        'disgenet_disease':disgenet_disease_v1.updateData,
        'clinVar_variant':clinVar_variant_v1.updateData,
        'protein_atlas':protein_atlas_v1.updateData,
        'miRTarBase':miRTarBase_v1.updateData,
        'hpo_phenotypic':hpo_phenotypic_v1.updateData,
        }

        return updates if allupdate else updates.get(self.modelname)

    def updateModel(self):
        '''
        this function is to update the specified mode 
        modelname ---the specified model's name,if == 'all',all model would be updated
        '''
        if self.modelname != 'all':

            if self.modelname not in models:

                print 'No model named {} '.format(self.modelname)
                sys.exit()

            else:
                try:
                    update_fun = self.importModel(allupdate=False)

                    print '-'*50
                    update_fun()

                except Exception,e:
                    print e
        else:

            update_funs =self.importModel(allupdate=True)

            n = 1
            for model,fun in update_funs.items():

                try:
                    print '*'*80
                    print n,model,'\n'
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

        db  = conn.get_database('mydb')

        self.db = db

    def selectFromModel(self,modelname,field,value,outputdir):

        col_names =[col_name for col_name in  self.db.collection_names() if col_name.startswith(modelname)]

        col_names.sort(key = lambda x:x.split('_')[1].strip())

        # default selcet from the newest edition
        col_name = col_names[-1]

        col = self.db.get_collection(col_name)

        docs = col.find({field:value})

        n = 0 

        for doc in docs:

            doc.pop('_id')

            with open(pjoin(outputdir,'{}_{}_{}_{}.json'.format(modelname,field,value,n)),'w') as wf:

                json.dump(doc,wf,indent=8)

                n += 1

def main():
    
    try:

        (opts,args) = getopt.getopt(sys.argv[1:],"hi:u:d:f:v:o:",['--help','--init=','--update=','--database=','--field=','--value=','--output='])

        (base,field,val,out) = ("","","","")

        for op,value in opts:

            man = manager(value)

            if op in ("-h","--help"):
                man.helper()

            elif op in ('-i','--init'):
                man.initModel()

            elif op in ('-u','--update'):
                man.updateModel()

            elif op in ('-d','--database'):
                base = value

            elif op in ('-f','--field'):
                field = value
                
            elif op in ('-v','--value'):
                val = value

            elif op in ('-o','--output'):
                out = value

            if base and field and val and out:
                man = query()
                man.selectFromModel(base,field,val,out)

    except getopt.GetoptError:

        sys.exit()

if __name__ == '__main__':

    main()
