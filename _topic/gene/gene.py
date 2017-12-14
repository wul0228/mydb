#!/usr/bin/env python
# --coding:utf-8--
# date:2017/12/13
# author:wuling
# emai:ling.wu@myhealthgene.com

'''
this model is set  to construct a topic base
'''
import sys
reload(sys)
sys.path.append('../')
sys.setdefaultencoding = ('utf-8')
from config import *
from share import *

current_path = psplit(os.path.abspath(__file__))[0]

class dbMap(object):

    """docstring for dbMap"""

    def __init__(self, ):

        super(dbMap, self).__init__()
    
        (db,db_cols,db_colnames) = initDB()

        self.db = db

        self.db_cols = db_cols

        self.db_colnames = db_colnames

        self.hgnc_col = self.db_cols.get('hgnc_gene')

        map_path = pjoin(current_path,'map_{}'.format(today))

        if not os.path.exists(map_path):
            os.mkdir(map_path)

        self.map_path = map_path

    def ncbi2hgnc(self):

        ncbi_col = self.db_cols.get('ncbi_gene')

        ncbi_docs = ncbi_col.find({})

        ncbi2hgnc = dict()

        n = 0

        for ncbi_doc in ncbi_docs:

            gene_id = ncbi_doc.get('GeneID')

            hgnc_docs= self.hgnc_col.find({"entrez_id":gene_id})

            if hgnc_docs:

                if gene_id not in ncbi2hgnc:

                    ncbi2hgnc[gene_id] = list()

                for hgnc_doc in hgnc_docs:

                    hgnc_symbol = hgnc_doc.get("symbol")

                    ncbi2hgnc[gene_id].append(hgnc_symbol)

            n += 1

            print 'ncbi doc',n

        with open(pjoin(self.map_path,'./ncbi2hgnc.json'),'w') as wf:
            json.dump(ncbi2hgnc,wf,indent=8)

        print 'ncbi2hgnc completed'

        return ncbi2hgnc

    def ensembl2hgnc(self):

        ensembl_col = self.db_cols.get('ensembl_gene')

        ensembl_docs = ensembl_col.find({})

        ensembl2hgnc = dict()

        n = 0

        for ensembl_doc in ensembl_docs:

            gene_id = ensembl_doc.get('gene_id')

            hgnc_docs= self.hgnc_col.find({"ensembl_gene_id":gene_id})

            if hgnc_docs:

                if gene_id not in ensembl2hgnc:

                    ensembl2hgnc[gene_id] = list()

                for hgnc_doc in hgnc_docs:

                    hgnc_symbol = hgnc_doc.get("symbol")

                    ensembl2hgnc[gene_id].append(hgnc_symbol)

            n += 1

            print 'ensembl doc',n

        with  open(pjoin(self.map_path,'./ensembl2hgnc.json'),'w') as wf:
            json.dump(ensembl2hgnc,wf,indent=8)

        print 'ensembl2hgnc completed'

        return ensembl2hgnc

    def uniprot2hgnc(self):

        hgnc_docs =self.hgnc_col.find({})

        uniprot2hgnc = dict()

        n = 0
        
        for hgnc_doc in hgnc_docs:

            hgnc_symbol = hgnc_doc.get("symbol")

            gene_id= hgnc_doc.get("uniprot_ids")

            if not gene_id:

                continue

            if gene_id not in uniprot2hgnc:

                uniprot2hgnc[gene_id] = list()

            uniprot2hgnc[gene_id].append(hgnc_symbol)

            n += 1

            print 'hgnc  doc',n

        with  open(pjoin(self.map_path,'./uniprot2hgnc.json'),'w') as wf:
            json.dump(uniprot2hgnc,wf,indent=8)

        print 'uniprot2hgnc completed'

        return uniprot2hgnc

    def go2hgnc(self):

        go_col = self.db_cols.get('go_gene')

        go_docs = go_col.find({})

        go2hgnc = dict()

        n = 0
        
        for go_doc in go_docs:

            go_ids = go_doc.get('GO',{}).keys()

            if not go_ids:

                continue

            for go_id in go_ids:

                if go_id not in go2hgnc:

                    go2hgnc[go_id] = list()

            gene_id = go_doc.get("DB_Object_ID")

            hgnc_docs= self.hgnc_col.find({"uniprot_ids":gene_id})

            if hgnc_docs:

                for hgnc_doc in hgnc_docs:

                    hgnc_symbol = hgnc_doc.get("symbol")

                    for go_id in go_ids:

                        go2hgnc[go_id].append(hgnc_symbol)

            for key,val in go2hgnc.items():

                go2hgnc[key] = list(set(val))

            n += 1

            print 'go  doc',n

        with  open(pjoin(self.map_path,'./go2hgnc.json'),'w') as wf:
            json.dump(go2hgnc,wf,indent=8)

        print 'go2hgnc completed'

        return go2hgnc

    def proteinatlas2hgnc(self):
        '''
        proteinatlas's identifier is ensembl id 
        '''
        proteinatlas_col = self.db_cols.get('protein_atlas')

        proteinatlas_docs = proteinatlas_col.find({})

        proteinatlas2hgnc = dict()

        n = 0

        for  proteinatlas_doc in  proteinatlas_docs:

            gene_id = proteinatlas_doc.get('Ensembl')

            hgnc_docs= self.hgnc_col.find({"ensembl_gene_id":gene_id})

            if hgnc_docs:

                if gene_id not in proteinatlas2hgnc:

                    proteinatlas2hgnc[gene_id] = list()

                for hgnc_doc in hgnc_docs:

                    hgnc_symbol = hgnc_doc.get("symbol")

                    proteinatlas2hgnc[gene_id].append(hgnc_symbol)


            n += 1

            print 'proteinatlas doc',n

        with  open(pjoin(self.map_path,'./proteinatlas2hgnc.json'),'w') as wf:
            json.dump(proteinatlas2hgnc,wf,indent=8)

        print 'proteinatlas2hgnc completed'

        return proteinatlas2hgnc

    def kegg2hgnc(self):

        kegg_col = self.db_cols.get('kegg_pathway')

        kegg_docs = kegg_col.find({})

        kegg2hgnc = dict()

        n  = 0

        for doc in kegg_docs:

            path_id = doc.get('path_id')

            genes = doc.get('gene',{}).keys()

            if genes:

                for gene_id in genes:

                    hgnc_docs= self.hgnc_col.find({"entrez_id":gene_id})

                    if hgnc_docs:

                        if path_id not in kegg2hgnc:

                            kegg2hgnc[path_id] = list()

                        for hgnc_doc in hgnc_docs:

                            hgnc_symbol = hgnc_doc.get("symbol")

                            kegg2hgnc[path_id].append(hgnc_symbol)
            n += 1
            print 'kegg2hgnc doc',n

        with  open(pjoin(self.map_path,'./kegg2hgnc.json'),'w') as wf:
            json.dump(kegg2hgnc,wf,indent=8)

        print 'kegg2hgnc completed'

        return kegg2hgnc     

    def reactom2hgnc(self):

        reactom_col = self.db_cols.get('reactom_pathway')

        reactom_docs = reactom_col.find({})

        reactom2hgnc = dict()

        n  = 0

        allgeneid = list()

        for doc in reactom_docs:

            path_id = doc.get('stId')

            genes = doc.get('nodes',{}).get('EntityWithAccessionedSequence',{})

            if genes:

                gene_ids = [val.get('identifier') for dbid,val in genes.items()]

                if 'ENST00000387405' in gene_ids:
                    print path_id

                allgeneid += gene_ids
                    # if hgnc_docs:

                    #     if path_id not in kegg2hgnc:

                    #         kegg2hgnc[path_id] = list()

                    #     for hgnc_doc in hgnc_docs:

                    #         hgnc_symbol = hgnc_doc.get("symbol")

                    #         kegg2hgnc[path_id].append(hgnc_symbol)
            # n += 1
            # print 'kegg2hgnc doc',n
        allgeneid = list(set(allgeneid))

        with open('./rreactom_allgeneid.json','w') as wf:
            json.dump(allgeneid,wf,indent=4)
        # with  open(pjoin(self.map_path,'./kegg2hgnc.json'),'w') as wf:
        #     json.dump(kegg2hgnc,wf,indent=8)

        # print 'kegg2hgnc completed'

        # return kegg2hgnc         

    def wiki2hgnc(self):
        pass
        
    def clinvar2hgnc(self):
        pass

    def disgenet2hgnc(self):
        pass

    def miRTarBase2hgnc(self):
        pass

    def hpo2hgnc(self):
        pass

    def mapping(self):

        proteinatlas2hgnc = self.proteinatlas2hgnc()
        ncbi2hgnc = self.ncbi2hgnc()
        ensembl2hgnc = self.ensembl2hgnc()
        go2hgnc = self.go2hgnc()
        uniprot2hgnc =  self.uniprot2hgnc()
        kegg2hgnc = selif.kegg2hgnc()

        pairs = {'ncbi_gene':ncbi2hgnc,'ensembl_gene':ensembl2hgnc,'go_gene':go2hgnc,'uniprot_gene':uniprot2hgnc,
                    'protein_atlas':proteinatlas2hgnc,'kegg_pathway':kegg2hgnc}

        main_keys = {'ncbi_gene':'GeneID','ensembl_gene':'gene_id','go_gene':'id','uniprot_gene':'DB_Object_ID','protein_atlas':'Ensembl'
                            'kegg_pathway':'path_id'}

        hgnc2dbs = dict()

        for key,val in pairs.items():

            for _id,hgncs in val.items():
                
                if not hgncs:
                    continue

                for hgnc in hgncs:

                    if hgnc not in hgnc2dbs:
                        hgnc2dbs[hgnc]= dict()

                    if key not in hgnc2dbs[hgnc]:

                        hgnc2dbs[hgnc][key] = dict()

                        hgnc2dbs[hgnc][key]['field'] = main_keys[key]

                        hgnc2dbs[hgnc][key]['vals'] = list()

                    hgnc2dbs[hgnc][key]['vals'].append(_id)

        savefilepath = pjoin(self.map_path,'./hgnc2dbs.json')
        savemappath = pjoin(self.map_path,'map.log')

        with open(savefilepath,'w') as wf:
            json.dump(hgnc2dbs,wf,indent=8)

        with open(savemappath,'w') as wf:
            json.dump(self.db_colnames,wf,indent=8)

        print 'hgnc symbol:',len(hgnc2dbs)   # 41182 # hgnc db 41283

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

        hgnc2dbs = json.load(open(self.filepath))

        map_log = json.load(open(self.maplogpath))

        mydb = self.conn.get_database('mydb')

        model_cols = dict()

        for model,colname in map_log.items():

            model_cols[model] = mydb.get_collection(colname)

        n = 0

        for symbol,symbol_val in hgnc2dbs.items():

            for model,model_val in symbol_val.items():

                if model == 'uniprot_gene':
                    
                    model = 'go_gene'

                    print symbol,model

                    col = model_cols.get(model)

                    field = model_val['field']
                    vals = model_val.get('vals')

                    for val in vals:
                        docs = col.find({field:val})

                        for doc in docs:

                            doc.pop('_id')

                            self.col.update(
                                {'symbol':symbol},
                                {'$push':{'{}'.format(model):doc}},
                                True
                            )
            n += 1
            # print n,symbol

def update():

    man = dbMap()
    (filepath,maplogpath) = man.mapping()

    man = dbCreate(filepath,maplogpath)
    man.build()

def select(queryvalue='SNRPD3'):
    '''
    this function is set to select data from mongodb
    args:
    querykey -- a specified field in database
    queryvalue -- a specified value for a specified field in database
    '''
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
    # man = dbMap()
    # man.reactom2hgnc()
    # filepath = '/home/user/project/dbproject/mydb_v1/_topic/gene/map_171214092209/hgnc2dbs.json'
    # maplogpath = '/home/user/project/dbproject/mydb_v1/_topic/gene/map_171214092209/map.log'
    # man = dbCreate(filepath,maplogpath)
    # man.build()

    # hgnc2dbs = json.load(open('/home/user/project/dbproject/mydb_v1/_topic/gene/map_171214092209/hgnc2dbs.json'))
    # kegg2hgnc = json.load(open('/home/user/project/dbproject/mydb_v1/_topic/gene/map_171214092209/kegg2hgnc.json'))

    # for path_id,hgncs in kegg2hgnc.items():

    #     for hgnc in hgncs:

    #         if hgnc not in hgnc2dbs:
    #             hgnc2dbs[hgnc] = dict()

    #         if 'kegg_pathway' not in hgnc2dbs[hgnc]:
    #             hgnc2dbs[hgnc]['kegg_pathway'] = dict()
    #             hgnc2dbs[hgnc]['kegg_pathway']['field'] = 'path_id'
    #             hgnc2dbs[hgnc]['kegg_pathway']['vals'] = list()

    #         hgnc2dbs[hgnc]['kegg_pathway']['vals'].append(path_id)


    # with open('./hgnc2dbs.json','w') as wf:
    #     json.dump(hgnc2dbs,wf,indent=4)

