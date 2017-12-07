#!/usr/bin/env python
# --coding:utf-8--
# date: 20171123
# author:wuling
# emai:ling.wu@myhealthgene.com

'''
this model set  to download,extract,standard insert and select gene data from ncbi
'''

import sys
sys.path.append('../')
sys.setdefaultencoding = ('utf-8')
from share import *
from config import *  

__all__ = ['downloadData','extractData','updateData','selectData']

version  = 1.0

model_name = psplit(os.path.abspath(__file__))[1]

(ncbi_gene_model,ncbi_gene_raw,ncbi_gene_store,ncbi_gene_db,ncbi_gene_map) = buildSubDir('ncbi_gene')

# main code

def downloadOne(ncbi_gene_ftp_infos,filename,rawdir):
    '''
    this function is to download  one file under  a given remote dir 
    args:
    ftp -- a ftp cursor for a specified
    filename --  the name of file need download
    rawdir -- the directory to save download file
    '''
    while  True:

        try:

            ftp = connectFTP(**ncbi_gene_ftp_infos)

            mt =  ftp.sendcmd('MDTM {}'.format(filename)).replace(' ','')

            savefilename = '{}_{}_{}.gz'.format(filename.rsplit('.',1)[0].strip(),mt,today)

            remoteabsfilepath = pjoin(ncbi_gene_ftp_infos['logdir'],'{}'.format(filename))

            save_file_path = ftpDownload(ftp,filename,savefilename,rawdir,remoteabsfilepath)

            print filename,'done'

            return (save_file_path,mt)

        except:

            ftp = connectFTP(**ncbi_gene_ftp_infos)

def downloadData(redownload = False):
    '''
    this function is to download the raw data from ncbi gene FTP WebSite
    args:
    redownload-- default False, check to see if exists an old edition before download
                       -- if set to true, download directly with no check
    rawdir-- the directory to save download file
    '''
    if  not redownload:

        (choice,existNcbiFile) = lookforExisted(ncbi_gene_raw,'gene')

        if choice != 'y':
            return

    if redownload or not existNcbiFile or  choice == 'y':

        rawdir = pjoin(ncbi_gene_raw,'gene_{}'.format(today))

        createDir(rawdir)

        func = lambda x:downloadOne(ncbi_gene_ftp_infos,x,rawdir)

        # dodownload gene_group ,gene_neibors, gene_pubmed ,gene_info
        multiProcess(func,ncbi_gene_filenames,size=16)

        #download gene_expression 
        ncbi_gene_ftp_infos['logdir'] =  '/gene/DATA/expression/Mammalia/Homo_sapiens/'

        # get expression filenames
        ftp = connectFTP(**ncbi_gene_ftp_infos)

        gene_expression_filenames = ftp.nlst()
        
        # download
        multiProcess(func,gene_expression_filenames,size=16)

    # rawdir = '/home/user/project/dbproject/mygene_v1/ncbi/dataraw/gene_171124124057'

    if not os.path.exists(pjoin(ncbi_gene_model,'ncbi_gene.log')):

        initLogFile('ncbi_gene',model_name,ncbi_gene_model,rawdir=rawdir)

    update_file_heads =dict()

    for filename in listdir(rawdir):

        head = filename.split('_213')[0].strip()

        update_file_heads[head] = filename


    with open(pjoin(ncbi_gene_db,'gene_{}.files'.format(today)),'w') as wf:
        json.dump(update_file_heads,wf,indent=2)

    print  'datadowload completed !'

    filepaths = [pjoin(rawdir,filename) for filename in rawdir]

    return (filepaths,today)

def extractData(filepaths,version):

    file_index = dict()

    for filepath in filepaths:

        filename = psplit(filepath)[1].strip()

        if filename.startswith('gene_info'):

            process = gene_parser(filepath,version)

            process.gene_info()

            print filename,'done'

    for filepath in filepaths:

        filename = psplit(filepath)[1].strip()

        if not filename.startswith('gene_info'):

            process = gene_parser(filepath,version)

            head_fun = {
            'gene_group':process.gene_group,
            'gene2pubmed':process.gene_pubmed,
            'gene_neighbors':process.gene_neighbors,
            'PRJ':process.gene_expression
            }

            for head in head_fun.keys():
                if filename.startswith(head):
                    fun = head_fun.get(head)
                    fun()

            print filename,'done'

        print 'extract and insert complete '

        return (filepaths,version)

def updateData(insert=False):

    ncbi_gene_log = json.load(open('./ncbi_gene.log'))

    rawdir = pjoin(ncbi_gene_raw,'gene_update_{}'.format(today))

    expression_ftp_info = copy.deepcopy(ncbi_gene_ftp_infos)

    expression_ftp_info['logdir'] = ncbi_gene_expression_path

    ftp = connectFTP(**expression_ftp_info)

    gene_expression_filenames = ftp.nlst()

    new = False

    filenames = ncbi_gene_filenames + gene_expression_filenames

    for filename in filenames[::-1]:

        if filename.startswith('gene'):

           ftp = connectFTP(**ncbi_gene_ftp_infos)
           download_ftp_infos = ncbi_gene_ftp_infos

        else:
            # if no link with ftp for a while ,it may be break 
            ftp = connectFTP(**expression_ftp_info)
            download_ftp_infos = expression_ftp_info

        mt = ftp.sendcmd('MDTM {}'.format(filename)).strip()

        if mt != ncbi_gene_log.get(filename)[-1][0].strip():

            new = True

            createDir(rawdir)
            
            downloadOne(download_ftp_infos,filename,rawdir)

            ncbi_gene_log[filename].append((mt,today,model_name))

            print  '{} \'s new edition is {} '.format(filename,mt)

        else:
            print  '{} is the latest !'.format(filename)

    if new:

        with open(pjoin(ncbi_gene_model,'ncbi_gene.log'),'w') as wf:

            json.dump(ncbi_gene_log,wf,indent=2)

        (latest_file,version) = createNewVersion(ncbi_gene_raw,ncbi_gene_db,rawdir,'gene_',today)
        
        if insert:

            insertUpdatedData(ncbi_gene_raw,latest_file,'gene_',version,extractData)

def selectData(querykey = 'GeneID',value='1'):
    '''
    this function is set to select data from mongodb
    args:
    querykey -- a specified field in database
    queryvalue -- a specified value for a specified field in database
    '''
    conn = MongoClient('127.0.0.1', 27017 )

    db = conn.mygene

    colnamehead = 'ncbi'

    dataFromDB(db,colnamehead,querykey,queryvalue=None)

class dbMap(object):

    #class introduction

    def __init__(self,version):

        self.version = version

        conn = MongoClient('127.0.0.1',27017)

        db = conn.get_database('mygene')

        col = db.get_collection('ncbi_{}'.format(self.version))

        docs = col.find({})

        self.docs  = docs

    def mapGeneID2Symbol(self):
        '''
        this function is to create a mapping relation between GeneID with Symbol
        '''
        id_symbol = dict()

        symbol_id = dict()

        n = 0

        for doc in self.docs:

            gene_id = doc.get('GeneID')

            symbol = doc.get('Symbol')

            # add Synonyms
            synonyms = doc.get('Synonyms')

            if synonyms:

                synonyms.append(symbol) 

                symbols =list(set( [ i for i in synonyms if i != '-']))

            else:
                symbols = [symbol,]

            id_symbol[gene_id] = symbols

            for s in symbols:

                if s not in symbol_id:

                    symbol_id[s] = list()
                
                symbol_id[s].append(gene_id)

            n += 1

            print n,gene_id,symbol

        with open(pjoin(ncbi_gene_map,'geneid2symbol_{}.json'.format(self.version)),'w') as wf:
            json.dump(id_symbol,wf,indent=2)
        with open(pjoin(ncbi_gene_map,'symbol2geneid_{}.json'.format(self.version)),'w') as wf:
            json.dump(symbol_id,wf,indent=2)

    def mapping(self):

        self.mapGeneID2Symbol()


class gene_parser(object):
    
    def __init__(self,filepath,version):

        rawdir = psplit(filepath)[0].strip()

        filename = psplit(filepath)[1].strip()

                # gunzip file
        if filename.endswith('.gz'):

            command = 'gunzip  {}'.format(filepath)
            
            os.popen(command)

            filename = filename.replace('.gz','')

        tsvfile = open(pjoin(rawdir,filename))

        self.rawdir = rawdir

        self.filename = filename

        self.tsvfile = tsvfile

        self.version = version

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mygene')
        # db.authenticate('wul','wul@408')
        col = db.get_collection('ncbi_gene_{}'.format(version))

        self.col = col

        print 'init completed !'

    def gene_info(self):

        print '-'*50

        n =0

        gene_num = 0

        for line in self.tsvfile:

            data = line.strip().split('\t')

            if n == 0:
                keys =[ key.replace('#','').strip().replace('.','*').replace(' ','&') for key  in  data]
 
            else:
                data = dict(zip(keys,data))

                tax_id = data.get('tax_id')

                if tax_id != '9606':
                    continue

                # alter synonyms,dbxrefs,other_designations and feature_type
                '''
                Synonyms  list   '|',dbXrefs dict '|',':',Other_designations list '|',feature_type dict '|',':'
                '''
                data['Synonyms'] = data['Synonyms'].strip().split('|')
                data['Other_designations'] = data['Other_designations'].strip().split('|')

                dbXrefs = data['dbXrefs']

                if dbXrefs  != '-':
                    data['dbXrefs']  = dict([tuple(refs.split(':',1)) for refs in dbXrefs.split('|')])

                self.col.insert(data)

                gene_num += 1

                print 'gene_info',n,gene_num,data['GeneID']

            n += 1

        print 'gene_info completed'

    def gene_group(self):

        print '-'*50

        n =0

        # gene_group = dict()
        for line in self.tsvfile:

            data = line.strip().split('\t')
    
            if n == 0:
                keys =[ key.replace('#','').strip().replace('.','*').replace(' ','&') for key  in  data]

            else:
                data = dict(zip(keys,data))

                # filter line not human's
                other_tax_id = data['Other_tax_id']

                tax_id = data['tax_id']

                gene_id = data['GeneID']

                if  tax_id == '9606':
    
                    group = {
                    'relationship':data['relationship'],
                    'Other_tax_id':other_tax_id,
                    'Other_GeneID':data['Other_GeneID']
                     }

                elif other_tax_id == '9606':
    
                    gene_id = data['Other_GeneID']

                    group = {
                    'relationship':data['relationship'],
                    'Other_tax_id':tax_id,
                    'Other_GeneID':gene_id,
                     }

                else:
                    continue

                print 'gene_group',n,gene_id

                self.col.update({'GeneID':gene_id},{'$push':{'gene_group':group}})

            n += 1

        print 'gene_group completed'

    def gene_neighbors(self):

        print '-'*50

        n =0

        for line in self.tsvfile:

            data = line.strip().split('\t')
    
            if n == 0:
                keys =[ key.replace('#','').strip().replace('.','*').replace(' ','&') for key  in  data]

            else:

                data = dict(zip(keys,data))

                tax_id = data['tax_id']

                assembly = data['assembly']

                # filter tax_id
                if tax_id != '9606' or assembly == '-':
                    continue
                
                gene_id = data['GeneID']
         
                # delete tax_id assembly gene_id
                neighbors = dict()

                for key,val in data.items():

                    if key in ['tax_id','assembly','GeneID']:

                        continue
           
                    key = key.strip().replace(' ','&').replace('.','*')

                    if val.count('|'):

                        val = val.strip().split('|')

                    neighbors[key] = val     

                    assembly = assembly.strip().replace(' ','&').replace('.','*')

                
                self.col.update({'GeneID':gene_id},{'$push':{'gene_neighbors.{}'.format(assembly):neighbors}})

                print 'gene_neighbors',n,gene_id #112312

            n += 1

        print 'gene_neighbors completed'

    def gene_pubmed(self):

        print '-'*50

        n =0

        for line in self.tsvfile:

            data = line.strip().split('\t')
    
            if n == 0:
                keys =[ key.replace('#','').strip().replace('.','*').replace(' ','&') for key  in  data]

            else:
                data = dict(zip(keys,data))

                # filter line not human's

                tax_id = data['tax_id']

                if  tax_id != '9606':
                    continue

                gene_id = data['GeneID']
                pubmed_id = data['PubMed_ID']

                self.col.update({'GeneID':gene_id},{'$push':{'PubMed_ID':pubmed_id}})

                print 'gene_pubmed',n ,gene_id

            n += 1

        print 'gene_pubmed completed'

    def gene_expression(self):

        print '-'*50

        project_desc = self.filename.split('_',1)[0].strip()

        xmlfile = self.tsvfile

        sample_source = dict()

        # can not parse by xml or lxml ,so deal as a txt file
        aset = dict()

        n = 0

        for line in xmlfile:

            line = line.strip()

            if not line or line == '<doc>':
                continue

            elif line == '</doc>' or line == '</doc><doc>':

                _id  = aset.get('id')

                if _id.startswith('metadata'):

                    doc_type = 'sample_infos'

                    sra_id = aset.get("sra_id")
                    sample_id = aset.get("sample_id")
                    source_name = aset.get("source_name").replace(',','_')

                    sample_source[sample_id] = {'source_name':source_name,'sra_id':sra_id}

                elif _id.split('_',1)[1].strip().startswith('SAM'):

                    doc_type = 'gene_infos'

                    gene_id = aset.get('gene')

                    sample_id = aset.get("sample_id")
                    sample = sample_source.get(sample_id)
                    if sample:
                        source_name = sample.get('source_name')
                                    #aset['source_name'] = source_name
                        
                        sra_id = sample.get('sra_id')
                        if sra_id :
                            aset['sra_id'] = sra_id

                    aset.pop('id')
                    aset.pop('gene')
                    aset.pop("project_desc")
                    if 'source_name' in aset:
                        aset.pop('source_name')

                    print 'gene_expression',n,gene_id
                    self.col.update({'GeneID':gene_id},{'$push':{'project.{}.{}'.format(project_desc,source_name):aset}})

                # dict clear
                aset = dict()

            else:
                key_val = line.split('<field name="')[1].split('</field>')[0].strip()
                (key,val) = tuple(key_val.split('">'))
                aset[key] = val

            n += 1
    
        print self.filename,'completed'

        # PRJEB2445_GRCh38.p2_107_expression.xml    / 11191815
        # PRJNA280600_GRCh38.p7_108_expression.xml completed /  6510332
        # PRJNA270632_GRCh38.p7_108_expression.xml /17754426
        # PRJEB4337_GRCh38.p7_108_expression.xml /37083157
            
def main():

    modelhelp = model_help.replace('&'*6,'NCBI_GENE').replace('#'*6,'ncbi_gene')

    funcs = (downloadData,extractData,updateData,selectData,dbMap,ncbi_gene_store)

    getOpts(modelhelp,funcs=funcs)
        
if __name__ == '__main__':
    main()
    

