#!/usr/bin/env python
# --coding:utf-8--
# date: 2017/11/28
# author:wuling
# emai:ling.wu@myhealthygene.com

#this model set  to download,extract,standard insert and select gene data from ensembl

import sys
sys.path.append('../')
sys.setdefaultencoding = ('utf-8')
from share import *
from config import *  

__all__ = ['downloadData','extractData','standarData','insertData','updateData','selectData']

version  = 1.0

model_name = psplit(os.path.abspath(__file__))[1]

(ensembl_gene_model,ensembl_gene_raw,ensembl_gene_store,ensembl_gene_db,ensembl_gene_map) = buildSubDir('ensembl_gene')

# main code
def downloadOne(ensembl_gene_ftp_infos,filename,rawdir):
    '''
    this function is to download  one file under  a given remote dir 
    args:
    ftp -- a ftp cursor for a specified
    filename --  the name of file need download
    rawdir -- the directory to save download file
    '''
    while  True:

        try:

            ftp = connectFTP(**ensembl_gene_ftp_infos)

            mt =  ftp.sendcmd('MDTM {}'.format(filename))

            print 'mt',mt

            savefilename = '{}_{}_{}.gz'.format(filename.rsplit('.',1)[0].strip(),mt,today).replace(' ','')

            remoteabsfilepath = pjoin(ensembl_gene_ftp_infos['logdir'],'{}'.format(filename))

            print filename,'start'

            save_file_path = ftpDownload(ftp,filename,savefilename,rawdir,remoteabsfilepath)

            print filename,'done'

            return (save_file_path,mt)

        except:

            ftp = connectFTP(**ensembl_gene_ftp_infos)

def downloadData(redownload = False):

    '''
    this function is to download the raw data from ensembl gene FTP WebSite
    args:
    redownload-- default False, check to see if exists an old edition before download
                       -- if set to true, download directly with no check
    rawdir-- the directory to save download file
    '''
    if  not redownload:

        (choice,existensemblFile) = lookforExisted(ensembl_gene_raw,'gene')

        if choice != 'y':
            return

    if redownload or not existensemblFile or  choice == 'y':

        rawdir = pjoin(ensembl_gene_raw,'gene_{}'.format(today))

        createDir(rawdir)

        # download gtfGRch38 file
        ensembl_gene_ftp_infos['logdir'] = ensembl_gtfGRch38_ftp_path

        downloadOne(ensembl_gene_ftp_infos,filename_gtfGRch38,rawdir)

        ensembl_gene_ftp_infos['logdir'] = ensembl_regulatorGRch38_ftp_path

        func = lambda x:downloadOne(ensembl_gene_ftp_infos,x,rawdir)

        # download regulatorGRch38 files
        multiProcess(func,filenames_regulatorGRch38,size=16)

          # download gtfGRch37 file
        ensembl_gene_ftp_infos['logdir'] = ensembl_gtfGRch37_ftp_path

        downloadOne(ensembl_gene_ftp_infos,filename_gtfGRch37,rawdir)

        ensembl_gene_ftp_infos['logdir'] = ensembl_regulatorGRch37_ftp_path

        func = lambda x:downloadOne(ensembl_gene_ftp_infos,x,rawdir)

        # download regulatorGRch38 files
        multiProcess(func,filenames_regulatorGRch37,size=16)

    # create log file
    rawdir = '/home/user/project/dbproject/mygene_v1/ensembl_gene/dataraw/gene_171127151418'
    if not os.path.exists(pjoin(ensembl_gene_model,'ensembl_gene.log')):

        initLogFile('ensembl_gene',model_name,ensembl_gene_model,rawdir=rawdir)

    # create every version files included
    update_file_heads =dict()

    for filename in listdir(rawdir):

        head = filename.split('_213')[0].strip()

        update_file_heads[head] = filename

    with open(pjoin(ensembl_gene_db,'gene_{}.files'.format(today)),'w') as wf:
        json.dump(update_file_heads,wf,indent=2)

    print  'datadowload completed !'

    # generate filepaths to next step extractData
    filepaths = [pjoin(rawdir,filename) for filename in rawdir]

    return (filepaths,today)

def extractData(filepaths,version):

    #  deal the gtf file to  insert basic info
    for filepath in filepaths:

        name = psplit(filepath)[1].strip().split('_sapiens.')[1].split('_213')[0].strip()

        grch = name.split('.',1)[0]

        process = gene_parser(filepath,version,grch)

        if  name.endswith('chr.gtf'):

            process.gtf()

        print name,'completed'

    # inset the regulatory and motif data according to the range of gene start and  end 
    for filepath in filepaths:

        name = psplit(filepath)[1].strip().split('_sapiens.')[1].split('_213')[0].strip()

        # skip gtf file
        if  name.endswith('chr.gtf'):
            continue

        grch = name.split('.',1)[0]

        process = gene_parser(filepath,version,grch)

        if name.count('Regulatory_Build'):

            process.regulatory()

        elif name.count('motiffeatures'):

            process.motiff()

        print name,'completed'

def updateData(insert=False):

    ensembl_gene_log = json.load(open('./ensembl_gene.log'))

    rawdir = pjoin(ensembl_gene_raw,'gene_update_{}'.format(today))

    new = False

    for  file,ftpsite in ensembl_file_ftplogdir.items():

        ftp_infos = copy.deepcopy(ensembl_gene_ftp_infos)

        ftp_infos['logdir'] = ftpsite

        ftp = connectFTP(**ftp_infos)

        filenames = ftp.nlst()

        for filename in filenames:

            if filename.count(ensembl_file_mark.get(file)):

                mt = ftp.sendcmd('MDTM {}'.format(filename))

                if mt != ensembl_gene_log.get(file)[-1][0]:

                    new = True

                    createDir(rawdir)
                
                    downloadOne(ftp_infos,filename,rawdir)   

                    ensembl_gene_log[file].append((mt,today,model_name))

                    print  '{} \'s new edition is {} '.format(filename,mt)

                else:
                    print  '{} is the latest !'.format(filename)

        print '~'*50

    if new:

        with open(pjoin(ensembl_gene_model,'ensembl_gene.log'),'w') as wf:

            json.dump(ensembl_gene_log,wf,indent=2)

        (latest_file,version) = createNewVersion(ensembl_gene_raw,ensembl_gene_db,rawdir,'gene_',today)
        
        if insert:

            insertUpdatedData(ensembl_gene_raw,latest_file,'gene_',version,extractData)

def selectData(querykey = 'gene_id',value='ENSG00000243485'):

    '''
    this function is set to select data from mongodb
    args:
    querykey -- a specified field in database
    queryvalue -- a specified value for a specified field in database
    '''
    conn = MongoClient('127.0.0.1', 27017 )

    db = conn.mygene

    colnamehead = 'ensembl'

    dataFromDB(db,colnamehead,querykey,queryvalue=None)

class dbMap(object):

    #class introduction

    def __init__(self,version):

        self.version = version

        conn = MongoClient('127.0.0.1',27017)

        db = conn.get_database('mygene')

        col = db.get_collection('ensembl_{}'.format(self.version))

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

            gene_id = doc.get('gene_id')

            symbol = doc.get('gene_name')

            id_symbol[gene_id] = symbol

            if symbol not in symbol_id:

                symbol_id[symbol] = list()
            
            symbol_id[symbol].append(gene_id)

            n += 1

            print n,gene_id,symbol

        for key,val in symbol_id.items():
            if len(val) >= 2:
                print key,val

        with open(pjoin(ensembl_gene_map,'geneid2symbol_{}.json'.format(self.version)),'w') as wf:
            json.dump(id_symbol,wf,indent=2)
        with open(pjoin(ensembl_gene_map,'symbol2geneid_{}.json'.format(self.version)),'w') as wf:
            json.dump(symbol_id,wf,indent=2)

    def mapping(self):

        self.mapGeneID2Symbol()


class gene_parser(object):
    """docstring for gene_parser"""
    def __init__(self, filepath,version,grch):

        self.filepath = filepath

        self.version = version

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mygene')

        col = db.get_collection('ensembl_gene_{}'.format(version))

        self.col = col

        rawdir = psplit(filepath)[0].strip()

        filename = psplit(filepath)[1].strip()

        # gunzip file
        if filename.endswith('.gz'):

            command = 'gunzip  {}'.format(filepath)
            
            os.popen(command)

            filename = filename.replace('.gz','')

        file = open(pjoin(rawdir,filename))

        self.file = file

        self.filename = filename

        self.rawdir = rawdir

        self.grch = grch

        self.ensembl = 'ensembl_{}'.format(self.grch)

        print 'init completed'

    def gtf(self):
        
        print '-'*50

        n = 0

        for line in self.file:

            if line.startswith('#'):
                continue

            # the front of line ,delimited by tab and have no key, and the latter delimited by ; with the format key space value("") 

            front_keys = ['chr','data_source','entry','start','end','score','strand','fields']

            front = line.split('gene_id')[0].strip().split('\t')

            front_dic = dict([(key,val) for key,val in zip(front_keys,front)])

            # transform the string to int in start and end
            front_dic['start' ] = int(front_dic['start' ])
            front_dic['end' ] = int(front_dic['end' ])

            latter = [i.strip() for i in line.strip().split('gene_id')[1].strip().split(';') if i ]

            latter_dic = dict([(i.split(' ')[0],i.split(' ')[1].replace('"','')) for i in latter[1:] ])

            gene_id = latter[0].replace('"','')

            entry = front_dic.get('entry')

            if  entry == 'gene':

                gene_name = latter_dic.get('gene_name')

                latter_dic.update(front_dic)

                for key in ['data_source','entry','score','fields','gene_name',]:

                    latter_dic.pop(key)

                self.col.update(
                    {'gene_id':gene_id},
                    {'$set':{
                                'gene_name':gene_name,
                                '{}'.format(self.ensembl):latter_dic
                                }
                    },
                    True
                )
                # clear gene_aset
                gene_aset = dict()

            elif entry == 'transcript':

                latter_dic['transcript_start'] = front_dic.get('start')
                latter_dic['transcript_end'] = front_dic.get('end')

                transcript_id = latter_dic.get('transcript_id')

                for key in latter_dic.keys():
                    if key.startswith('gene'):
                        latter_dic.pop(key)

                latter_dic.pop('transcript_id')

                self.col.update({'gene_id':gene_id},{'$set':{'{}.transcript.{}'.format(self.ensembl,transcript_id):latter_dic}})

            elif entry in ['Selenocysteine','five_prime_utr','three_prime_utr']:

                transcript_id = latter_dic.get('transcript_id')

                entry_aset = dict()

                entry_aset['start'] = front_dic.get('start')
                entry_aset['end'] = front_dic.get('end')

                self.col.update({'gene_id':gene_id},{'$set':{'{}.transcript.{}.{}'.format(self.ensembl,transcript_id,entry):entry_aset}})

            elif entry == 'exon':

                transcript_id = latter_dic.get('transcript_id')

                for key in latter_dic.keys():

                    if not key.startswith('exon'):
                        latter_dic.pop(key)

                latter_dic['exon_start'] = front_dic.get('start')
                latter_dic['exon_end'] = front_dic.get('end')
                exon_number = latter_dic["exon_number"]
                latter_dic.pop('exon_number')

                self.col.update({'gene_id':gene_id},{'$set':{'{}.transcript.{}.exon.{}'.format(self.ensembl,transcript_id,exon_number):latter_dic}})

            elif entry == 'CDS':

                cds_aset = dict()

                transcript_id = latter_dic.get('transcript_id')

                exon_number = latter_dic["exon_number"]

                cds_aset['cds_start'] = front_dic.get('start')

                cds_aset['cds_end'] = front_dic.get('end')

                cds_aset['protein_id'] = latter_dic['protein_id']

                cds_aset['protein_version'] = latter_dic['protein_version']

                self.col.update({'gene_id':gene_id},{'$push':{'{}.transcript.{}.exon.{}.cds'.format(self.ensembl,transcript_id,exon_number):cds_aset}})
            
            elif entry == 'start_codon' or entry == 'stop_codon':

                transcript_id = latter_dic.get('transcript_id')

                exon_number = latter_dic["exon_number"]

                self.col.update({'gene_id':gene_id},{'$set':{'{}.transcript.{}.exon.{}.{}'.format(self.ensembl,transcript_id,exon_number,entry):(front_dic.get('start'),front_dic.get('end'))}})
 
            else:
                print '================================',entry

            latter_dic = dict()
            
            n += 1
            print self.grch,'gtf line',n,entry,gene_id

    def regulatory(self):

         n = 0

         for line in self.file:

            front = line.split('ID=')[0].strip().split('\t')

            Chr = front[0].strip()

            latter = [i.strip() for i in line.strip().split('ID=')[1].strip().split(';') if i ]

            ID = latter[0].strip()

            latter_dic = dict([(i.split('=')[0],i.split('=')[1]) for i in latter[1:] ])

            latter_dic.update({'ID':ID})

            start = int(latter_dic['bound_start'])
            end = int(latter_dic['bound_end'])

            latter_dic['bound_end'] = end
            latter_dic['bound_start'] = start

            self.col.update(
                {
                '{}.chr'.format(self.ensembl):Chr,
                '{}.end'.format(self.ensembl):{'$gt':end},
                '{}.start'.format(self.ensembl):{'$lt':start}
                },
                {'$push':{'{}.regulatory'.format(self.ensembl):latter_dic}},
                upsert=False,
                multi=True
                )
            n += 1

            print self.grch,'regulatory line',n

    def motiff(self):

        n = 0

        for line in self.file:

            front_keys = ['chr','data_source','entry','start','end','score','strand','fields']

            front = line.split('binding_matrix=')[0].strip().split('\t')
            
            front_dic = dict([(key,val) for key,val in zip(front_keys,front)])

            latter = [i.strip() for i in line.strip().split('binding_matrix=')[1].strip().split(';') if i ]

            latter_dic = dict([(i.split('=')[0],i.split('=')[1]) for i in latter[1:] ])

            latter_dic['binding_matrix'] = latter[0]

            Chr = front_dic['chr' ]
            motiff_start = int(front_dic['start' ])
            motiff_end = int(front_dic['end' ])
            latter_dic['score' ] = float(front_dic['score' ])
            # transform the string to int in start and end
            latter_dic['motif_start' ] = motiff_start
            latter_dic['motif_end' ] = motiff_end
            latter_dic['strand' ] = front_dic['strand' ]

            self.col.update({
                '{}.chr'.format(self.ensembl):Chr,
                '{}.end'.format(self.ensembl):{'$gt':motiff_end},
                '{}.start'.format(self.ensembl):{'$lt':motiff_start}
                },
                {'$push':{'{}.motif'.format(self.ensembl) :latter_dic}},
                upsert=False,
                multi=True
                )      
            n += 1

            print self.grch,'motif line',n

def main():

    modelhelp = model_help.replace('&'*6,'ENSEMBL_GENE').replace('#'*6,'ensembl_gene')

    funcs = (downloadData,extractData,updateData,selectData,dbMap,ensembl_gene_store)

    getOpts(modelhelp,funcs=funcs)
        
if __name__ == '__main__':
    main()