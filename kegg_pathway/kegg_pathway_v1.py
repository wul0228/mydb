#!/usr/bin/env python
# --coding:utf-8--
# date: 2017/11/30
# author:ling.wu
# emai:ling.wu@myhealthpathway.com

#this model set  to download,extract,standard insert and select pathway data from kegg pathway

import sys
sys.path.append('../')
sys.setdefaultencoding = ('utf-8')
from share import *
from config import *  


__all__ = ['downloadData','extractData','standarData','insertData','updateData','selectData']

version  = 1.0

model_name = psplit(os.path.abspath(__file__))[1]

(kegg_pathway_model,kegg_pathway_raw,kegg_pathway_store,kegg_pathway_db,kegg_pathway_map) = buildSubDir('kegg_pathway')

log_path = pjoin(kegg_pathway_model,'kegg_pathway.log')

# main code
def downloadData(redownload = False,rawdir = None):

    '''
    this function is to download the raw data from kegg web
    args:
    redownload-- default False, check to see if exists an old edition before download
                       -- if set to true, download directly with no check
    rawdir-- the directory to save download file
    '''
    if  not redownload:

        (choice,existkeggFile) = lookforExisted(kegg_pathway_raw,'pathway')

        if choice != 'y':
            return

    if redownload or not existkeggFile or  choice == 'y':

        if not rawdir:

            rawdir = pjoin(kegg_pathway_raw,'pathway_{}'.format(today))

        process = kegg_parser(today)

        mt = process.getMt()

        # download pathway file
        options = webdriver.ChromeOptions()

        prefs = {'profile.default_content_settings.popups':0,'download.default_directory':rawdir}

        options.add_experimental_option('prefs',prefs)

        driver = webdriver.Chrome(chrome_options=options)

        driver.get('http://www.kegg.jp/kegg-bin/get_htext?hsa00001')

        # download = driver.find_element_by_link_text('Download htext')

        download = driver.find_element_by_link_text('Download json')

        download.click()

        sleep(5)

        driver.close()

    if not os.path.exists(log_path):

        with open(log_path,'w') as wf:
            json.dump({'hsa00001.json':[(mt,today,model_name),]},wf,indent=2)

    print  'datadowload completed !'
    
    filepath = pjoin(rawdir,'hsa00001.json')

    sleep(5)

    return (filepath,today)

def extractData(filepath,version):

    pathway_raw = psplit(filepath)[0].strip()

    pathway_rawdirname = psplit(pathway_raw)[1].strip()

    pathway_store = pjoin (kegg_pathway_store,pathway_rawdirname)
    
    createDir(pathway_store)

    process = kegg_parser(version)

    all_path = process.pathway(filepath,pathway_store)

    print 'all_path',len(all_path)

    # get all  (XML FILE)
    func = lambda x:process.pathway_relation(x,pathway_raw)

    path_ids = [i.get("path_id") for i in all_path]

    print 'path ids :',len(path_ids)

    multiProcess(func,path_ids,size=50)

    print 'get xml file completed'

    # extract entry  reaction and relations data from xml file
    relation_filepaths = [pjoin(pathway_raw,filename) for filename in listdir(pathway_raw) if filename.endswith('.xml') ]

    func = lambda x:process.pathway_standar(x,pathway_store)

    multiProcess(func,relation_filepaths,size = 20)

    print 'relation extract and insert completed'

    return(pathway_store,version)

def updateData(insert=False):

    kegg_pathway_log = json.load(open(log_path))

    rawdir = pjoin(kegg_pathway_raw,'pathway_update_{}'.format(today))

    process = kegg_parser(today)

    mt = process.getMt()

    if mt != kegg_pathway_log.get('hsa00001.json')[-1][0]:

        createDir(rawdir)

        (filepath,version)= downloadData(redownload=True,rawdir=rawdir)

        extractData(filepath,version)

        _map = dbMap(version)

        _map.mapping()

        kegg_pathway_log['hsa00001.json'].append((mt,today,model_name))

        with open(log_path,'w') as wf:
            json.dump(kegg_pathway_log,wf,indent=2)

        print  '{} \'s new edition is {} '.format('kegg_pathway',mt)

        bakeupCol('kegg_pathway_{}'.format(version),'kegg_pathway')
        
    else:
        print  '{} {} is the latest !'.format('kegg_pathway',mt)

def selectData(querykey = 'path_id',value='00010'):
    '''
    this function is set to select data from mongodb
    args:
    querykey -- a specified field in database
    queryvalue -- a specified value for a specified field in database
    '''
    conn = MongoClient('127.0.0.1', 27017 )

    db = conn.mydb

    colnamehead = 'kegg_pathway'

    dataFromDB(db,colnamehead,querykey,queryvalue=None)

class kegg_parser(object):

    def __init__(self,version):

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        colname = 'kegg_pathway_{}'.format(version)

        col = db.get_collection(colname)

        self.db = db

        self.col = col

        self.version = version

        self.colname = colname

    def getMt(self):

        headers = {'User_Agent':'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Ubuntu Chromium/37.0.2062.120 Chrome/37.0.2062.120 Safari/537.36'}

        web = requests.get(kegg_pathway_web,headers=headers,verify=False)

        soup = bs(web.content,'lxml')

        mt = soup.select('body')[0].text.split('Last updated: ')[1].strip().replace(' ','&').replace(',','#')

        return mt

    def pathway(self,filepath,pathway_store):

        jsonfile = json.load(open(filepath))

        filname = jsonfile.get('name')

        childrens = jsonfile.get('children')

        all_path = list()

        n = 0
        for path_class in childrens:
            
            path_class_name = path_class.get('name')
            path_class_children = path_class.get('children')

            # print path_class_name

            for path_subclass in  path_class_children:
                path_subclass_name = path_subclass.get('name')
                path_subclass_children = path_subclass.get('children')

                # print '             ',path_subclass_name

                n  += len(path_subclass_children)

                for path in path_subclass_children:
                    
                    path_name_info = path.get('name')
                    path_gene_info = path.get('children')

                    # print '                             ',path_name_info

                    path_id = path_name_info.split(' ',1)[0].strip()
                    path_name = path_name_info.rsplit('[PATH')[0].replace(path_id,'').strip()

                    path = {
                    'path_id':path_id,
                    'path_name':path_name,
                    'path_class':path_class_name,
                    'path_subclass':path_subclass_name
                    }

                    if path_gene_info:

                        genes = dict()

                        for gene in path_gene_info:

                            gene_id,gene_symbol,gene_name,ko_entry,ko_entry,definition= ('','','','','','')

                            gene_name_info = gene.get('name')

                            ko_head= gene_name_info.split('\t')[0]

                            if gene_name_info.count('\t'):
                                ko_tail= gene_name_info.split('\t')[1]
                            else:
                                ko_tail = ''
                            
                            if ko_head.count(';'):
                                (gene_id,gene_symbol) = tuple(ko_head.split(';',1)[0].strip().split(' ',1))
                                gene_name = ko_head.split(';',1)[1].strip()

                            else:
                                gene_id = ko_head.split(' ',1)[0]
                                gene_name = ko_head.split(' ',1)[1].strip()

                            gene={
                            'gene_symbol':gene_symbol,
                            'gene_name':gene_name,
                            }
                            if ko_tail:
                                (ko_entry,ko_name) = tuple(ko_tail.split(';',1)[0].strip().split(' ',1))
                                definition = ko_tail.split(';',1)[1].strip()
                                gene.update({
                                    'ko_entry':ko_entry,
                                    'ko_name':ko_name,
                                    'definition':definition
                                    })
                            else:
                                print 'no ko ',gene_name_info

                            genes[gene_id] = gene

                            path.update({'gene':genes})

                    all_path.append(path)

        print 'pathways:' ,n

        with open(pjoin(pathway_store,'hsa00001.json'),'w') as wf:
            json.dump(all_path,wf,indent=2)

        self.col.insert(all_path)

        print 'pathway extract and insert completed'
        
        
        return all_path

    def pathway_relation(self,path_id,rawdir):

        # download relation xml file
        url = 'http://rest.kegg.jp/get/hsa{}/kgml'.format(path_id)

        savefile_path = pjoin(rawdir,'relation_hsa{}.xml'.format(path_id))

        wget = 'wget --retry-connrefused  -O {} {}'.format(savefile_path,url)

        os.popen(wget)

    def pathway_standar(self,savefile_path,storedir):

        path_id = psplit(savefile_path)[1].split('.xml')[0].split('_hsa')[1].strip()

        # parser relation in xml file
        xmlfile = open(savefile_path)

        try:
            dictfile = parse(xmlfile)
        except:
            # print path_id,'relation xml file is blank'
            return
  
        pathway_info  = dictfile.get('pathway')

        kgml = dict()

        entry = pathway_info.get('entry')
        # entry_keys = reduce(lambda x,y: set(x) | set(y) ,[e.keys() for e in entry])

        if entry:

            entrydict= dict()

            for e in entry:

                entrydict.update(
                    {
                        e.get('@id'):{
                            'entry_name':e.get('@name'),
                            'entry_type':e.get('@type'),
                            'entry_reaction':e.get('@reaction'),
                        }
                    }
                )

            kgml['path_entry'] = entrydict

        reaction = pathway_info.get('reaction')
        # reaction_keys = reduce(lambda x,y: set(x) | set(y) ,[e.keys() for e in reaction])
        if reaction:

            reactiondict = dict()

            for rxn in reaction:

                reactiondict.update(
                    {
                        rxn.get('@id'):{
                            'reaction_name':rxn.get('@name'),
                            'reaction_type':rxn.get('@type'),
                            'reaction_substrate':{
                                'entry_id':rxn.get('@substrate',{}).get('@id'),
                                'entry_name':rxn.get('@substrate',{}).get('@name')
                                },
                            'reaction_product':{
                                'entry_id':rxn.get('@product',{}).get('@id'),
                                'entry_name':rxn.get('@product',{}).get('@name')
                                }
                        }
                    }
                )

            kgml['path_reaction'] = reactiondict

        relation = pathway_info.get('relation')

        if relation:

            # relation_keys = reduce(lambda x,y: set(x) | set(y) ,[e.keys() for e in relation])

            relationdict = dict()

            # just one relation 
            if isinstance(relation,dict):
                relation = [relation,]

            for re in relation:

                entry1 = re.get('@entry1')
                entry2 = re.get('@entry2')
                relation_type = re.get('@type')
                relation_subtype = re.get('subtype',{})

                entry1_name = entrydict.get(entry1).get('entry_name')
                entry2_name = entrydict.get(entry2).get('entry_name')

                relation_key = entry1_name + '_' + entry2_name

                # subtype  may be have more than one 
                relation_subtype_value = list()

                if isinstance(relation_subtype,dict):

                    relation_subtype_value = [{

                        'subtype_name':relation_subtype.get('@name'),

                        'subtype_value':relation_subtype.get('@value'),

                        },]

                elif isinstance(relation_subtype,list):

                    for subtype in relation_subtype:

                        relation_subtype_value.append({

                        'subtype_name':subtype.get('@name'),

                        'subtype_value':subtype.get('@value'),
                        }
                    )

                if relation_type not in relationdict:

                    relationdict[relation_type] = list()

                relationdict[relation_type].append({
                    'entry1':{
                        'id':entry1,
                        'name':entry1_name,
                        'type':entrydict.get(entry1).get('entry_type')
                        },
                    'entry2':{
                        'id':entry2,
                        'name':entry2_name,
                        'type':entrydict.get(entry2).get('entry_type')
                        },
                    'relation_subtype':relation_subtype_value
                    }
                )

            kgml['path_relation'] = relationdict

        with open(pjoin(storedir,'kgml_{}.json'.format(path_id)),'w') as wf:
            json.dump(kgml,wf,indent=2)

        self.col.update({'path_id':path_id},{'$set':kgml})

class dbMap(object):

    #class introduction

    def __init__(self,version):

        self.version = version

        conn = MongoClient('127.0.0.1',27017)

        db = conn.get_database('mydb')

        col = db.get_collection('kegg_pathway_{}'.format(self.version))

        docs = col.find({})

        self.docs  = docs

    def mapID2Name(self):

        id_name = dict()

        name_id = dict()

        for doc in self.docs:

            path_id = doc.get('path_id')

            path_name = doc.get('path_name')

            id_name.update({path_id:path_name})

            if path_name not in name_id:
                
                name_id[path_name] = list()

            name_id[path_name].append(path_id)

        with open(pjoin(kegg_pathway_map,'id2name_{}.json'.format(self.version)),'w') as wf:
            json.dump(id_name,wf,indent=2)

        with open(pjoin(kegg_pathway_map,'name2id_{}.json'.format(self.version)),'w') as wf:
            json.dump(name_id,wf,indent=2)

    def mapping(self):

        self.mapID2Name()

def main():

    modelhelp = model_help.replace('&'*6,'KEGG_PATHWAY').replace('#'*6,'kegg_pathway')

    funcs = (downloadData,extractData,updateData,selectData,dbMap,kegg_pathway_store)

    getOpts(modelhelp,funcs=funcs)
        
if __name__ == '__main__':

    main()
