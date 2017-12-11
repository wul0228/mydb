#!/usr/bin/env python
# --coding:utf-8--
# date: 2017/12/05
# author:wuling
# emai:ling.wu@myhealthgene.com

#this model set  to download,extract,standard insert and select pathway data from wiki pathway

import sys
sys.path.append('../')
sys.setdefaultencoding = ('utf-8')
from share import *
from config import *  

__all__ = ['downloadData','extractData','standarData','insertData','updateData','selectData']

version  = 1.0

model_name = psplit(os.path.abspath(__file__))[1]

(wiki_pathway_model,wiki_pathway_raw,wiki_pathway_store,wiki_pathway_db,wiki_pathway_map) = buildSubDir('wiki_pathway')

log_path = pjoin(wiki_pathway_model,'wiki_pathway.log')

# main code
def downloadData(redownload = False):

    '''
    this function is to download the raw data from wiki web
    args:
    redownload-- default False, check to see if exists an old edition before download
                       -- if set to true, download directly with no check
    rawdir-- the directory to save download file
    '''
    if  not redownload:

        (choice,existwikiFile) = lookforExisted(wiki_pathway_raw,'pathway')

        if choice != 'y':
            return

    if redownload or not existwikiFile or  choice == 'y':

        process = wiki_parser(today)

        (down_url,mt) = process.getMt()

        # download file
        unzipdir = process.wget(down_url,mt,wiki_pathway_raw)

    # create log file
    if not os.path.exists(log_path):

        with open('./wiki_pathway.log','w') as wf:
            json.dump({
                'wiki_pathway':[(mt,today,model_name),]
                },wf,indent=2)

    print  'datadowload completed !'

    filepaths = [pjoin(unzipdir,filename) for filename in unzipdir]

    return (filepaths,today)

def extractData(filepaths,version):

    rawdirname = psplit(psplit(filepaths[0])[0])[1].strip()

    storedir = pjoin(wiki_pathway_store,rawdirname)

    createDir(storedir)

    process = wiki_parser(version)

    for filepath in filepaths:

        # if filepath.count('Hs_EBV_LMP1_signaling_WP262_86888.gpml'):

        #     process.gpml(filepath,storedir)


        process.gpml(filepath,storedir)

    print 'extract an insert completed!'

    return (storedir,version)

def updateData():

    wiki_pathway_log = json.load(open(log_path))

    rawdir = pjoin(wiki_pathway_raw,'pathway_update_{}'.format(today))

    latest = wiki_pathway_log.get('wiki_pathway')[-1][0].strip()

    process = wiki_parser(today)

    (dowload_url,mt) = process.getMt()

    if mt != latest:

        createDir(rawdir)

        filepaths,version  = downloadData(redownload = True)

        extractData(filepaths,version)

        wiki_pathway_log['wiki_pathway'].append((mt,today,model_name))

        with open(log_path,'w') as wf:

            json.dump(wiki_pathway_log,wf,indent=2)

        print  '{} \'s new edition is {} '.format('wiki_pathway',mt)

        bakeupCol('wiki_pathway_{}'.format(version),'wiki_pathway')

    else:

        print  '{} {} is the latest !'.format('wiki_pathway',mt)

def selectData(querykey = 'name',value='EBV LMP1 signaling'):
    '''
    this function is set to select data from mongodb
    args:
    querykey -- a specified field in database
    queryvalue -- a specified value for a specified field in database
    '''
    conn = MongoClient('127.0.0.1', 27017 )

    db = conn.mydb

    colnamehead = 'wiki_pathway'

    dataFromDB(db,colnamehead,querykey,queryvalue=None)


class dbMap(object):

    #class introduction

    def __init__(self,version):

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        col = db.get_collection('wiki_pathway_{}'.format(version))

        self.col = col

        self.version = version

    def mapXX2XX(self):

        docs =self.col.find({})

        for doc in docs:
            pass

    def mapping(self):

        self.mapXX2XX()

class wiki_parser(object):

    def __init__(self, version):

        self.version = version

        conn = MongoClient('localhost',27017)

        db = conn.get_database('mydb')

        col = db.get_collection('wiki_pathway_{}'.format(self.version))

        self.col = col

    def getMt(self):

        headers = {'User-Agent':'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Ubuntu Chromium/37.0.2062.120 Chrome/37.0.2062.120 Safari/537.36',}

        web = requests.get(wiki_pathway_download,headers = headers,verify=False)

        soup = bs(web.content,'lxml')

        down = soup.find(text = 'Current version: ')

        p =down.findParent('p')

        a = p.findChild('a')

        mt = a.text.split('(')[0].strip()

        href = a.attrs.get('href')

        dowload_url  = href + 'gpml/wikipathways-{}-gpml-Homo_sapiens.zip'.format(mt)
        
        mt = '213' + mt

        return  (dowload_url,mt)

    def wget(self,url,mt,rawdir):

        # download file
        filename = url.rsplit('/',1)[1].strip().rsplit('.',1)[0].strip()

        savename = '{}_{}_{}.zip'.format(filename,mt,today)

        storefilepath = pjoin(rawdir,savename)

        command = 'wget -O {} {}'.format(storefilepath,url)

        os.popen(command)

        # unzip file
        unzipdir = pjoin(wiki_pathway_raw,'pathway_{}_{}'.format(mt,today))

        createDir(unzipdir)

        unzip = 'unzip {} -d {}'.format(storefilepath,unzipdir)

        os.popen(unzip)

        # remove file
        remove = 'rm {}'.format(storefilepath)

        os.popen(remove)

        return unzipdir

    def gpml(self,filepath,storedir):
            '''
            this  function is to parse the wiki human pathway infos 
            1. get the name and organism
            2. get the catogory and description in comment
            3. get the datanode infos and classfied to GeneProduct,Complex,Metabolite and so on according to the type info
            4. get the relation between groupId and groupRef
            4. get the interaction  infos ,specially  the point infos
            '''
            print '~'*50

            filename = psplit(filepath)[1].strip().split('.gpml')[0].strip()

            print filename

            jsonfile = parse(open(filepath))

            # create a dict to store pathway
            pathwaydic = dict()

            # a file just contain a pathway
            pathway = jsonfile.get('Pathway')

            # q. get the name and organism
            pathwaydic.update(
                {
                'name':pathway.get('@Name'),
                'organism':pathway.get('@Organism')
                })

            #2. get the catogory and description in comment
            # sometimes  there is just one record of comment (dict), 
            # multi comment usally is dict but unicode included sometime
            comment = pathway.get('Comment')

            if comment:

                if  isinstance(comment,dict): 

                    Source = comment.get("@Source")

                    if Source in ['WikiPathways-category','WikiPathways-description']:

                        pathwaydic.update({Source:comment.get("#text")})

                elif isinstance(comment,list):
                  
                    for i in comment:

                        if isinstance(i,dict):

                            Source = i.get("@Source")

                            # get all category, one or more catagory infos
                            if Source == 'WikiPathways-category':

                                if  Source not in pathwaydic:

                                    pathwaydic.update({Source:[i.get("#text"),]})

                                else:

                                    pathwaydic[Source].append(i.get("#text"))

                            if Source =='WikiPathways-description':

                                # get the longest description,one or more description infos 
                                if  Source not in pathwaydic:

                                    pathwaydic.update({Source:i.get("#text")})

                                else:

                                    text = i.get("#text")

                                    # compare the length ,return the longest
                                    compare_len = lambda x,y : x if len(x)>=len(y) else y

                                    pathwaydic[Source] = compare_len(pathwaydic[Source],text)
                else:

                    pass

            #3. get the datanode infos and classfied to GeneProduct,Complex,Metabolite and so on according to the type info
            # a pathway must have more than one node
            datanode = pathway.get('DataNode')

            (nodedic,nodeinfo,groupId_nodeId) = self.gpml_datanode(datanode)

            pathwaydic.update(nodedic)

            # pathwaydic.update({'groupId_nodeId':groupId_nodeId})
            # get the relation between groupId with graphId
            group = pathway.get('Group')

            graphId_groupId = self.gpml_group(group)

            label = pathway.get('Label')

            graphId_label = self.gpml_label(label)

            # a pathway must have more than one interaction

            interaction = pathway.get('Interaction')

            anchorgraphId_pointgraphId = self.gpml_anchor(interaction)

            allgraphId_nodegraphId = self.gpml_allgraphId_nodegraphId(nodeinfo,groupId_nodeId,graphId_groupId,anchorgraphId_pointgraphId,graphId_label)

            if interaction:

                if isinstance(interaction,dict):interaction = [interaction,]

                inters = list()

                print 'len(interaction)',len(interaction)

                for inter in interaction:

                    points = inter.get('Graphics',{}).get('Point')

                    # points must have more than 2 point that have 2 GraphId
                    if points and len(points) >= 2:

                        adic = dict()

                        adic['ArrowHead'] = dict()

                        adic['ArrowHead']['output'] = list()
                        adic['ArrowHead']['input'] = list()

                        # n = 1 #entry number

                        m = 1 #group number

                        for p in points:

                            GraphId = p.get('@GraphRef')

                            if not GraphId:

                                continue

                            ArrowHead = p.get('@ArrowHead')

                            if  ArrowHead:

                                adic['ArrowHead']['type'] = ArrowHead

                            # look for nodeinfo
                            group_num = 'group{}'.format(m)

                            # adic[group_num] = dict()

                            entry = allgraphId_nodegraphId.get(GraphId)

                            # print '*****',GraphId,entry

                            if entry and isinstance(entry,(unicode,list)):

                               # adic[group_num]['entry{}'.format(n)] = entry 
                               adic[group_num] = entry 
                                
                            elif entry and  isinstance(entry,dict):

                                adic[group_num]= entry.values()

                            else:
                                pass
                            
                            if ArrowHead:
                                adic['ArrowHead']['output'] .append(group_num)
                            else:
                                adic['ArrowHead']['input'].append(group_num)

                            m += 1

                            if  not ArrowHead:

                                adic['ArrowHead']['type'] = 'line'

                            del ArrowHead

                            # modify output 

                        if len(adic)  >=3:

                            adic = self.gpml_arrowhead(adic)

                            inters.append(adic)

                        else:
                            pass

                    if inters:

                         pathwaydic.update({'Interaction':inters})

                with open('./test.json','w') as wf:
                    json.dump(pathwaydic,wf,indent=8)

                print 'len(pathwaydic["Interaction"])',len(pathwaydic.get("Interaction",{}))

            with open(pjoin(storedir,'{}.json'.format(filename)),'w' ) as wf:
                json.dump(pathwaydic,wf,indent=2)

            self.col.insert(pathwaydic)

    def gpml_datanode(self,datanode):

        # create a dict to store node info ,classfied to node type class
        nodedic = dict()

        # create a dict to store node info, with GraphId as the key
        nodeinfo = dict()

        # create a dict  to store groupref2nodeids
        groupId_nodeId = dict()

        for node in datanode:
            #----------------node graphid--------------------
            node_graphid = node.get('@GraphId')

            nodeinfo[node_graphid] = dict()

            #----------------node groupref--------------------
            node_groupref = node.get('@GroupRef')

            if node_groupref not in groupId_nodeId:

                groupId_nodeId[node_groupref] = list()

            groupId_nodeId[node_groupref].append(node_graphid)

            #----------------node type--------------------
            node_type = node.get('@Type')
            # some node hava no type
            if not node_type:

                node_type = 'Other'

            nodeinfo[node_graphid]['type'] = node_type

            # gene metabbolit complex pathway geneproduct protein
            if node_type not in nodedic:
                
                nodedic[node_type] = dict()

            node_lable =  node.get('@TextLabel').replace(' ','&').replace('.','*').strip()

            nodeinfo[node_graphid]['name'] = node_lable

            if node_lable not in nodedic[node_type]:

                nodedic[node_type][node_lable] = dict()

            node_comment = node.get('Comment')

            if node_comment:

                if isinstance(node_comment,unicode):
                    nodedic[node_type][node_lable]['comment'] = node_comment
                    nodeinfo[node_graphid]['comment'] = node_comment

                elif isinstance(node_comment,dict):
                    node_comment = [node_comment,]

                if isinstance(node_comment,list):

                    nodedic[node_type][node_lable]['comment'] = dict()
                    nodeinfo[node_graphid]['comment'] = dict()
                    for c in node_comment:
                        if isinstance(c,dict):
                            nodedic[node_type][node_lable]['comment'].update({c.get('@Source'):c.get("#text")})
                            nodeinfo[node_graphid]['comment'].update({c.get('@Source'):c.get("#text")})
                        else:
                             nodeinfo[node_graphid]['comment'].update({'other':c})

            node_xref = node.get("Xref")

            if node_xref:

                if isinstance(node_xref,dict): node_xref=[node_xref,]

                nodedic[node_type][node_lable]['xref'] = dict()
                
                nodeinfo[node_graphid]['xref'] = dict()

                for xref in node_xref:
                    db = xref.get('@Database')
                    _id = xref.get('@ID')

                    if db:
                        nodedic[node_type][node_lable]['xref'].update({db:_id})

                        nodeinfo[node_graphid]['xref'].update({db:_id})
        
        # with open('./groupId_nodeId.json','w') as wf:
        #     json.dump(groupId_nodeId,wf,indent=2)

        # with open('./nodeinfo.json','w') as wf:
        #     json.dump(nodeinfo,wf,indent=2)
            
        return (nodedic,nodeinfo,groupId_nodeId)

    def gpml_group(self,group):

        graphId_groupId = dict()

        if group:

            if isinstance(group,dict): group=[group,]

            for g in group:

                group_groupId =g.get('@GroupId')
                group_graphId = g.get('@GraphId')

                # a graphId to 2 or more group_groupId??
                graphId_groupId[group_graphId] = group_groupId

        # with open('./graphId_groupId.json','w') as wf:
        #     json.dump(graphId_groupId,wf,indent=2)

        return graphId_groupId

    def gpml_label(self,label):

        graphId_label = dict()

        if label:

            if isinstance(label,dict): label=[label,]

            for l in label:
                label_graphId =l.get('@GraphId')
                label_name = l.get('@TextLabel').strip()

                graphId_label[label_graphId] = label_name

        return graphId_label

    def gpml_anchor(self,interaction):

        if interaction:

            if isinstance(interaction,dict):interaction = [interaction,]

            anchorgraphId_pointgraphId = dict()

            for index,inter in enumerate(interaction):

                anchors = inter.get('Graphics',{}).get('Anchor')

                if anchors:

                    if isinstance(anchors,dict):anchors = [anchors,]

                    for anchor in anchors:

                        anchor_graphId = anchor.get('@GraphId')

                        anchorgraphId_pointgraphId[anchor_graphId] = list()

                        points = inter.get('Graphics',{}).get('Point')

                        # if points and len(points) >= 2:
                        if points:
                       
                            for p in points:

                                GraphId = p.get('@GraphRef')

                                if not GraphId:

                                    continue

                                anchorgraphId_pointgraphId[anchor_graphId].append(GraphId)
    
            return anchorgraphId_pointgraphId

    def gpml_arrowhead(self,adic):

        ArrowHead = adic.get('ArrowHead')

        _input = ArrowHead.get('input')

        _output = ArrowHead.get('output')

        arrow_type = ArrowHead.get('type')

        if not arrow_type:

            print adic
        if (_input and not _output) or  arrow_type.lower() == 'line':

            effect_type = '   [ -------- ]   '.join(_input)

        elif _output and not _input:

            effect_type = '   [ <------> ]   '.join(_output)

        elif _input and _output:

            if arrow_type.strip() == 'Arrow':
                arrow_type = " --------> " 

            if arrow_type.strip() == 'TBar': 
                arrow_type = ' --------| '

            effect_type = ','.join(_input)   + '   [ ' + arrow_type + ' ]   ' + ','.join(_output)
            # print effect_type
        else:
            pass

        adic.pop('ArrowHead')

        adic.update({'effect_type':effect_type})

        return adic

    def gpml_allgraphId_nodegraphId(self,nodeinfo,groupId_nodeId,graphId_groupId,anchorgraphId_pointgraphId,graphId_label):

        # print len(nodeinfo)
        # print len(graphId_groupId)
        # print len(anchorgraphId_pointgraphId)
        # print len(graphId_label)

        allgraphId_nodegraphId = dict()

        for node,val in nodeinfo.items():

            allgraphId_nodegraphId.update({node:[val,]})

        for graphId,groupId in graphId_groupId.items():

            allgraphId_nodegraphId[graphId] = list()

            group_nodes = groupId_nodeId.get(groupId)

            if group_nodes:

                for node_id in group_nodes:

                    node_info = nodeinfo.get(node_id)

                    if node_info:

                        allgraphId_nodegraphId[graphId].append(node_info)

        if anchorgraphId_pointgraphId:

            for graphId,graphIds in anchorgraphId_pointgraphId.items():

                allgraphId_nodegraphId[graphId] = dict()

                # print graphId

                for _gid in graphIds:

                    _gid_node = allgraphId_nodegraphId.get(_gid)

                    if _gid_node:

                        allgraphId_nodegraphId[graphId].update({
                            _gid:_gid_node
                            })

        allgraphId_nodegraphId.update(graphId_label)

        # with open('allgraphId_nodegraphId.json','w') as wf:
# 
        #     json.dump(allgraphId_nodegraphId,wf,indent=2)

        # print len(allgraphId_nodegraphId)

        return allgraphId_nodegraphId

def main():

    modelhelp = model_help.replace('&'*6,'WIKI_PATHWAY').replace('#'*6,'wiki_pathway')

    funcs = (downloadData,extractData,updateData,selectData,dbMap,wiki_pathway_store)

    getOpts(modelhelp,funcs=funcs)
        
if __name__ == '__main__':

    main()

    # rawdir = '/home/user/project/dbproject/mydb_v1/wiki_pathway/dataraw/pathway_21320171116_171205101037'

    # filepaths = [pjoin(rawdir,filename) for filename in listdir(rawdir)]

    # extractData(filepaths,'171205101037')
    # updateData()