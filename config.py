

#!/usr/bin/env python
# ---coding:utf-8---
# date:20171123
# author:wuling
# emai:ling.wu@myhealthgene.com

# this model is set to
#+++++++++++++++++++++++++ packages ++++++++++++++++++++++++++++++++++++++#
import os , sys, getopt, tsv, json, time,copy,requests ,pprint,xlrd
from template import *
from ftplib import FTP
from time import sleep
from lxml import etree as et
from xmltodict import parse
from datetime import datetime
from selenium import webdriver
from pymongo import MongoClient
from bs4 import BeautifulSoup as bs
from multiprocessing.dummy import Pool as ThreadPool

#+++++++++++++++++++++++++ simplify  method+++++++++++++++++++++++++++++++++#
listdir = os.listdir

pjoin = os.path.join

psplit = os.path.split

pexists = os.path.exists
#++++++++++++++++++++++++universal consttant +++++++++++++++++++++++++++++++++#

mydb_path =psplit(os.path.abspath(__file__))[0]

now  = datetime.now().strftime('%y%m%d')

today = datetime.now().strftime('%y%m%d%H%M%S')

#~~~~~~~~~~~~~~~~~~~NCBI gene~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ncbi_gene_ftp_infos = {
    'host' : 'ftp.ncbi.nlm.nih.gov' ,
    'user':'anonymous',
    'passwd' : '',
    'logdir' : '/gene/DATA/'
    }
ncbi_gene_data_ftp_path = '/gene/DATA/'

ncbi_gene_expression_path = '/gene/DATA/expression/Mammalia/Homo_sapiens/'

ncbi_gene_filenames = ['gene_group.gz','gene_neighbors.gz','gene2pubmed.gz','gene_info.gz']

#~~~~~~~~~~~~~~~~~~~Ensemble gene~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ensembl_gene_ftp_infos = {
'host' : 'ftp.ensembl.org' ,
'user':'anonymous',
'passwd' : '',
'logdir' : '/pub/'
}

ensembl_gtfGRch38_ftp_path = '/pub/current_gtf/homo_sapiens/'

filename_gtfGRch38 = 'Homo_sapiens.GRCh38.90.chr.gtf.gz'

ensembl_regulatorGRch38_ftp_path ='/pub/current_regulation/homo_sapiens/'

filenames_regulatorGRch38 = [
'homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20161111.gff.gz',
'homo_sapiens.GRCh38.motiffeatures.20161111.gff.gz'
]

ensembl_gtfGRch37_ftp_path = '/pub/grch37/current/gtf/homo_sapiens/'

filename_gtfGRch37 = 'Homo_sapiens.GRCh37.87.chr.gtf.gz'

ensembl_regulatorGRch37_ftp_path = '/pub/grch37/current/regulation/homo_sapiens/'

filenames_regulatorGRch37 = [
'homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20161117.gff.gz',
'homo_sapiens.GRCh37.motiffeatures.20161117.gff.gz'
]

ensembl_file_mark = {
'GRCh38':'chr.gtf',
'GRCh37':'chr.gtf',
'GRCh38_regularory':'regulatory_features',
'GRCh37_regularory':'regulatory_features',
 'GRCh38_motif':'motiffeatures',
 'GRCh37_motif':'motiffeatures',
 }

ensembl_file_ftplogdir={
    'GRCh38':ensembl_gtfGRch38_ftp_path,
    'GRCh38_regularory':ensembl_regulatorGRch38_ftp_path,
    'GRCh38_motif':ensembl_regulatorGRch38_ftp_path,
    'GRCh37':ensembl_gtfGRch37_ftp_path,
    'GRCh37_regularory':ensembl_regulatorGRch37_ftp_path,
    'GRCh37_motif':ensembl_regulatorGRch37_ftp_path
    }
#~~~~~~~~~~~~~~~~~~~GO gene~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
go_gene_ftp_infos =  {
'host' : 'ftp.ebi.ac.uk' ,
'user':'anonymous',
'passwd' : '',
'logdir' : '/pub/databases/GO/goa/HUMAN/'
}


go_web = 'http://www.ebi.ac.uk/QuickGO/'

go_gene_filenames = [
'goa_human.gpa.gz',
# 'goa_human.gpi.gz'
]

go_obo_ftp_infos =  {
'host' : 'ftp.geneontology.org' ,
'user':'anonymous',
'passwd' : '',
'logdir' : '/pub/go/ontology/'
}

go_obo_filenames=[
'go.obo',
]

go_gene_obo_link = 'http://purl.obolibrary.org/obo/go.obo'

#~~~~~~~~~~~~~~~~~~~hgnc  gene~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
hgnc_gene_ftp_infos =  {
'host' : 'ftp.ebi.ac.uk' ,
'user':'anonymous',
'passwd' : '',
'logdir' : '/pub/databases/genenames/new/tsv/'
}

hgnc_genename_filename = 'hgnc_complete_set.txt'
#~~~~~~~~~~~~~~~~~~~kegg pathway~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# kegg_pathway_downloadurl = 'http://www.kegg.jp/kegg-bin/get_htext?hsa00001'
kegg_pathway_downloadurl = 'http://www.kegg.jp/kegg-bin/download_htext?htext=hsa00001&format=htext&filedir='
kegg_pathway_web = 'http://www.kegg.jp/kegg-bin/get_htext?ko00000.keg'

#~~~~~~~~~~~~~~~~~~~reactom pathway~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# reactome_download_web1 = 'https://reactome.org/download/current/'
# reactome_download_web2 = 'https://reactome.org/download/current/interactors/'

# reactom_pathway_urls ={
# # protein-protein interactor
# 'reactome.all_species.interactions.tab-delimited':'https://reactome.org/download/current/interactors/reactome.all_species.interactions.tab-delimited.txt',
# # gene pathway
# 'UniProt2Reactome':'https://reactome.org/download/current/UniProt2Reactome.txt',
# # 'https://reactome.org/download/current/NCBI2Reactome.txt',
# # 'https://reactome.org/download/current/ChEBI2Reactome.txt',
# # 'https://reactome.org/download/current/Ensembl2Reactome.txt',
# # 'https://reactome.org/download/current/miRBase2Reactome.txt',
# 'ReactomePathwaysRelation':'https://reactome.org/download/current/ReactomePathwaysRelation.txt',
# }

# reactomfile_mt = {
# 'UniProt2Reactome':'',
# 'ReactomePathwaysRelation':'',
# 'reactome.all_species.interactions.tab-delimited':''
# }
reactome_download_web1 = 'https://reactome.org/download/current/diagram/'
reactome_download_web2 = 'https://reactome.org/download/current/'
#~~~~~~~~~~~~~~~~~~~wiki pathway~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

wiki_pathway_download = 'https://www.wikipathways.org/index.php/Download_Pathways'

#~~~~~~~~~~~~~~~~~~~clinvar varient~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
clinVar_varient_ftp_infos = {
'host' : 'ftp.ncbi.nlm.nih.gov' ,
'user':'anonymous',
'passwd' : '',
'logdir' : '/pub/clinvar/tab_delimited/'
}

clinVar_varient_filename = 'variant_summary.txt.gz'

#~~~~~~~~~~~~~~~~~~~disgenet~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

disgenet_download_web = 'http://www.disgenet.org/web/DisGeNET/menu/downloads#gdascurated'

disgenet_download_url = 'http://www.disgenet.org/ds/DisGeNET/results/all_gene_disease_pmid_associations.tsv.gz'

#~~~~~~~~~~~~~~~~~~~human protein atlas~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

protein_atlas_download_web = 'https://www.proteinatlas.org/about/download'
#~~~~~~~~~~~~~~~~~~~miRTarbase ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


miRTarbase_homepage = 'http://mirtarbase.mbc.nctu.edu.tw/'
miRTarbase_download_web = 'http://mirtarbase.mbc.nctu.edu.tw/php/download.php'


#~~~~~~~~~~~~~~~~~~~hpo phenotypic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# hpo_web = 'http://compbio.charite.de/hpoweb/showterm?id=HP:0000118'
hpo_download_web = 'http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/'

hpo_download_urls = [
'http://purl.obolibrary.org/obo/hp.obo',
'http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/OMIM_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt',
'http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ORPHA_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt',
'http://compbio.charite.de/jenkins/job/hpo.annotations/lastSuccessfulBuild/artifact/misc/phenotype_annotation.tab'
]