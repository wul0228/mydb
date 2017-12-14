#!/usr/bin/env python
# --coding:utf-8--
# date:2017/12/13
# author:wuling
# emai:ling.wu@myhealthgene.com

from template import *
import os,json,sys,getopt
from datetime import datetime
from pprint import pprint
from pymongo import MongoClient
from multiprocessing.dummy import Pool as ThreadPool

pjoin = os.path.join

psplit = os.path.split

listdir = os.listdir

today = datetime.now().strftime('%y%m%d%H%M%S')

# chose model to build  a topic-base
include_model = ['ncbi_gene','ensembl_gene','go_gene','protein_atlas']