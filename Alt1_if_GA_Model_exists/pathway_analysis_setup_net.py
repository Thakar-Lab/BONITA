#import python packages
import pkg_resources
#pkg_resources.require("networkx==1.11")
import networkx as nx
import operator
from sets import Set
import scipy.stats as stat
import requests
import argparse as argparse
from scipy.stats import variation
import numpy as np
import csv as csv
import pickle
# import other pieces of our software
import networkConstructor as nc
from utils import readFpkmData
import binascii

# download and prepare graph for finding the rules
def retrieveGraph(pathway, geneDict):
		#with open('WP23_edgeList_edited.txt_net_2.rules.graphml', 'r') as file:
		with open(str(pathway), 'r') as file:
    			data = file.read().replace('\n', '')
		graph=nx.parse_graphml(data)
		names=nx.get_node_attributes(graph, "gene_symbol")
		dicty1={key.encode('ascii'): value.encode('utf-8') for key, value in names.iteritems()}
		graph=nx.relabel_nodes(graph, dicty1, copy=True)
		#coder='WP23_edgeList_edited.txt_net_2.rules'
		coder = str(pathway)[:8]
		nx.write_gpickle(graph,coder+'.gpickle')
		# save the removed nodes and omics data values for just those nodes in the particular pathway
		pathwaySampleList=[{} for q in range(len(geneDict[list(graph.nodes())[0]]))]
		for noder in graph.nodes():
			for jn in range(len(pathwaySampleList)):
				pathwaySampleList[jn][noder]=geneDict[noder][jn]
		pickle.dump( pathwaySampleList, open( coder+"_sss.pickle", "wb" ) )
if __name__ == '__main__':

	# read in options
	parser = argparse.ArgumentParser(prog='BONITA')
	parser.set_defaults(verbose=False, mode='PA',sep=',')
	parser.add_argument("-v", action="store_true", dest="verbose",  help="output ongoing iterations to screen [default off]")
	parser.add_argument("-m", "--mode", metavar="mode", help="What BONITA functions should be run?")
	parser.add_argument("-sep", "--sep", metavar="seperator", help="How are columns in datafile specified")
	parser.add_argument("-t", action='store_const',const='\t', dest="sep",help="Tab delimited?")
	parser.add_argument("data")
	parser.add_argument("pathways") #'WP23_edgeList_edited.txt_net_2.rules.graphml' # 'filtered.c2.cp.kegg.v3.0.symbols.gmt'
	results = parser.parse_args()
	dataName=results.data
	pathwaysName=results.pathways
	verbose=results.verbose
	mode=results.mode


	#print(set(np.array([list(graph.node[n].keys()) for n in graph.nodes()]).flatten()))
	sss, geneDict, cvDict=readFpkmData(dataName, results.sep) # read in data
	pickle.dump( sss, open( 'sss.pickle', "wb" ) ) # save data in correct format for runs
	retrieveGraph(pathways, geneDict) # generate gpickles needed for pathway analysis