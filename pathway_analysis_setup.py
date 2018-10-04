#import python packages
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
import utils as utils

# read in fpkm data in a csv and construct correct data formats
def readFpkmData(dataName, delmited):
	with open(dataName) as csvfile:
		data=[]
		reader = csv.reader(csvfile, delimiter=delmited)
		for row in reader:
			data.append(row)
	sampleList=[]
	geneDict={}
	cvDict={}
	for j in range(1,len(data[1])):
		sampleList.append({})
	for i in range(1,len(data)):
		tempDatalist=[]
		maxdata=0
		for j in range(1,len(data[i])):
			currentDataPoint=float(data[i][j])
			tempDatalist.append(currentDataPoint)
		maxdata=max(tempDatalist)
		cvDict[data[i][0]]=variation(tempDatalist)
		if maxdata==0:
			maxdata=1.
		geneDict[data[i][0]]=[itemer/maxdata for itemer in tempDatalist]
		for j in range(0,len(data[i])-1):
			sampleList[j][str.upper(data[i][0])]=float(data[i][1])/maxdata
	return sampleList, geneDict, cvDict

def read_gmt(filename):
	gmt_dict={}
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		newline=line.split('\t')
		gmt_dict[newline[0]]=Set(newline[2:])
	return gmt_dict

# find list of pathways with at least four genes found in data
def find_overlaps(filename,geneDict):
	overlapsets=[] # list of pathways with enough overlaps
	genes=Set(geneDict.keys())
	keggDict=read_gmt(filename)
	for key in keggDict.keys():
		if len(genes.intersection(keggDict[key]))>4: # ensure there are at least 4 nodes in both pathway and detected genes
			overlapsets.append(key)
			print(key)
			print(len(genes.intersection(keggDict[key])))
	return overlapsets

# download and prepare graph for finding the rules
def retrieveGraph(name,aliasDict,dict1,dict2, cvDict, geneDict):
	print(name)
	
	# use KEGG API to figure out what the pathway code is
	namelist=name.split('_')
	namelist.pop(0)
	requester='http://rest.kegg.jp/find/pathway/'+namelist.pop(0)
	for item in namelist:
		requester=requester+'+'+item
	r=requests.get(requester)
	genes=Set(geneDict.keys())
	lines=r.text
	# parse file that comes from KEGG
	if len(lines.split('\n')[0].split(':'))>1:
		code=lines.split('\n')[0].split(':')[1][3:8] # KEGG number of overlapped pathway
		graph=nx.DiGraph()
		# download and integrate human and generic versions of pathway
		coder=str('ko'+code) 
		nc.uploadKEGGcodes([coder], graph, dict2)
		coder=str('hsa'+code)
		nc.uploadKEGGcodes_hsa([coder], graph,dict1, dict2)
		# check to see if there is a connected component, simplify graph and print if so
		if len(list(nx.connected_component_subgraphs(graph.to_undirected() )))>0:
			newgraph = max(nx.connected_component_subgraphs(graph.to_undirected()), key=len)
			newOverlap=genes.intersection(set(newgraph.nodes()))
			graph, addLaterNodes=simplifyNetworkpathwayAnalysis(graph, cvDict)
			print('nodes: ',str(len(graph.nodes())),',   edges:',str(len(graph.edges())))
			if len(newOverlap)>4: # save the graph if there is a connected component of at least length 4
				nx.write_graphml(graph,coder+'.graphml')
				nx.write_gpickle(graph,coder+'.gpickle')
				# save the removed nodes and omics data values for just those nodes in the particular pathway
				pathwaySampleList=[{} for q in range(len(geneDict[list(graph.nodes())[0]]))]
				for noder in graph.nodes():
					for jn in range(len(pathwaySampleList)):
						pathwaySampleList[jn][noder]=geneDict[noder][jn]
				pickle.dump( pathwaySampleList, open( coder+"_sss.pickle", "wb" ) )
				pickle.dump( addLaterNodes, open( coder+"_addLaterNodes.pickle", "wb" ) )
	else:
		print('not found:')
		print(requester)
		print(lines)

# identify pathways and complete setup for simulation
def findPathways(cvDict,gmtName, geneDict):
	aliasDict, dict1, dict2={}, {}, {} # set up dicts for reading KEGG files
	# read in kegg gene symbol dictionaries
	nc.parseKEGGdicthsa('inputData/hsa00001.keg',aliasDict,dict1)
	nc.parseKEGGdict('inputData/ko00001.keg',aliasDict,dict2)
	namelist=find_overlaps(gmtName,cvDict) # find list of pathways with overlaps with the genes from omics data
	print('num of overlap nodes: ' + str(len(namelist)))
	for name in namelist:
		retrieveGraph(name,aliasDict,dict1,dict2, cvDict, geneDict) # find and store gpickles for graphs found

# collapse unnecessary nodes for easier rule determination
def simplifyNetworkpathwayAnalysis(graph, ss):
	#network simplification algorithm. 
	# # 1. remove nodes with no input data
	# # 2. remove edges to nodes from complexes they are a part of 
	# # 4. remove nodes with only 1 input
	# # 5. remove nodes with only 1 output
	# # 6. remove self edges

	addBackNodes=[]

	# 1. remove self edges
	for edge in graph.edges():
		if edge[0]==edge[1]:
			graph.remove_edge(edge[0],edge[1])

	# 1. remove nodes with no input data
	removeNodeList= [x for x in graph.nodes() if  not x  in ss.keys()]
	for rm in removeNodeList:
		for start in graph.predecessors(rm):
			for finish in graph.successors(rm):
				edge1=graph.get_edge_data(start,rm)['signal']
				edge2=graph.get_edge_data(rm,finish)['signal']
				inhCount=0
				if edge1=='i':
					inhCount=inhCount+1
				if edge2=='i':
					inhCount=inhCount+1
				if inhCount==1:
					graph.add_edge(start,finish,signal='i')
				else:
					graph.add_edge(start,finish,signal='a')
		graph.remove_node(rm)


	#print(graph.nodes())

	# 2. remove dependence of nodes on complexes that include that node
	for node in graph.nodes():
		predlist=graph.predecessors(node)
		for pred in predlist:
			if '-' in pred:
				genes=pred.split('-')
				flag=True
				for gene in genes:
					if not gene in predlist:
						flag=False
				if flag:
					graph.remove_edge(pred,node)
		

	# # 4. rewire nodes that have only one upstream node
	# #print(len(graph.nodes()))
	# removeNodeList= [x for x in graph.nodes() if (len(graph.predecessors(x))==1) ]
	# for rm in removeNodeList:
	# 	before=graph.predecessors(rm)[0]
	# 	for after in graph.successors(rm):
	# 		edge1=graph.get_edge_data(before,rm)['signal']
	# 		edge2=graph.get_edge_data(rm,after)['signal']
	# 		inhCount=0
	# 		if edge1=='i':
	# 			inhCount=inhCount+1
	# 		if edge2=='i':
	# 			inhCount=inhCount+1
	# 		if inhCount==1:
	# 			graph.add_edge(before,after,signal='i')
	# 		else:
	# 			graph.add_edge(before,after,signal='a')
	# 	graph.remove_node(rm)
	# 	addBackNodes.append([rm,before])

	# # 5. rewire nodes that have only one downstream node
	# removeNodeList= [x for x in graph.nodes() if (len(graph.successors(x))==1) ]
	# for rm in removeNodeList:
	# 	if len(graph.successors(x))==1:
	# 		finish=graph.successors(rm)[0]
	# 		for start in graph.predecessors(rm):
	# 			edge1=graph.get_edge_data(start,rm)['signal']
	# 			edge2=graph.get_edge_data(rm,finish)['signal']
	# 			inhCount=0
	# 			if edge1=='i':
	# 				inhCount=inhCount+1
	# 			if edge2=='i':
	# 				inhCount=inhCount+1
	# 			if inhCount==1:
	# 				graph.add_edge(start,finish,signal='i')
	# 			else:
	# 				graph.add_edge(start,finish,signal='a')
	# 		graph.remove_node(rm)
	# 		addBackNodes.append([rm,finish])
	# 6. remove self edges
	for edge in graph.edges():
		if edge[0]==edge[1]:
			graph.remove_edge(edge[0],edge[1])
	return graph, addBackNodes

if __name__ == '__main__':
	# read in options
	parser = argparse.ArgumentParser(prog='BONITA') 
	parser.set_defaults(verbose=False, mode='PA',sep=',')
	parser.add_argument("-v", action="store_true", dest="verbose",  help="output ongoing iterations to screen [default off]")
	parser.add_argument("-m", "--mode", metavar="mode", help="What BONITA functions should be run?")	
	parser.add_argument("-sep", "--sep", metavar="seperator", help="How are columns in datafile specified")	
	parser.add_argument("-t", action='store_const',const='\t', dest="sep",help="Tab delimited?")	
	parser.add_argument("data")
	parser.add_argument("pathways") # 'filtered.c2.cp.kegg.v3.0.symbols.gmt'
	results = parser.parse_args()
	dataName=results.data
	pathwaysName=results.pathways
	verbose=results.verbose
	mode=results.mode

	sss, geneDict, cvDict=readFpkmData(dataName, results.sep) # read in data
	pickle.dump( sss, open( 'sss.pickle', "wb" ) ) # save data in correct format for runs
	findPathways(cvDict, pathwaysName, geneDict) # generate gpickles needed for pathway analysis
