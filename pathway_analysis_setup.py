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
from bioservices import KEGG
import urllib2

# import other pieces of our software
import networkConstructor as nc
from utils import readFpkmData
import string

def xstr(s):
	"""Handle NoneType in strings"""
	if s is None:
		return("Unknown")
	return str(s)
def makeSIF(pathway, keggObject):
	print("Pathway: ")
	print(pathway.encode('utf-8'))
	activationRelations = ["activation", "binding/association", "phosphorylation", "indirect effect", "dissociation"] #Change
	inhibitionRelations = ["inhibition", "dephosphorylation", "dissociation", "ubiquitination"] #Change
	try:
		res = keggObject.parse_kgml_pathway(pathway) #Change
	except:
		print('Failed graph download: '+ str(pathway))
		return []
	sif = []
	#print(res)
	#res2 = {k1.encode("utf-8"):[{k.encode("utf-8"): v.encode("utf-8") for k,v in entry.iteritems()} for entry in v1] for k1,v1 in res.iteritems() }
	#print(res2)
	for rel in res['relations']:
		# types can be PPrel (protein-protein interaction only)
		#print(rel)
		if rel['name'].encode("utf-8") in inhibitionRelations:
			Id1 = rel['entry1']
			Id2 = rel['entry2']
			type1 = res['entries'][[x['id'] for x in res['entries']].index(Id1)]['type']
			type2 = res['entries'][[x['id'] for x in res['entries']].index(Id2)]['type']
			if type1!='gene' or type2!='gene':
				continue
			name1 = xstr(res['entries'][[x['id'] for x in res['entries']].index(Id1)]['gene_names']).split(",")[0]
			name2 = xstr(res['entries'][[x['id'] for x in res['entries']].index(Id2)]['gene_names']).split(",")[0]
			sif.append([name1, -1, name2])
		else:
			Id1 = rel['entry1']
			Id2 = rel['entry2']
			type1 = res['entries'][[x['id'] for x in res['entries']].index(Id1)]['type']
			type2 = res['entries'][[x['id'] for x in res['entries']].index(Id2)]['type']
			#print(type1)
			#print(type2)
			#if type1.encode("utf-8")!='gene' or type2.encode("utf-8")!='gene':
			#	continue
			print(res['entries'])[[x['id'] for x in res['entries']].index(Id1)]
			print(res['entries'][[x['id'] for x in res['entries']].index(Id1)]['gene_names'])
			name1 = xstr(res['entries'][[x['id'] for x in res['entries']].index(Id1)]['gene_names']).split(",")[0]
			name2 = xstr(res['entries'][[x['id'] for x in res['entries']].index(Id2)]['gene_names']).split(",")[0]
			sif.append([str.upper(name1).replace('.',''), 1, str.upper(name2).replace('.',"")])
	sif = np.array([np.array(sif1) for sif1 in sif])
	#print(sif)
	return sif

def sif_to_digraph(pathwayGenes, sif):
	node1, node2, edges = set(sif[0:,0]), set(sif[0:,2]), sif[:, [0,2]]
	nodes=node1.union(node2)
	del node1, node2
	nodes=list(nodes)
	G=nx.DiGraph()
	G.add_nodes_from(nodes)
	for edge in sif:
		node1 = str(edge[0])
		node2 = str(edge[2])
		if (node1 == node2):
			#don't add self edges
			continue
		elif (node1 not in pathwayGenes) or (node2 not in pathwayGenes):
			#don't add nodes that are not in the set of filtered genes
			continue
		elif (node1 == node2):
			#don't add self-loops
			continue
		else:
			G.add_edge(str(edge[0]), str(edge[2]), signal=int(edge[1]))
	#graph post-processing
	#remove singletons/isolates
	G.remove_nodes_from(list(nx.isolates(G)))
	#To do: remove complexes, remove dependences of a node on complexes that include that node (which is a form of self-loop)
	return(G)

def sif_to_digraph2(sif):
	node1, node2, edges = set(sif[0:,0]), set(sif[0:,2]), sif[:, [0,2]]
	nodes=node1.union(node2)
	del node1, node2
	nodes=list(nodes)
	G=nx.DiGraph()
	G.add_nodes_from(nodes)
	for edge in sif:
		node1 = str(edge[0])
		node2 = str(edge[2])
		if (node1 == node2):
			#don't add self edges
			continue
		elif (node1 == node2):
			#don't add self-loops
			continue
		else:
			G.add_edge(str(edge[0]), str(edge[2]), signal=int(edge[1]))
	#graph post-processing
	#remove singletons/isolates
	G.remove_nodes_from(list(nx.isolates(G)))
	#To do: remove complexes, remove dependences of a node on complexes that include that node (which is a form of self-loop)
	return(G)

# read in file with pathway genes and names
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
			nx.write_graphml(graph,coder+'_before.graphml')
			graph=simplifyNetworkpathwayAnalysis(graph, cvDict)
			nx.write_graphml(graph,coder+'.graphml')
			if len(genes.intersection(graph.nodes()))>1:
				nx.write_gpickle(graph,coder+'.gpickle')
				print('nodes: ',str(len(graph.nodes())),',   edges:',str(len(graph.edges())))
				# save the removed nodes and omics data values for just those nodes in the particular pathway
				pathwaySampleList=[{} for q in range(len(geneDict[list(graph.nodes())[0]]))]
				for noder in graph.nodes():
					for jn in range(len(pathwaySampleList)):
						pathwaySampleList[jn][noder]=geneDict[noder][jn]
				pickle.dump( pathwaySampleList, open( coder+"_sss.pickle", "wb" ) )
	else:
		print('not found:')
		print(requester)
		print(lines)

def find_pathways_organism(cvDict, preDefList = [],writeGraphml=True,  organism="hsa"):
	aliasDict, koDict, orgDict = {}, {}, {} # set up empty dictionaries for converting codes
	nc.parseKEGGdict('inputData/ko00001.keg',aliasDict,koDict) # parse the dictionary of ko codes
	try: # try to retrieve and parse the dictionary containing organism gene names to codes conversion
		url=urllib2.urlopen('http://rest.kegg.jp/list/'+organism)
		text=url.readlines()
		# reads KEGG dictionary of identifiers between numbers and actual protein names and saves it to a python dictionary
		for line in text:
			line_split=line.split('\t')
			k=line_split[0].split(':')[1]
			nameline=line_split[1].split(';')
			name=nameline[0]
			if ',' in name:
				nameline=name.split(',')
				name=nameline[0]
				for entry in range(1,len(nameline)):
					aliasDict[nameline[entry].strip()]=name.upper()
			orgDict[k]=name
	except:
		print('Could not get library: ' + organism )
	k = KEGG() # read KEGG from bioservices
	k.organism=organism
	minOverlap=5
	if len(preDefList) ==0:
		pathwayList=list(k.pathwayIds)
	else:
		pathwayList= list(preDefList)

	# set up a converter to retain only numbers from KEGG pathway codes
	allChars=string.maketrans('','')
	noDigits=allChars.translate(allChars, string.digits)

	genes=set(cvDict.keys()) # find the list of genes included in dataset

	for x in pathwayList:
		x=x.replace("path:","")
		code=str(x)
		code= code.translate(allChars, noDigits) # eliminate org letters
		coder=str('ko'+code) # add ko
		graph=nx.DiGraph() # open a graph object
		nc.uploadKEGGcodes([coder], graph, koDict) # get ko pathway
		coder=str(organism+code) # set up with org letters
		nc.uploadKEGGcodes_hsa([coder], graph,orgDict, koDict) # get org pathway
		# check to see if there is a connected component, simplify graph and print if so
		if len(list(nx.connected_component_subgraphs(graph.to_undirected() )))>1: # if there is more than a 1 node connected component, run BONITA
			nx.write_graphml(graph,coder+'_before.graphml')
			if len(genes.intersection(graph.nodes()))> minOverlap: # if there are 5 genes shared
				graph=simplifyNetworkpathwayAnalysis(graph, cvDict) # simplify graph to nodes in dataset
				nx.write_graphml(graph,coder+'.graphml') # write graph out
				nx.write_gpickle(graph,coder+'.gpickle') # write graph out
				print('nodes: ',str(len(graph.nodes())),',   edges:',str(len(graph.edges())))
				# save the removed nodes and omics data values for just those nodes in the particular pathway
				pathwaySampleList=[{} for q in range(len(geneDict[list(graph.nodes())[0]]))]
				for noder in graph.nodes():
					for jn in range(len(pathwaySampleList)):
						pathwaySampleList[jn][noder]=geneDict[noder][jn]
				pickle.dump( pathwaySampleList, open( coder+"_sss.pickle", "wb" ) )

# identify pathways and complete setup for simulation
def findPathwaysHuman(cvDict,gmtName, geneDict):
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
	# # 1. remove self edges
	# # 2. remove complexes and rewire components
	# # 3. remove nodes with no input data
	# # 4. remove dependence of nodes on complexes that include that node

	# 1. remove self edges
	for edge in graph.edges():
		if edge[0]==edge[1]:
			graph.remove_edge(edge[0],edge[1])

	# 2.  remove complexes and rewire components
	removeNodeList= [x for x in graph.nodes() if  '-' in x]
	for rm in removeNodeList:
		for start in graph.predecessors(rm):
			edge1=graph.get_edge_data(start,rm)['signal']
			if edge1=='i':
				for element in rm.split('-'):
					graph.add_edge(start,element,signal='i')
			else:
				for element in rm.split('-'):
					graph.add_edge(start,element,signal='a')
		for finish in graph.successors(rm):
			edge2=graph.get_edge_data(rm,finish)['signal']
			if edge2=='i':
				for element in rm.split('-'):
					graph.add_edge(element,finish,signal='i')
			else:
				for element in rm.split('-'):
					graph.add_edge(element,finish,signal='a')
		graph.remove_node(rm)

	# 3. remove nodes with no input data
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

	# 4. remove dependence of nodes on complexes that include that node
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

	for edge in graph.edges():
		if edge[0]==edge[1]:
			graph.remove_edge(edge[0],edge[1])
	return graph

#Upload KEGG codes modified for human pathways
def uploadKEGGcodes_org(codelist, graph, orgDict, KEGGdict):
#queries the KEGG for the pathways with the given codes then uploads to graph. Need to provide the KEGGdict so that we can name the nodes with gene names rather than KO numbers
	for code in codelist:
		try:
			url=urllib2.urlopen('http://rest.kegg.jp/get/'+code+'/kgml')
		except:
			print('could not read code: ' + code )
			continue
		text=url.readlines()
		readKEGGorg(text, graph, orgDict, KEGGdict)
		#print(code)
		#print(graph.nodes())

def readKEGGorg(lines, graph, orgDict, KEGGdict):
	#read all lines into a bs4 object using libXML parser
	soup = BeautifulSoup(''.join(lines), 'xml')
	groups = {} # store group IDs and list of sub-ids
	id_to_name = {} # map id numbers to names

	for entry in soup.find_all('entry'):
		print(entry)
		entry_split= entry['name'].split(':')
		if len(entry_split)>2:
			if entry_split[0]=='hsa' or entry_split[0]=='ko':
				if entry_split[0]=='hsa':
					useDict=hsaDict
				else:
					useDict=KEGGdict
				nameList=[]
				entry_name=''
				namer=entry_split.pop(0)
				namer=entry_split.pop(0)
				namer=namer.split()[0]
				entry_name=entry_name+useDict[namer] if namer in useDict.keys() else entry_name+namer
				for i in range(len(entry_split)):
					nameList.append(entry_split[i].split()[0])
				for namer in nameList:
					entry_name=entry_name+'-'+useDict[namer] if namer in useDict.keys() else entry_name+'-'+namer
				entry_type = entry['type']
			else:
				entry_name=entry['name']
				entry_type = entry['type']
		else:
			if entry_split[0]=='hsa':
				entry_name=entry_split[1]
				entry_type = entry['type']
				entry_name = hsaDict[entry_name] if entry_name in hsaDict.keys() else entry_name
			elif entry_split[0]=='ko':
				entry_name=entry_split[1]
				entry_type = entry['type']
				entry_name = KEGGdict[entry_name] if entry_name in KEGGdict.keys() else entry_name
			elif entry_split[0]=='path':
				entry_name=entry['name']
				entry_type='path'
			else:
				entry_name=entry['name']
				entry_type = entry['type']

		entry_id = entry['id']
		id_to_name[entry_id] = entry_name

		if entry_type == 'group':
			group_ids = []
			for component in entry.find_all('component'):
				group_ids.append(component['id'])
			groups[entry_id] = group_ids
		else:
			graph.add_node(entry_name, {'name': entry_name, 'type': entry_type})

	for relation in soup.find_all('relation'):
		(color, signal) = ('black', 'a')

		relation_entry1 = relation['entry1']
		relation_entry2 = relation['entry2']
		relation_type = relation['type']

		subtypes = []

		for subtype in relation.find_all('subtype'):
			subtypes.append(subtype['name'])

		if ('activation' in subtypes) or ('expression' in subtypes):
			color='green'
			signal='a'
		elif 'inhibition' in subtypes:
			color='red'
			signal='i'
		elif ('binding/association' in subtypes) or('compound' in subtypes):
			color='purple'
			signal='a'
		elif 'phosphorylation' in subtypes:
			color='orange'
			signal='a'
		elif 'dephosphorylation' in subtypes:
			color='pink'
			signal='i'
		elif 'indirect effect' in subtypes:
			color='cyan'
			signal='a'
		elif 'dissociation' in subtypes:
			color='yellow'
			signal='i'
		elif 'ubiquitination' in subtypes:
			color='cyan'
			signal='i'
		else:
			print('color not detected. Signal assigned to activation arbitrarily')
			print(subtypes)
			signal='a'

		#given: a node ID that may be a group
		#returns: a list that contains all group IDs deconvoluted
		def expand_groups(node_id):
			node_list = []
			if node_id in groups.keys():
				for component_id in groups[node_id]:
					node_list.extend(expand_groups(component_id))
			else:
				node_list.extend([node_id])
			return node_list

		entry1_list = expand_groups(relation_entry1)
		entry2_list = expand_groups(relation_entry2)

		for (entry1, entry2) in it.product(entry1_list, entry2_list):
			node1 = id_to_name[entry1]
			node2 = id_to_name[entry2]
			graph.add_edge(node1,node2, color=color, subtype='/'.join(subtypes), type=relation_type, signal=signal)

	for node in graph.nodes():
		if graph.degree(node)==0:
			graph.remove_node(node)

if __name__ == '__main__':
	# read in options
	parser = argparse.ArgumentParser(prog='BONITA')
	parser.set_defaults(verbose=False, mode='PA',sep=',', org='human', pathways='None', gmt='None')
	parser.add_argument("-v", action="store_true", dest="verbose",  help="output ongoing iterations to screen [default off]")
	parser.add_argument("-m", "--mode", metavar="mode", help="What BONITA functions should be run?")
	parser.add_argument("-sep", "--sep", metavar="seperator", help="How are columns in datafile specified")
	parser.add_argument("-t", action='store_const',const='\t', dest="sep",help="Tab delimited?")
	parser.add_argument("-org", "--org", metavar="org", help="How are columns in datafile specified")
	parser.add_argument("-paths", "--paths", dest="pathways", help="File with list of pathways to be considered each on one line")
	#parser.add_argument("pathways") # 'filtered.c2.cp.kegg.v3.0.symbols.gmt'
	parser.add_argument("-gmt", "--gmt", metavar="gmt", help="GMT file with human pathways from msigDB")
	parser.add_argument("data")

	results = parser.parse_args()
	dataName=results.data
	gmtName=results.gmt
	verbose=results.verbose
	mode=results.mode
	org=results.org
	paths=results.pathways

	sss, geneDict, cvDict=readFpkmData(dataName, results.sep) # read in data
	#pickle.dump( sss, open( 'sss.pickle', "wb" ) ) # save data in correct format for runs
	if org=='human':
		if gmtName=='None':
			print('Please provide either a specific organism for which all of KEGG should be searched using \"-org\" or specify a gmt of specific human pathways using \"paths\"')
		else:
			findPathwaysHuman(cvDict, gmtName, geneDict) # generate gpickles needed for pathway analysis
	else:
		print(org)
		print(paths)
		if paths=='None':
			find_pathways_organism( cvDict, organism=org, writeGraphml=True)
		else:
			inputfile = open(paths, 'r')
			lines = inputfile.readlines()
			pathList= []
			for line in lines:
				for element in line.split(','):
					pathList.append(element.strip())
			find_pathways_organism(cvDict,organism=org,  preDefList=pathList,  writeGraphml=True)