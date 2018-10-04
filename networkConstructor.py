#import external code
import operator
import networkx as nx
import re
import urllib2 
import csv 
import itertools as it
import sys
from bs4 import BeautifulSoup
from random import randint 
# import our code
import simulation as sim
import utils as utils
#definitions from BioPAX level 3 reference manual (http://www.biopax.org/mediawiki/index.php?title=Specification)
#these biopax classes are iteracted over in the biopax methods
edge_classes = ['Interaction', 'GeneticInteraction', 'MolecularInteraction', 'TemplateReaction', 'Control', 'Catalysis', 'TemplateReactionRegulation', 'Modulation', 'Conversion', 'ComplexAssembly', 'BiochemicalReaction', 'Degradation', 'Transport', 'TransportWithBiochemicalReaction']
node_classes = ['PhysicalEntity', 'Complex', 'Dna', 'DnaRegion', 'Protein', 'Rna', 'RnaRegion', 'SmallMolecule', 'Gene']
graph_classes = ['Pathway']

#Upload KEGG codes modified for human pathways
def uploadKEGGcodes_hsa(codelist, graph, hsaDict, KEGGdict):
#queries the KEGG for the pathways with the given codes then uploads to graph. Need to provide the KEGGdict so that we can name the nodes with gene names rather than KO numbers
	for code in codelist:
		url=urllib2.urlopen('http://rest.kegg.jp/get/'+code+'/kgml')
		text=url.readlines()
		readKEGGhsa(text, graph, hsaDict, KEGGdict)
		#print(code)
		#print(graph.nodes())

def readKEGGhsa(lines, graph, hsaDict, KEGGdict):
	#read all lines into a bs4 object using libXML parser
	soup = BeautifulSoup(''.join(lines), 'xml')
	groups = {} # store group IDs and list of sub-ids
	id_to_name = {} # map id numbers to names

	for entry in soup.find_all('entry'):
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

def parseKEGGdict(filename):
	#makes a dictionary to convert ko numbers from KEGG into real gene names
	#this is all file formatting. it reads a line, parses the string into the gene name and ko # then adds to a dict that identifies the two.
	dict={}
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		if line[0]=='D':
			kline=line.split('      ')
			kstring=kline[1]
			kline=kstring.split('  ')
			k=kline[0]
			nameline=line.replace('D      ', 'D')
			nameline=nameline.split('  ')
			namestring=nameline[1]
			nameline=namestring.split(';')
			name=nameline[0]
			dict[k]=name
	return dict

def parseKEGGdicthsa(filename, aliasDict, dict1):
	#reads KEGG dictionary of identifiers between orthology numbers and actual protein names and saves it to a python dictionary
	#extends above method to also keep track of gene names that are identical so we can recover those from input data as well
	#again, read in lines from the file, save ko number and gene name, identify them in the dictionary. 
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		if line[0]=='D':
			kline=line.split('      ')
			kstring=kline[1]
			kline=kstring.split(' ')
			k=kline[0]
			nameline=kline[1]
			nameline=nameline.split(';')
			name=nameline[0]
			if ',' in name:
				nameline=name.split(',')
				name=nameline[0]
				for entry in range(1,len(nameline)):
					aliasDict[nameline[entry].strip()]=name
			dict1[k]=name
	return dict1

def parseKEGGdict(filename, aliasDict, dict1):
	#reads KEGG dictionary of identifiers between orthology numbers and actual protein names and saves it to a python dictionary
	#extends above method to also keep track of gene names that are identical so we can recover those from input data as well
	#again, read in lines from the file, save ko number and gene name, identify them in the dictionary. 
	inputfile = open(filename, 'r')
	lines = inputfile.readlines()
	for line in lines:
		if line[0]=='D':
			kline=line.split('      ')
			kstring=kline[1]
			kline=kstring.split('  ')
			k=kline[0]
			nameline=line.replace('D      ', 'D')
			nameline=nameline.split('  ')
			namestring=nameline[1]
			nameline=namestring.split(';')
			name=nameline[0]
			if ',' in name:
				nameline=name.split(',')
				name=nameline[0]
				for entry in range(1,len(nameline)):
					aliasDict[nameline[entry].strip()]=name
			dict1[k]=name
	return dict1

def readKEGG(lines, graph, KEGGdict):
	#read all lines into a bs4 object using libXML parser
	soup = BeautifulSoup(''.join(lines), 'xml')
	groups = {} # store group IDs and list of sub-ids
	id_to_name = {} # map id numbers to names

	for entry in soup.find_all('entry'):
		entry_name = entry['name'].split(':')[1] if ':' in entry['name'] else entry['name']
		entry_name = KEGGdict[entry_name] if entry_name in KEGGdict.keys() else entry_name
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

def uploadKEGGfiles(filelist, graph, foldername, KEGGdict):
	#upload KEGG files from a particular folder.
	#just provide the folder, file names as a list, and the graph you want to put things in.
	#iteratively calls teh readKEGG method
	for file in filelist:
		inputfile = open(foldername+'/'+file, 'r')
		lines = inputfile.readlines()
		readKEGG(lines, graph, KEGGdict)

def uploadKEGGfolder(foldername, graph, KEGGdict):
#uploads an entire folder of KEGG files given by foldername to the netx graph
#just picks up names from the folder and passes them to the uploadKEGGfiles method
#be sure not to pass in a folder that has files other than the KEGG files in it
	filelist= [ f for f in listdir(foldername+'/') if isfile(join(foldername,join('/',f))) ]
	uploadKEGGfiles(filelist, graph, foldername, KEGGdict)
	
def uploadKEGGcodes(codelist, graph, KEGGdict):
#queries the KEGG for the pathways with the given codes then uploads to graph. Need to provide the KEGGdict so that we can name the nodes with gene names rather than KO numbers
	for code in codelist:
		url=urllib2.urlopen('http://rest.kegg.jp/get/'+code+'/kgml')
		text=url.readlines()
		readKEGG(text, graph, KEGGdict)

# creates graphs as gpickles for experiments using IFNG networks
def ifngStimTestSetup():
	aliasDict={}
	dict1={}
	parseKEGGdicthsa('inputData/hsa00001.keg',aliasDict,dict1)
	dict2={}
	parseKEGGdict('inputData/ko00001.keg',aliasDict,dict2)
	# read in list of codes then load them into network
	inputfile = open('inputData/ID_filtered_IFNGpathways2.txt', 'r')
	lines = inputfile.readlines()

	# iterate through and create gpickle of graph for all IFNG networks
	for code in lines:
		graph=nx.DiGraph()
		coder=str('ko'+code[:-1])
		uploadKEGGcodes([coder], graph, dict2)
		coder=str('hsa'+code[:-1])
		uploadKEGGcodes_hsa([coder], graph,dict1, dict2)
		# if the code has an interesting logic rule to find, run the algorithm... 
		print(coder)
		checker=False
		for x in graph.in_degree():
			if x>1:
				checker=True
		if(checker):
			nx.write_graphml(graph,coder+'.graphml')
			nx.write_gpickle(graph,coder+'.gpickle')
			print(coder)