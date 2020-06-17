# import necessary python packages
import numpy as numpy
from random import random
import csv as csv
import networkx as nx
from scipy.stats import variation
import simulation as sim
from ctypes import *


# read rpm or fpkm data into format necessary for BONITA simulations and pathway analysis
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
		for j in range(1,len(data[i])):
			tempDatalist.append(float(data[i][j]))
		maxdata=numpy.max(tempDatalist)
		cvDict[data[i][0]]=variation(tempDatalist)
		if maxdata==0:
			maxdata=1.
		geneDict[data[i][0]]=[temperDataPoint/maxdata for temperDataPoint in tempDatalist]
		for j in range(0,len(data[i])-1):
			sampleList[j][str.upper(data[i][0])]=float(data[i][j+1])/maxdata
	return sampleList, geneDict, cvDict

# writes rules as a network
def Get_expanded_network(rules,equal_sign='*='):
	'''
	The code is written by Gang Yang, Department of Physics, Penn State University if not specified.

  Return the expanded network for a given Boolean network model.
  The Boolean network model is a DiGraph object in the output format of form_network().
  The Boolean network model can be generated through form_network function by reading a text file in the Booleannet format.
  The Boolean rules will first be converted to a disjuctive normal form before generating the expanded network.

  Parameters
  ----------
  Gread     : the given Boolean network model
  prefix='n': prefix to encode the node name to avoid one node's name is a part of another node's name
  suffix='n': suffix to encode the node name to avoid one node's name is a part of another node's name
              e.g. node name '1' will become 'n1n' in the returned result
  equal_sign: the equal sign of the rule in the returned result, whose default value follows the Booleannet format

  Returns
  -------
  The expanded network for the given Boolean network model.
	'''
	composite_nodes=[]
	G_expand=nx.DiGraph()
	for line in rules:
		child, update_rule=line.split(equal_sign) #correctly annootate child, rule
		update_rule=update_rule # remove white space from parents
		if update_rule[0]=='(' and update_rule[-1]==')': # remove parens from parent
			update_rule=update_rule[1:-1]
		if child[0]=='~':# figure out parity of node
			normal_child=child[1:].strip()
		else:
			normal_child=child[:].strip()
		if 'or' in update_rule:
			parents=update_rule.split(' or ')
		else:
			parents=[update_rule]
		parents.sort()
		for parent in parents:
			parent=parent.replace('not ','~').replace('(','').replace(')','').strip()
			if 'and' in parent:
				composite_node=parent.replace(' and ','_').strip()
				composite_nodes.append(composite_node)
				G_expand.add_edge(composite_node,child)
				for component in composite_node.split('_'):
					G_expand.add_edge(component.strip(),composite_node)
			elif not parent==child:
				G_expand.add_edge(parent,child)
	for node in G_expand.nodes():
		if node[0]=='~' and not '_' in node:
			G_expand.add_edge( node[1:], node)
	nx.set_node_attributes(G_expand,'Display Name',{k:' ' if k in composite_nodes else k for k in G_expand.nodes() })
	nx.set_node_attributes(G_expand,'andNode',{k: 1 if k in composite_nodes else 0 for k in G_expand.nodes()})

	edgedict={}
	for edge in G_expand.edges():
		edgedict[edge]='a'
	nx.set_edge_attributes(G_expand, 'signal', edgedict)

	for node in G_expand.nodes():
		if node[0]=='~' and not '_' in node:
 			for downstream in G_expand.successors(node):
				G_expand.add_edge( node[1:], downstream, signal='i')
			G_expand.remove_node(node)
	return G_expand.copy()

# make empty list representing no knockouts or knockins
def setupEmptyKOKI(samples):
	knockoutLists=[]
	knockinLists=[]
	for q in range(samples):
		knockoutLists.append([])
		knockinLists.append([])
	return knockoutLists, knockinLists

# use dictionaries of values at each node for each sample to construct initial value list for the model
def genInitValueList(newSampleList,model):
	newInitValueList=[]
	for j in range(0,len(newSampleList)):
		newInitValueList.append([])
	for j in range(0,len(model.nodeList)):
		for k in range(0,len(newSampleList)):
			ss=newSampleList[k]
			if  model.nodeList[j] in newSampleList[0]:
				newInitValueList[k].append(ss[model.nodeList[j]])
			else:
				newInitValueList[k].append(0.5)
	return newInitValueList

def findEnd(node,model):
	if node==len(model.nodeList)-1:
		end= model.size
	else:
		end=model.individualParse[node+1]
	return end

def genRandBits(individualLength): #makes a random bitstring
	arr = numpy.random.randint(2, size=(int(individualLength),))
	return list(arr) 

def bitList(n, x):
	templist=[1 if digit=='1' else 0 for digit in bin(n)[::-1]]
	while len(templist)<x:
		templist.append(0)
	while(len(templist))>x:
		templist.pop()
	return templist

def loadFpkms(filename): #loads data from fpkms tab delimited csv file
	with open(filename) as csvfile:
		data=[]
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			data.append(row)
		return data

def sortFpkms(data): #puts fpkms data into appropriate lists of steady state dictionaries following normalization to largest value (as denominator)
	sss=[]
	for j in range(1,len(data[1])):
		sss.append({})
	for i in range(1,len(data)):
		maxdata=0
		for j in range(1,len(data[i])):
			if float(data[i][j])>maxdata:
				maxdata=float(data[i][j])
		if maxdata==0:
			maxdata=1
		for j in range(0,len(data[i])-1):
			sss[j][str.upper(data[i][0])]=float(data[i][1])/maxdata
	return sss

def synthesizeInputs(graph,samples): # generates synthetic completely random inputs for simulation to steady state
	sss=[]
	for i in range(0,samples):
		sss.append({})
	for node in graph.nodes():
		for i in range(0,samples):
			sss[i][node]=random()
	return sss

def writeModel(individual, model):
	#iterate over nodes to generate a BooleanNet representation for the entire model
	addString=''
	for i in range(0,len(model.nodeList)):
		addString=addString+writeNode(i,individual[model.individualParse[i]:model.individualParse[i+1]], model)
		addString=addString+'\n'
	return addString[:-1]

def writeBruteNode(currentNode,individual,model):
	padindividual=[0 for x in range(0,model.individualParse[currentNode][0])]
	padindividual.extend(individual)
	return (writeNode(currentNode, padindividual,model))

def writeNode(currentNode,nodeIndividual, model):
	#write out evaluation instructions in BooleanNet format. 
	# This follows the exact same code as updateNode (for switch=0), but writes a string instead of actually updating the values of the nodes
	andNodes=model.andNodeList[currentNode] # find the list of shadow and nodes we must compute before computing value of current nodes
	andNodeInvertList=model.andNodeInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	writenode=''+model.nodeList[currentNode]+'*=' # set up the initial string to use to write node

	print(currentNode)
	print(andNodes)
	print(andNodeInvertList)
	print(nodeIndividual)
	if model.andLenList[currentNode]==0 or sum(nodeIndividual)==0:
		return writenode + ' ' + model.nodeList[currentNode] #if no inputs, maintain value
	elif len(andNodes)==1: 
		#if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
		value=''	
		#if only one input, then set to that number
		if andNodeInvertList[0][0]==0:
			value= value + model.nodeList[andNodes[0][0]]
		else:
			value= value+ 'not ' + model.nodeList[andNodes[0][0]]
		return writenode + value 
	else:
		#update nodes with more than one input

		# first deal with case of simple logic without need of linear regression
		orset=[]
		# go through list of possible shadow and nodes to see which ones actually contribute
		for andindex in range(len(nodeIndividual)):
			newval='('
			if nodeIndividual[andindex]==1:
				# if a shadow and contributes, compute its value using its upstream nodes
				if andNodeInvertList[andindex][0]:
					newval=newval+'not '
				newval=newval+model.nodeList[andNodes[andindex][0]]
				for addnode in range(1,len(andNodes[andindex])):
					newval= newval + ' and '
					if andNodeInvertList[andindex][addnode]:
						newval=newval+' not '
					newval=newval+model.nodeList[andNodes[andindex][addnode]]
				orset.append(newval +')')
			#combine the shadow and nodes with or operations
		writenode=writenode + orset.pop()
		for val in orset:
			writenode = writenode + ' or ' + val
		print(writenode)
		return writenode

def LiuNetwork1Builder():
	graph = nx.DiGraph()
	graph.add_edge('g','k', signal='a')	
	graph.add_edge('h','j', signal='a')
	graph.add_edge('j','c', signal='i')	
	graph.add_edge('f','k', signal='i')
	graph.add_edge('a','c', signal='a')
	graph.add_edge('b','d', signal='a')
	graph.add_edge('c','f', signal='a')	
	graph.add_edge('c','h', signal='a')
	graph.add_edge('d','f', signal='a')	
	graph.add_edge('d','g', signal='a')
	return graph

# write model from rule set
def makeModelRules(rules,sss,equal_sign='*='):
	graph=nx.DiGraph()
	andNodeList=[]
	nodeListTemp=[]
	for rule in rules:
		andNodeTemp=[]
		ruler=rule.strip('( )\t\n')
		startNode=ruler.split(equal_sign)[0].strip('( )\t')
		nodeListTemp.append(startNode)
		ruler=ruler.split(equal_sign)[1]
		if 'or' in ruler:
			rulers=ruler.split('or')
		else:
			rulers=[ruler]
		for ruler in rulers:
			andNode=[]
			if 'and' in ruler:
				andRules=ruler.split('and')
			else:
				andRules=[ruler]
			for andRule in andRules:
				temprule=andRule.strip('( )\t')
				if 'not' in andRule:
					graph.add_edge(temprule[3:].strip('( )\t'),startNode,attr_dict={'signal':'i'})
					andNode.append(temprule[3:].strip('( )\t'))
				else:
					andNode.append(temprule)
					graph.add_edge(temprule,startNode,attr_dict={'signal':'a'})
			andNode.sort()
			andNodeTemp.append(andNode)
		andNodeList.append(andNodeTemp)
	model=sim.modelClass(graph, sss, True)
	individual=[]
	for i in range(len(model.nodeList)):
		nodeTemp=nodeListTemp.index(model.nodeList[i])
		for j in range(0,model.individualParse[i+1]-model.individualParse[i]):
			tempAndNode=[model.nodeList[node] for node in model.andNodeList[i][j]]
			tempAndNode.sort()
			if  tempAndNode in andNodeList[nodeTemp]:
				individual.append(1)
			else:
				individual.append(0)
	return model, individual, graph

# test whether the simulation code is working properly
def simTest():
	sampleList, geneDict, cvDict=readFpkmData('testInput.txt', '\t')
	with open('testRules.txt') as csvfile:
		model, individual, graph= makeModelRules(csvfile.readlines(),sss=sampleList,equal_sign='=')
	boolValues1=genInitValueList(sampleList,model)
	boolValues2=[]
	updateBooler=cdll.LoadLibrary('./simulator.so')
	boolC=updateBooler.syncBool 
	params=sim.paramClass()
	model.initValueList=boolValues1
	model.updateCpointers()
	KOs,KIs=setupEmptyKOKI(len(sampleList))
	for j in range(len(boolValues1)):
		boolValues2.append(sim.NPsync(individual, model, params.cells, boolValues1[j], params, KOs[j], KIs[j], boolC, True))
	print(boolValues2)
	sampleList3, geneDict3, cvDict3=readFpkmData('testOutput.txt', '\t')
	boolValues3=genInitValueList(sampleList3,model)
	print(boolValues3)
	print(boolValues3==boolValues2)