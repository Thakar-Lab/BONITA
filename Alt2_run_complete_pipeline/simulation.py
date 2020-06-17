#import from other parts of ruleMaker
import utils as utils
#import python modules
from random import random, shuffle, randint
import operator
import networkx as nx
import itertools as itertool
import scipy.stats as regress
import numpy as np
import ctypes 
import math as math

class modelClass:
	def __init__(self,graph, sss, groundTruth): 
		#remove self loops from the graph
		for node in graph.nodes():
			repeat=True
			while(repeat):
				repeat=False
				if node in graph.successors(node):
					graph.remove_edge(node,node)
					repeat=True
		#set up empty lists and dicts for later
		permList=[]
		individualParse=[] # list of the number of shadow and nodes that contribute to each node, in order by index num
		andNodeList=[] #a list of the shadow nodes that represent and relations between incoming edge
		andNodeInvertList=[] # keeps track of which incoming nodes for each node need to be inverted
		andLenList=[] # keeps track of how many nodes are coming into each shadow AND node
		nodeList=graph.nodes()#define the node list simply as the nodes in the graph. 
		nodeDict={} #identifies names of nodes with their index in the node list- provide name, get index
		possibilityLister=[]
		possibilityInverter=[]
		predProbs=[]
		predCorrs=[]
		succnum=[]
		for i in range(0,len(nodeList)):
			nodeDict[nodeList[i]]=i #constructs the node dict so we can easily look up nodes
		counter=int(0) #keeps track of where we are in the generic individual
		initValueList=[] #starting states for nodes
		for j in range(0,len(sss)): #construct empty lists to stick values in later for intiial value list
			initValueList.append([])
		
		#find all possible combinations of upstream contributors for each node. These become the shadow And nodes
		for i in range(0,len(nodeList)):
			permtemp=[]
			predtemp=graph.predecessors(nodeList[i]) # get predecessors of node. 
			succnum.append(len(graph.successors(nodeList[i])))
			possibilitytemp=[nodeDict[predder] for predder in predtemp]
			possibilityLister.append(list(possibilitytemp))
			activity=[]
			predProbsTemp=[]
			for node in possibilitytemp:
				if graph.edge[nodeList[node]][nodeList[i]]['signal']=='a':
					activity.append(False)
				else:
					activity.append(True)
			possibilityInverter.append(list(activity))
			if len(predtemp)>3:
				if groundTruth:
					while len(predtemp)>3: #handle case where there are too many predecessors by truncation
						predtemp.pop(randint(0, len(predtemp)-1))
					preds=predtemp
				else:
					preds=[]
					# select the best predecessors by linear regression		
					slopetemp=[]
					predCorr=[]
					jarray=[ss[nodeList[i]] for ss in sss]
					for k in range(len(predtemp)):
						karray=[ss[predtemp[k]] for ss in sss]
						mi = regress.spearmanr(jarray,karray)
						if math.isnan(mi[0]):
							val1=0.
						else:
							val1=abs(mi[0])
						slopetemp.append(val1)
						predCorr.append(val1)
					# slopetemp2=[abs(temper) for temper in slopetemp]
					while len(preds)<3 and len(slopetemp)>0:
						max1=slopetemp.index(max(slopetemp))
						r=slopetemp.pop(max1)
						adder=predtemp.pop(max1)
						preds.append(adder)
						if len(permtemp)<2:
							permtemp.append(nodeDict[adder])
					predProbsTemp=predCorr
					print('preds')
					print(nodeList[i])
					print(preds)
			else:
				preds=predtemp
			for j in range(0,len(preds)):
				preds[j]=nodeDict[preds[j]]
			# the followign section constructs a list of possible node orders
			# this is accomblished by finding all possible subsets of the list of predecessor nodes
			withNones = zip(preds, itertool.repeat('empty'))
			possibilities=list(itertool.product(*withNones))
			for j in range(0,len(possibilities)):
				possibilities[j]=list(possibilities[j])
				while 'empty' in possibilities[j]:
					possibilities[j].remove('empty')
				while [] in possibilities[j]:
					possibilities[j].remove([])
			while [] in possibilities:
				possibilities.remove([])
			predProbs.append(predProbsTemp)
			# create a list of the activities of each node and store alongside the contributors to each and node for easy reference later
			activities=[] #list to store activities of nodes (a vs i)
			for sequence in possibilities:
				activity=[]
				for node in sequence:
					if graph.edge[nodeList[node]][nodeList[i]]['signal']=='a':
						activity.append(False)
					else:
						activity.append(True)
				activities.append(activity)
			andNodeList.append(possibilities)
			andNodeInvertList.append(activities)
			andLenList.append(len(possibilities))
			permList.append(list(permtemp))
			# construct the list of lengths of possibilties for each node, add to the counter that keeps track of how many bits are necessary
			individualParse.append(counter)
			counter=counter+len(possibilities)
		self.size=counter
		individualParse.append(counter)
		self.rvalues=predProbs # stores the correlations
		self.individualParse=individualParse #index of start value of current node on the individual
		self.andNodeList=andNodeList # shadow and node inputs
		self.andNodeInvertList=andNodeInvertList # keeps track of which incoming nodes for each node need to be inverted
		self.andLenList=andLenList # keeps track of length of above inputOrderList for each node
		self.nodeList=nodeList #define the node list simply as the nodes in the graph. 
		self.nodeDict=nodeDict #identifies names of nodes with their index in the node list.. provide name, get index
		self.initValueList=initValueList #puts an empty and correctly structured initValueList together for later population. 
		self.possibilityList=possibilityLister
		self.possibilityInverter=possibilityInverter
		self.permList=permList
		self.successorNums=succnum
		self.nodeNum=len(nodeList)

		
	def update_upstream(self, node, newUpstreams):
		withNones = zip(newUpstreams, itertool.repeat('empty'))
		possibilities=list(itertool.product(*withNones))
		for j in range(0,len(possibilities)):
			possibilities[j]=list(possibilities[j])
			while 'empty' in possibilities[j]:
				possibilities[j].remove('empty')
			while [] in possibilities[j]:
				possibilities[j].remove([])
		while [] in possibilities:
			possibilities.remove([])
		# create a list of the activities of each node and store alongside the contributors to each and node for easy reference later
		activities=[] #list to store activities of nodes (a vs i)
		for sequence in possibilities:
			activity=[]
			for node1 in sequence:
				if self.possibilityInverter[self.possibilityList[node].index(node1)]=='a':
					activity.append(False)
				else:
					activity.append(True)
			activities.append(activity)
		self.andNodeList[node]=possibilities
		self.andNodeInvertList[node]=activities
	# setup C pointers with correct lengths to pass to simulation software in C
	def updateCpointers(self):
		tempandnoder=[]
		tempandinverter=[]
		for currentNode in range(1000):
			tempAndNodes=[]
			tempandNodeInvertList=[]
			if currentNode<len(self.nodeList):
				tempAndNodes=[xi+[-1]*(3-len(xi)) for xi in self.andNodeList[currentNode]]
				tempandNodeInvertList=[xi+[-1]*(3-len(xi)) for xi in self.andNodeInvertList[currentNode]]
			while(len(tempAndNodes)<7):
				tempAndNodes.append([0,0,0])
				tempandNodeInvertList.append([0,0,0])
			tempandnoder.append(tempAndNodes)
			tempandinverter.append(tempandNodeInvertList)
		self.andNodeInvert=np.array(tempandinverter, dtype=np.intc, order='C')
		self.andNodes=np.array(tempandnoder, dtype=np.intc, order='C')

class paramClass:
	def __init__(self,):    
		self.mutModel=.25
		self.cells=1000
		self.samples=5
		self.generations=120 # generations to run
		self.popSize=24 #size of population
		self.mu= 24 #individuals selected
		self.lambd= 24 #children produced
		self.iters=100 #number of simulations to try in asynchronous mode
		self.genSteps=100 # steps to find steady state with fake data
		self.simSteps=100 # number of steps each individual is run when evaluating
		self.crossoverProb=.6 # prob of crossing over a particular parent
		self.mutationProb=.4 # prob of mutating a particular parent
		self.rewire=False
		self.async=False # run in asynchronous mode
		self.verbose=True
		self.bitFlipProb=.1 # prob of flipping bits inside mutation
		self.sigmaNetwork=0
		self.sigmaNode=0
		self.trials=1
		
def updateBool(currentNode,oldValue,nodeIndividual, model):
	# we update node by updating shadow and nodes then combining them to update or nodes. 
	andNodes=model.andNodeList[currentNode] # find the list of shadow and nodes we must compute before computing value of current nodes
	andNodeInvertList=model.andNodeInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	#update nodes with more than one input
	# first deal with case of simple logic without need of linear regression
	counter =0
	orset=[]
	# go through list of possible shadow and nodes to see which ones actually contribute
	for andindex in range(len(nodeIndividual)):
		if nodeIndividual[andindex]==1:
			# if a shadow and contributes, compute its value using its upstream nodes
			# calculate value of first then use and to append rest in list of predecessors
			newval=oldValue[andNodes[andindex][0]]!=andNodeInvertList[andindex][0]
			for addnode in range(1,len(andNodes[andindex])):
				newval=newval and (oldValue[andNodes[andindex][addnode]]!=andNodeInvertList[andindex][addnode])
			orset.append(newval)
	#combine the shadow and nodes with or operations
	newval=orset.pop()
	for val in orset:
		newval= newval or val
	return int(newval)

def updateFuzzy(currentNode,oldValue,nodeIndividual, model):
	# we update node by updating shadow and nodes then combining them to update or nodes. 
	andNodes=model.andNodeList[currentNode] # find the list of shadow and nodes we must compute before computing value of current nodes
	andNodeInvertList=model.andNodeInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	#update nodes with more than one input
	# first deal with case of simple logic without need of linear regression
	counter = 0
	orset=[]
	# go through list of possible shadow and nodes to see which ones actually contribute
	for andindex in range(len(nodeIndividual)):
		if nodeIndividual[andindex]==1:
			# if a shadow and contributes, compute its value using its upstream nodes
			# calculate value of first then use and to append rest in list of predecessors
			if andNodeInvertList[andindex][0]:
				newval=1-oldValue[andNodes[andindex][0]]
			else:
				newval=oldValue[andNodes[andindex][0]]
			for addnode in range(1,len(andNodes[andindex])):
				if oldValue[andNodes[andindex][addnode]]:
					newval=min(newval,1-andNodeInvertList[andindex][addnode])
				else:
					newval=min(newval,andNodeInvertList[andindex][addnode])
			orset.append(newval)
	#combine the shadow and nodes with or operations
	newval=orset.pop()
	for val in orset:
		newval= max(newval , val)
	return newval

#run a simulation given a starting state in purely python... legacy and unused now
def runBool(individual, model,  simSteps, initValues, params, knockouts, knockins, async):
	if async:
		return asyncBool(individual, model, simSteps, initValues, params.iters, knockouts, knockins)
	else:
		return syncBool(individual, model, simSteps, initValues, knockouts, knockins)

# synchronous Boolean simulation function
def syncBool(individual, model, simSteps, initValues, knockouts, knockins):
	# do simulation. individual specifies the particular logic rules on the model. params is a generic holder for simulation parameters. 
	# set up data storage for simulation, add step 0
	newValue=list(initValues)
	nodeNum=len(model.nodeList)
	simData=np.zeros((nodeNum, simSteps))
	#iterate over number of steps necessary
	for step in range(0,simSteps):
		oldValue=list(newValue)
		for i in range(0,nodeNum):
			#find start and finish for each node to update from the individualParse list
			if i in knockouts:
				temp=0
			elif i in knockins:
				temp=1
			elif model.andLenList[i]==1:
				if oldValue[i]!=model.andNodeInvertList[i][0][0]:
					temp=1
				else:
					temp=0   
			elif model.andLenList[i]==0:
				temp=oldValue[i]
			else:
				if i==len(model.nodeList)-1:
					end= model.size
				else:
					end=model.individualParse[i+1]	 
				if sum(individual[model.individualParse[i]:end])==0:
					temp=oldValue[i]
				else:
					# temp=1
					temp=updateBool(i,oldValue,individual[model.individualParse[i]:end], model)
			newValue[i]=temp
			simData[i,step]=temp
	avg=[.1*np.count_nonzero(simData[i,simSteps-10:simSteps]) for i in range(nodeNum)]
	return avg

#run asyncrhnonoyus simulation and average it over iters trials
def asyncBool(individual, model, simSteps, initValues, iters, knockouts, knockins):
	sum1=[0 for x in range(0,len(initValues))]
	# run iterations with different orderings
	for i in range(0,params.iters):
		# do simulation. individual specifies the particular logic rules on the model. params is a generic holder for simulation parameters. 
		# set up data storage for simulation, add step 0
		newValue=list(initValues)
		# totalNodes=0... turn this on for information criterion

		# set up the sequence of nodes to be updated
		seq=range(0,len(model.nodeList))
		#iterate over number of steps necessary
		for step in range(0,simSteps):
			oldValue=list(newValue)
			#shuffle- async 
			shuffle(seq)
			for seq[i] in range(0,len(model.nodeList)):
				#find start and finish for each node to update from the individualParse list
				if seq[i] in knockouts:
					temp=0
				elif seq[i] in knockins:
					temp=1
				elif model.andLenList[seq[i]]==1:
					if model.andNodeInvertList[seq[i]][0][0]:
						temp=   1- oldValue[seq[i]]
					else:
						temp=  oldValue[seq[i]]
				elif model.andLenList[seq[i]]==0:
					temp=oldValue[seq[i]]
				else:
					if seq[i]==len(model.nodeList)-1:
						end= model.size
					else:
						end=model.individualParse[i+1]	 
					temp=updateBool(seq[i],newValue,individual[model.individualParse[seq[i]]:end],  model)
				newValue[seq[i]]=temp
		for j in range(0,len(sum1)):
			sum1[j]=sum1[j]+newValue[j]
	avgs=list(sum1)
	for i in range(0,len(sum1)):
		avgs[i]=sum1[i]/float(params.iters)
	return avgs

# init value generator for EBNs
def genEBNInitValues(individual, model,sampleProbs):
	#return [True if (random()<sampleProbs[node]) else False for node in range(0,len(sampleProbs))]
	initValues=np.zeros(1000,dtype=np.intc, order='C')
	for node in range(0,len(sampleProbs)):
		if random()<sampleProbs[node]:
			initValues[node]=1
	return initValues

def NP(individual, model, cells, sampleProbs, params, KOs, KIs, boolC):
	if params.async:
		return NPasync(individual, model, cells, sampleProbs, params, KOs, KIs, boolC )
	else:
		return NPsync(individual, model, cells, sampleProbs, params, KOs, KIs, boolC)

# NP simulation code for synchronous simulation... fast bc in C
def NPsync(individual, model, cells, sampleProbs, params, KOs, KIs, syncBoolC):
	cellArray=[]
	simSteps= 100
	
	# set up knockin and knockout lists
	knockins=np.zeros(len(model.nodeList),dtype=np.intc, order='C')
	knockouts=np.zeros(len(model.nodeList),dtype=np.intc, order='C')
	for knocker in KOs:
		knockouts[knocker]=1
	for knocker in KIs:
		knockins[knocker]=1
	# put objects in correct format for passing to C
	nodeIndividual=np.array(individual, dtype=np.intc, order='C')
	indLen=len(nodeIndividual)
	nodeNum=len(model.nodeList)
	individualParse=np.array(model.individualParse, dtype=np.intc, order='C')
	
	andLenList=np.array(model.andLenList, dtype=np.intc, order='C')

	# convert objects into C pointers
	nodeIndividual1=ctypes.c_void_p(nodeIndividual.ctypes.data)
	indLen1=ctypes.c_void_p(indLen)
	andNodes1=ctypes.c_void_p(model.andNodes.ctypes.data)
	individualParse1=ctypes.c_void_p(individualParse.ctypes.data)
	andLenList1=ctypes.c_void_p(andLenList.ctypes.data)
	andNodeInvertList1=ctypes.c_void_p(model.andNodeInvert.ctypes.data)
	nodeNum1=ctypes.c_void_p(nodeNum)
	simSteps1=ctypes.c_void_p(simSteps)
	knockouts1=ctypes.c_void_p(knockouts.ctypes.data)
	knockins1=ctypes.c_void_p(knockins.ctypes.data)
	for j in range(0,cells):
		# run simulation across cells
		initValues=genEBNInitValues(individual, model,sampleProbs) # get initial values for all nodes
		vals=np.zeros(len(model.nodeList), dtype=np.intc, order='C') # initiate output array
		initValues=np.array(initValues,dtype=np.intc, order='C') # setup input array
		initValues1=ctypes.c_void_p(initValues.ctypes.data) # put input array as C pointer
		valsubmit=ctypes.c_void_p(vals.ctypes.data) # put output array into C pointer
		# run simulation and average over last ten steps
		syncBoolC(valsubmit,nodeIndividual1, indLen1, nodeNum1, andLenList1, individualParse1, andNodes1, andNodeInvertList1, simSteps1, initValues1, knockouts1, knockins1)
		cellArray.append(.1*np.copy(vals))
	return [( (1.*np.sum(col))/cells) for col in zip(*cellArray)]

# NP simulation code for asynchronous simulation... slow because no C code yet
def NPasync(individual, model, cells, sampleProbs, params, KOs, KIs):
	cellArray=[]
	simSteps= 100
	# run simulation across many cells
	for j in range(0,cells):
		initValues=genEBNInitValues(individual, model,sampleProbs) # get initial values for all nodes
		asyncResult = asyncBool(individual, model, simSteps, initValues, params.iters, KOs, KIs) # run asynchronous simulation
		cellArray.append(asyncResult) # append results for this cell
	return [(1.*np.sum(col)) / cells for col in zip(*cellArray)]


# wrapper to run update for a single node in C... mainly used for testing, no longer used at all in code
def updateBool2(currentNode,oldValue,nodeIndividual,model):
	# load in library
	updateBooler=ctypes.cdll.LoadLibrary('./simulator.so')
	updateBool3=updateBooler.updateBool  
	# save old value and output arrays
	oldValue=np.array(oldValue,dtype=np.intc)
	nodeIndividual=np.array(nodeIndividual, dtype=np.intc)
	
	# set up and node arrays, convert to C pointers....
	# all of this is obsolete now and is in the model, but this is legacy code so left here for future reference
	indLen=len(nodeIndividual)
	tempAndNodes=[xi+[-1]*(3-len(xi)) for xi in model.andNodeList[currentNode]]
	while(len(tempAndNodes)<7):
		tempAndNodes.append([0,0,0])
	andNodes=np.array(tempAndNodes)
	
	tempandNodeInvertList=[xi+[-1]*(3-len(xi)) for xi in model.andNodeInvertList[currentNode]]
	while(len(tempandNodeInvertList)<7):
		tempandNodeInvertList.append([0,0,0])
	andNodeInvertList=np.array(tempandNodeInvertList)
	
	andLenList=np.array([len(lister) for lister in model.andNodeList[currentNode]], dtype=np.intc)

	currentNode=ctypes.c_void_p(currentNode)
	nodeIndividual=ctypes.c_void_p(nodeIndividual.ctypes.data)
	indLen=ctypes.c_void_p(indLen)
	andNodes=ctypes.c_void_p(andNodes.ctypes.data)
	andLenList=ctypes.c_void_p(andLenList.ctypes.data)
	oldValue=ctypes.c_void_p(oldValue.ctypes.data)
	andNodeInvertList=ctypes.c_void_p(andNodeInvertList.ctypes.data)

	value=updateBool3(currentNode,oldValue,nodeIndividual, indLen,andNodes, andLenList,andNodeInvertList)
	return value