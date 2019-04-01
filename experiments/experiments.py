# import python packages needed in this file
import pickle
from random import random, randint
import networkx as nx
import copy as copy
import operator
import argparse as argparse
from ctypes import *
# import other pieces of our software
import simulation as sim
import GA as ga
from utils import genInitValueList, synthesizeInputs, setupEmptyKOKI

def findNoiseValue(noiseNum):
	if noiseNum==1:
		noise=.01
	elif noiseNum==2:
		noise=.02
	elif noiseNum==3:
		noise=.05
	elif noiseNum==4:
		noise=.1
	elif noiseNum==5:
		noise=.5
	elif noiseNum==6:
		noise=.75
	elif noiseNum==7:
		noise=1.
	elif noiseNum==8:
		noise=2.
	return noise

def findSamples(sampleNum):
	if sampleNum==1:
		samples=2
	elif sampleNum==2:
		samples=3
	elif sampleNum==3:
		samples=4
	elif sampleNum==4:
		samples=5
	elif sampleNum==5:
		samples=10
	elif sampleNum==6:
		samples=15
	return samples
def findRPKNnoise(noiseNum):
	if noiseNum==1:
		noiseEdges=0
	elif noiseNum==2:
		noiseEdges=2
	elif noiseNum==3:
		noiseEdges=5
	elif noiseNum==4:
		noiseEdges=10
	elif noiseNum==5:
		noiseEdges=15
	elif noiseNum==6:
		noiseEdges=20
	elif noiseNum==7:
		noiseEdges=40
	return noiseEdges
# make a list of dictionaries giving values at each node from list of values across samples and a dictionary structure with random numbers
def genSampleList(output, sampleDict, samples, model):
	newSampleList=[]
	for k in range(0,samples):
		newSample=copy.deepcopy(sampleDict[k])
		for j in range(0,len(model.nodeList)):
			newSample[model.nodeList[j]]=output[k][j]
		newSampleList.append(newSample)
	return newSampleList

# run an experiment comparing 
def runExperiment(graph, name, samples, noise, edgeNoise, params):
	#creates a model, runs simulations, then tests reverse engineering capabilities of models in a single function
	#samples is the number of different initial conditions to provide per trial
	#graph specifies the network we are testing. 
	# does everything except params

	# load in C function
	#updateBooler=ctypes.cdll.LoadLibrary('./testRun.so')
	updateBooler=cdll.LoadLibrary('./simulator.so')
	boolC=updateBooler.syncBool 
	params.sample=samples

	sampleList=synthesizeInputs(graph,samples) # get empty list of inputs

	model=sim.modelClass(graph,sampleList, True) # generate empty model
	model.updateCpointers()
	individual=ga.genBits(model) #generate random set of logic rules to start with


	if edgeNoise > 0:
		individual[1]=	[0,1,1,1,1,1,0,0,1,0,1,1,1]


	initModel=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]
	knockoutLists, knockinLists= setupEmptyKOKI(samples)
		
	# generate some simulated samples
	output=ga.runProbabilityBooleanSims(individual[1], model, samples, params.cells, params, knockoutLists, knockinLists, boolC)
	
	# add noise in omics data
	if noise>0:
		multiplier=findNoiseValue(noise)
		for sample in output:
			for i in range(len(sample)):
				sample[i]=min(max(0,sample[i]+multiplier*(random()*2-1)),1)

	# add noise in RPKN
	if edgeNoise > 0:
		newgraph=graph.copy()
		edgelist=newgraph.edges()
		nodelist=newgraph.nodes()
		for newer in range(edgeNoise): # add edgeNoise FP edges
			rand1=randint(0,len(nodelist)-1)
			rand2=randint(0,len(nodelist)-1)
			edgeCandidate=(nodelist[rand1],nodelist[rand2])
			while edgeCandidate in edgelist or edgeCandidate[0]==edgeCandidate[1]:
				rand1=randint(0,len(nodelist)-1)
				rand2=randint(0,len(nodelist)-1)
				edgeCandidate=(nodelist[rand1],nodelist[rand2])
			if random()<.5:
				activity1='a'
			else:
				activity1='i'
			print(edgeCandidate)
			newgraph.add_edge(nodelist[rand1],nodelist[rand2], signal=activity1)
			edgelist.append((nodelist[rand1],nodelist[rand2]))

		print(edgelist)
		print(newgraph.edges())
	else:
		newgraph=graph
	
	# output the initial generated data
	pickle.dump( output, open( name+"_input.pickle", "wb" ) )

	# copy simulated data into right format
	newSampleList=genSampleList(output, sampleList, samples, model)
	testModel=sim.modelClass(newgraph,newSampleList, False)
	testModel.updateCpointers()
	# put initial values into correct format, add to model
	newInitValueList=genInitValueList(newSampleList,testModel)
	testModel.initValueList=newInitValueList
	
	#find rules
	testModel, dev, bruteOut =ga.GAsearchModel(testModel, newSampleList, params, knockoutLists, knockinLists, name, boolC) # run GA
	bruteOut, equivalents, dev2 = ga.localSearch(testModel, bruteOut, newSampleList, params, knockoutLists, knockinLists, boolC) # run local search
	storeModel3=[(testModel.size), list(testModel.nodeList), list(testModel.individualParse), list(testModel.andNodeList) , list(testModel.andNodeInvertList), list(testModel.andLenList),	list(testModel.nodeList), dict(testModel.nodeDict), list(testModel.initValueList)]

	outputList=[individual[1],bruteOut,initModel, storeModel3, equivalents, dev2]
	pickle.dump( outputList, open( name+"_local1.pickle", "wb" ) )


def sampleTester(graph, name, sampleNum):
	samples=findSamples(sampleNum)
	params=sim.paramClass() # load in parameters
	runExperiment(graph, name, samples, 0., 0, params)

def omicsNoiseTester(graph, name, noise):
	params=sim.paramClass() # load in parameters
	runExperiment(graph, name, params.samples, noise,0, params)

def RPKNnoiseTester(graph, name, noiseNum):
	# runs experiment using graph and rule from Liu et al. 2016 along with additional false positive edges
	# loop over number of times we want to generate fake data and perform sequence of events
	params=sim.paramClass() # load in parameters
	noiseEdges=findRPKNnoise(noiseNum)
	runExperiment(graph, name, params.samples, 0. , noiseEdges, params)

def transformTest(graph,name,fileName):
	# can't fit a rule to only one node
	if len(graph.nodes())<2:
		print('not enough overlap')
		return
	
	# load in C function
	updateBooler=cdll.LoadLibrary('./simulator.so')
	boolC=updateBooler.syncBool 

	# load data, params, make empty knockout and knockin lists (no KO or KI in transform tests)
	sampleDict = constructBinInput(fileName)
	params=sim.paramClass()

	# generate turn sample dict into sample list (list of dicts instead of dict of lists)
	keyList=sampleDict.keys()
	sampleList=[{} for i in range(len(sampleDict[keyList[0]]))]
	for i in range(len(sampleList)):
		for key in keyList:
			if key in graph.nodes():
				sampleList[i][key]=sampleDict[key][i]
	
	knockoutLists, knockinLists= setupEmptyKOKI(len(sampleList))

	# generate model encompassing graph and samples to do rule inference on
	model=sim.modelClass(graph,sampleList, False)
	model.updateCpointers()
	# cpy data into correct order for simulation 
	newInitValueList=genInitValueList(sampleList,model)
	model.initValueList=newInitValueList
	print('setup successful')

	# find the rules
	model, dev1, bruteOut =ga.GAsearchModel(model, sampleList, params, knockoutLists, knockinLists, name, boolC)
	bruteOut, equivalents, dev2 = ga.localSearch(model, bruteOut, sampleList, params, knockoutLists, knockinLists, boolC)
	pickle.dump( [[dev1],[dev2],[bruteOut],[model]], open( name+"_output.pickle", "wb" ) )

def findPathways(geneDict):
	returnerlist=[]
	aliasDict={}
	dict1={}
	nc.parseKEGGdicthsa('inputData/hsa00001.keg',aliasDict,dict1)
	dict2={}
	nc.parseKEGGdict('inputData/ko00001.keg',aliasDict,dict2)
	namelist=pa.find_overlaps('filtered.c2.cp.kegg.v3.0.symbols.gmt',geneDict)
	print('num of overlap nodes')
	print(len(namelist))
	for name in namelist:
		returnerlist.append(pa.retrieveGraph(name,aliasDict,dict1,dict2, geneDict))
	return returnerlist

def constructBinInput(filename):
	binput = open(filename, 'r')
        ssDict = {}
	lines = binput.readlines()
	for line in lines:
		geneCount = line.split('\t')
		i = len(geneCount) 
		ssDict[geneCount[0]] = [float(q) for q in geneCount[1:i]]
	return(ssDict)

if __name__ == '__main__':

	# read in arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("graph")
	parser.add_argument("noiseNum")
	parser.add_argument("iterNum")
	results = parser.parse_args()
	graphName=results.graph
	noiseNum=int(results.noiseNum)
	iterNum=int(results.iterNum)
	# save name and print
	outname=graphName[:-8]+results.noiseNum+'_'+results.iterNum
	print(outname) #

	# read graph and run experiment
	graph = nx.read_gpickle(graphName)
	print(len(graph.nodes()))
	omicsNoiseTester(graph,outname,noiseNum)