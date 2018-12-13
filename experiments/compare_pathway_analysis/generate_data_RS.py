# import python packages needed in this file
import pickle
from random import random, randint, shuffle
import networkx as nx
import copy as copy
import operator
import argparse as argparse
from ctypes import *
import csv as csv
import numpy as numpy
from scipy.stats import variation
from math import log, exp
from sets import Set

# import other pieces of our software
import simulation as sim
import GA as ga
from pathway_analysis_score_nodes import calcImportance
from utils import genInitValueList, synthesizeInputs, setupEmptyKOKI
from random import random, seed, shuffle, gauss
from pathway_analysis_score_pathways import scorePathway

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
		geneDict[data[i][0]]=tempDatalist
		for j in range(0,len(data[i])-1):
			sampleList[j][str.upper(data[i][0])]=float(data[i][j+1])/maxdata
	return sampleList, geneDict, cvDict


def outputData(control, experiment, genelist,filename, geneDict):
	count=0
	pathset=Set(experiment[0].keys())
	datatable=[]
	datatable.append(['','c1','c2','c3','c4','c5','e1','e2','e3','e4','e5'])
	geneset=Set(genelist)
	for gene in geneset.intersection(pathset):
		temp=[gene]
		maxValue = numpy.max(geneDict[gene])
		temp.extend([ maxValue*sampleLister[gene] for sampleLister in control])
		temp.extend([ maxValue*sampleLister[gene] for sampleLister in experiment])
		datatable.append(temp)
	for gene in geneset.difference(pathset):
		temp=[gene]
		temp.extend([geneDict[gene][j] for j in range(10)])
		datatable.append(temp)
	with open(filename, "wb") as f:
		writer = csv.writer(f)
		writer.writerows(datatable)

# put together the output of a simulation into a sampleList
def compileOuts(output,sampleList, sampleNum, model):
	newSamplelists=[]
	for k in range(0,sampleNum):
		newSamplelist=copy.deepcopy(sampleList[k])
		for j in range(0,len(model.nodeList)):
			newSamplelist[model.nodeList[j]]=float(output[k][j])
		newSamplelists.append(newSamplelist)
	return newSamplelists

def PAtester(graph, name):
	controls=5
	experimentals=5
	true=10
	false=true
	params=sim.paramClass()
	sampleLists,geneDicts,cvDicts= [],[],[]
	# import starting points
	for i in range(1,11):
		sampleList, geneDict, cvDict=readFpkmData('neg_binom_gen_'+str(i)+'.csv', ',') # read in data
		sampleLists.append(sampleList)
		geneDicts.append(geneDict)
		cvDicts.append(cvDict)
	
	for j in range(10): # iterate over imported starting points
		# generate model
		# loop over number of times we want to generate fake data and perform sequence of events
		# generate Boolean model for this trial
		controlSampleList=sampleLists[j][0:5]
		genelist=geneDicts[j].keys()
		for perturbation in [0,5,10,15,20]: 
			tSampleList=list(sampleLists[j][5:10])
			perturbationSize=2.**(-.1*perturbation)
			for i in range(5):
				# generate values across samples
				for node in graph.nodes():
					if len(graph.predecessors(node))==0:
						tSampleList[i][node]=min(max(0,sampleLists[j][i+5][node]*(perturbationSize)),1)
			outputData(controlSampleList, tSampleList, genelist,name+str(perturbation)+'_true_'+str(j)+'.csv', geneDicts[j])

if __name__ == '__main__':	
	import time	
	start_time = time.time()
	parser = argparse.ArgumentParser()
	parser.add_argument("graph")
	results = parser.parse_args()
	graphName=results.graph
	name=graphName[:-8]+'_'
	graph = nx.read_gpickle(graphName)
	PAtester(graph, name)
	print("--- %s seconds ---" % (time.time() - start_time))