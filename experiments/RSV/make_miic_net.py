#import python packages
import os as os
import pickle
import argparse as argparse
import numpy as np
import csv
import math
from sets import Set
from random import randint
import scipy
import pandas as pd
import requests
import networkx as nx
import gc as gc
from ctypes import *

# import other pieces of our software
from GA import GAsearchModel, localSearch
from pathway_analysis_setup import readFpkmData
from pathway_analysis_score_nodes import calcImportance
from utils import genInitValueList, setupEmptyKOKI, writeModel, Get_expanded_network
from simulation import paramClass, modelClass

# write out graphs with relative abundance, importance scores on them
def outputGraphMiic(graph, rules, pathImportances):
	# write original graph with annotations
	original=graph.copy()
	nx.set_node_attributes(original,'Display Name',{k: k for k in original.nodes()})
	nx.set_node_attributes(original,'andNode',{k: 0 for k in original.nodes()})
	nx.set_node_attributes(original,'IS',{k: float(pathImportances[k]) for k in original.nodes()})
	
	# write graph of rules with annotations
	ruleGraph=Get_expanded_network(rules.split('\n'),equal_sign='*=')
	nx.set_node_attributes(ruleGraph,'IS',{k: float(pathImportances[k]) if k in original.nodes() else 0. for k in ruleGraph.nodes()})
	
	# output graphs
	nx.write_graphml(original,'importance_score_graph.graphml')
	nx.write_graphml(ruleGraph,'importance_score_rules.graphml')

if __name__ == '__main__':
	import time
	start_time = time.time()
	
	# get arguments
	parser = argparse.ArgumentParser(prog='BONITA') 
	parser.set_defaults(verbose=False, mode='PA',sep=',')
	parser.add_argument("-v", action="store_true", dest="verbose",  help="output ongoing iterations to screen [default off]")
	parser.add_argument("-sep", "--sep", metavar="seperator", help="How are columns in datafile specified")	
	parser.add_argument("-t", action='store_const',const='\t', dest="sep",help="Tab delimited?")	
	parser.add_argument("data")
	parser.add_argument("miic_graphml")
	results = parser.parse_args()
	dataName=results.data
	graphName=results.miic_graphml
	verbose=results.verbose
	
	# read in graph and fix node labels, edge labels
	graph=nx.read_graphml(graphName)
	names=nx.get_node_attributes(graph, 'name')
	dicty1={key.encode('utf-8'): value.encode('utf-8') for key, value in names.iteritems()}
	graph=nx.relabel_nodes(graph, dicty1, copy=True)
	signs=nx.get_edge_attributes(graph, 'sign')
	activities={}
	for e in graph.edges_iter():
		if signs[e]=='+':
			activities[e]= 'a' 
		elif signs[e]=='-':
			activities[e]='i'
		else:
			print('not found: '+signs[e])
	nx.set_edge_attributes(graph,'signal',activities)


	# read in C library
	updateBooler=cdll.LoadLibrary('./testRun.so')
	boolC=updateBooler.syncBool 
	sampleList, geneDict, cvDict=readFpkmData(dataName, results.sep)
	genes=Set(geneDict.keys())
	newOverlap=genes.intersection(set(graph.nodes()))

	# graph, addLaterNodes=PA.simplifyNetworkpathwayAnalysis(graph, geneDict)
	print(len(graph.nodes()))
	sampleList=[{} for q in range(len(geneDict[list(newOverlap)[0]]))]
	for noder in newOverlap:
		for jn in range(len(sampleList)):
			sampleList[jn][noder]=geneDict[noder][jn]
	print(sampleList)
	# zero out old arrays so memory can be recycled
	genes=[]
	geneDict=[]
	cvDict=[]
	gc.collect()

	# load params, model, KO and KI lists, initValueLists
	params=paramClass()
	model=modelClass(graph,sampleList, False)
	model.updateCpointers()
	samples=params.samples
	knockoutLists, knockinLists= setupEmptyKOKI(len(sampleList))
	newInitValueList=genInitValueList(sampleList,model)
	model.initValueList=newInitValueList

	# run GA, local search
	model1, dev, bruteOut =GAsearchModel(model, sampleList, params, knockoutLists, knockinLists, 'RSV_miic_out', boolC) # run GA
	bruteOut1, equivalents, dev2 = localSearch(model1, bruteOut, sampleList, params, knockoutLists, knockinLists, boolC) # run local search
	storeModel3=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]
	outputList=[bruteOut1,dev, storeModel3, equivalents, dev2]
	pickle.dump( outputList, open( 'RSV_miic_out'+"_local1.pickle", "wb" ) ) # output rules
	# calculate importance scores and output
	scores1=calcImportance(bruteOut1,params,model, sampleList,knockoutLists, knockinLists, boolC)
	pickle.dump( scores1, open( 'RSV_miic_out'+"_scores1.pickle", "wb" ) )
	ImportanceVals={}
	for node in range(len(storeModel3[1])): 
		ImportanceVals[storeModel3[1][node]]=float(scores1[node])

	# write out rules
	modeler=writeModel(bruteOut1,model1)
	rules=modeler.replace('not','~').split('\n')
	outputGraphMiic(graph, modeler, ImportanceVals)

	print("--- %s seconds ---" % (time.time() - start_time))