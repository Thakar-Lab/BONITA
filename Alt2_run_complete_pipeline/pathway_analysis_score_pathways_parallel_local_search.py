# import necessary modules
import argparse as argparse
import operator
import networkx as nx
import pickle
from ctypes import *
from simulation import paramClass, modelClass, NPsync
from utils import genInitValueList, setupEmptyKOKI, writeModel
from GA import GAsearchModel, localSearch
from parallel_local_search import modelHolder
from pathway_analysis_score_nodes import *
import time
#import python packages
import os as os
import pickle
import argparse as argparse
import numpy as np
from scipy.stats import variation
import csv
import math
from random import randint
import scipy
import pandas as pd
import requests
import networkx as nx
import operator
import math as math
# import other pieces of our software
import networkConstructor as nc
from utils import writeModel, Get_expanded_network

# write out graphs with relative abundance, importance scores on them
def outputGraphs(pathway, RAval, comparator, pathImportances, ERS):
	# write original graph with annotations
	original=pathway[3].copy()
	nx.set_node_attributes(original,'Display Name',{k: k for k in original.nodes()})
	nx.set_node_attributes(original,'andNode',{k: 0 for k in original.nodes()})
	nx.set_node_attributes(original,'RA',{k: RAval[k] for k in original.nodes()})
	nx.set_node_attributes(original,'ERS size',{k: ERS[k] for k in original.nodes()})
	# print(maximportance)
	maximportance=max([(pathImportances[v]) for v in pathImportances])
	if maximportance==0:
		maximportance=1.
	nx.set_node_attributes(original,'IS',{k: (float(pathImportances[k])/maximportance) for k in original.nodes()})
	nodeStorage=[]
	for node in original.nodes():
		nodeStorage.append(original.node[node])
	df=pd.DataFrame(nodeStorage)
	df.to_csv('temp_series1_net_scores.csv')
	# write graph of rules with annotations
	ruleGraph=Get_expanded_network(pathway[2][0].split('\n'),equal_sign='*=')
	nx.set_node_attributes(ruleGraph,'RA',{k: RAval[k] if k in original.nodes() else 0. for k in ruleGraph.nodes()})
	nx.set_node_attributes(ruleGraph,'IS', {k: (float(pathImportances[k])/maximportance) if k in original.nodes() else 0. for k in ruleGraph.nodes()})
	nx.write_graphml(original,comparator+'/'+pathway[0]+'.graphml')
	nx.write_graphml(ruleGraph,comparator+'/'+pathway[0]+'_rules.graphml')

# gets KEGG converter from pathway codes to names
def retrievePathKey():
	pathDict={}
	requester='http://rest.kegg.jp/list/pathway'
	r=requests.get(requester)
	lines=r.text
	for line in lines.split('\n'):
		pieces=line.split('\t')
		if(len(pieces[0])>8):
			pathDict[str(pieces[0][8:])]=str(pieces[1])
		else:
			print(pieces)
	return pathDict

# finds relative abundances across comparators
def makeRA(data,comparison,groups):
	RAdict={}
	group1=comparison[0]
	group2=comparison[1]
	for element in data:
		mean1=np.mean([data[element][temp] for temp in groups[group1]])
		mean2=np.mean([data[element][temp] for temp in groups[group2]])
		if mean1<=0 or mean2<=0:
			if mean1==0 and mean2==0:
				RAdict[element]=0.
				print('zero mean for '+ element+'.  means: '+str(mean1)+' '+str(mean2)+'. RA set to zero')
			else:
				print('zero mean for '+ element+'.  means: '+str(mean1)+' '+str(mean2)+'. Low replaced by .1')
				print(data[element])
				if mean1>mean2:
					RAdict[element]=math.log(mean1,2)-math.log(.1,2)
				else:
					RAdict[element]=math.log(.1,2)-math.log(mean2,2)
		else:
			differ=math.log(mean1,2)-math.log(mean2,2)
			RAdict[element]=differ
	return RAdict

# finds pathways that should be compared
def findPathwayList():
	pathways=[]
	codes=[]
	# looks for pathways that have gpickles generated originally
	for file in os.listdir("gpickles"):
		if file.endswith(".gpickle"):
			if os.path.isfile('pickles/'+file[:-8]+'_1_local1.pickle'):
				codes.append(file[:-8])
			else:
				print(file[:-8]+' has no output')
	print(codes)
	# for each of these pathways, we find the output of the rule determination and scoring procedures and put them together. 
	for code in codes:
		pathVals=[]
		rules=[]
		for i in range(1,6):
			[bruteOut1,dev,storeModel, storeModel3, equivalents, dev2]=pickle.Unpickler(open( 'pickles/'+code+'_'+str(i)+'_local1.pickle', "rb" )).load()
			model=modelHolder(storeModel3)
			pathVals.append(pickle.Unpickler(open( 'pickles/'+code+'_'+str(i)+'_scores1.pickle', "rb" )).load())
			rules.append(writeModel(bruteOut1, model))
		graph = nx.read_gpickle("gpickles/"+code+".gpickle")
		ImportanceVals={} # average importance vals over trials
		for node in range(len(storeModel[1])): 
			ImportanceVals[storeModel[1][node]]=float(np.mean([pathVals[i][node] for i in range(5)]))
		# add nodes removed during network simplification back in
		pathways.append([code,ImportanceVals, rules, graph])
	return pathways

def findPathwayList_justRavenPathway():
	pathways=[]
	codes=[]
	# looks for pathways that have gpickles generated originally
	# filename=+'IS_pickles/'moldelNum'_'+str(iterNum)+"_scores1.pickle"

	# for each of these pathways, we find the output of the rule determination and scoring procedures and put them together. 
	
	for code in ['IS_pickles/IS']:
		pathVals=[]
		rules=[]
		equivLengthAccumulated=[]
		for i in range(1,6): 
			[storeModel1, knockoutLists, knockinLists, individual, equivs]=pickle.Unpickler(open( code+'_'+str(i)+'_setup.pickle', "rb" )).load()
			model=modelHolder(storeModel1)
			values=[]
			for j in range(1,11):
				tempvals=pickle.Unpickler(open( code+'_'+str(i)+'_s_'+str(j)+'_scores1.pickle', "rb" )).load()
				values.append(tempvals)
			newvalues=[np.mean([values[j][k] for j in range(0,10)]) for k in range(0,len(values[0]))]
			variances=[variation([values[j][k] for j in range(0,10)]) for k in range(0,len(values[0]))]
			print(variances)
			pathVals.append(newvalues)
			individualOut=[]
			EquivLengths=[]
			for equiv in equivs:
				EquivLengths.append(len(equiv))
				individualOut.extend(equiv[0])
			rules.append(writeModel(individualOut, model))
			equivLengthAccumulated.append(EquivLengths)
		ERSsize={}
		graph = nx.read_gpickle("temp_series1_net"+".gpickle")
		ImportanceVals={} # average importance vals over trials
		for node in range(len(model.nodeList)): 
			ERSsize[model.nodeList[node]]= np.mean([equivLengthAccumulated[q][node] for q in range(len(equivLengthAccumulated))])
			ImportanceVals[model.nodeList[node]]=float(np.mean([pathVals[i][node] for i in range(5)]))
		# add nodes removed during network simplification back in
		pathways.append(['COVID-CONTROL',ImportanceVals, rules, graph,ERSsize])
	return pathways

# read in Omics data
def readFpkm(dataName,delmited):
	data=[]
	with open(dataName) as csvfile:
		data={}
		reader = csv.reader(csvfile, delimiter=delmited)
		firstline=reader.next()
		for row in reader:
			data[row[0]]=[float(row[k]) for k in range(1,len(row))]
	firstline.pop(0)
	
	# identify positions of each sample in the data
	colNums={}
	for item in range(len(firstline)):
		colNums[firstline[item]]=item
	return data, colNums

# read in contrasts to be used
def readDiffs(diffName,delmited):
	# open the set of differences to be considered
	with open(diffName) as csvfile:
		diffs=[]
		reader = csv.reader(csvfile, delimiter=delmited)
		for row in reader:
			diffs.append(row)
	return diffs

# read in matrix telling characteristics of each sample
def readMatrix(matrixName,delmited, colNums):
	# open and analyze matrix of group memberships
	with open(matrixName) as csvfile:
		matrix=[]
		reader = csv.reader(csvfile, delimiter=delmited)
		groupNames=reader.next()
		groupNames.pop(0)
		groups={}
		for name in groupNames:
			groups[name]=[]
		for row in reader:
			groupname=row.pop(0)
			for i in range(len(row)):
				if int(row[i])==1:
					groups[groupNames[i]].append(colNums[groupname])
	return groups

# do pathway analysis! store in one folder for each comparison
def analyze_pathways(diffName, matrixName, dataName, delmited):
	pathList=findPathwayList() # identify pathways under consideration
	data, colNums= readFpkm(dataName,delmited)	# read in fpkm data
	diffs= readDiffs(diffName,delmited) # read in difference to be considered
	groups=readMatrix(matrixName,delmited, colNums) # read design matrix

	# create an index of relative activities for all comparisons
	csvmaker=[]
	RAvals=[]
	comparisonStrings=[]
	# pathDict=retrievePathKey()
	for comparison in diffs:
		comparisonStrings.append(comparison[0]+'-'+comparison[1])
		RAvals.append(makeRA(data,comparison,groups))
		if not os.path.exists(comparison[0]+'-'+comparison[1]):
			os.makedirs(comparison[0]+'-'+comparison[1])
	CVdict={}
	for keyVal in data.keys():
		CVdict[keyVal]=np.std(data[keyVal])
	# iterate over pathways and calculate scores for pathway
	for pathway in pathList:
		print(pathway[0])
		# print out graphs with importance scores, rules, and relative abundances
		for RAval, comparator in zip(RAvals, comparisonStrings):
			outputGraphs(pathway, RAval, comparator, pathway[1])
		z_scores=[]
		# iterate over comparisons for each pathway and calculate z score
		for RAval in RAvals:
			z_scores.append(scorePathway(RAval,pathway[1], CVdict))
		pvals=scipy.stats.norm.sf(z_scores) # calculate p value
		# store p values
		tempdict={'pathway':pathDict[pathway[0][3:]],'code':pathway[0][3:] }
		for i in range(len(comparisonStrings)):
			tempdict[comparisonStrings[i]]=-math.log(pvals[i],10)
		csvmaker.append(tempdict)
	# output data
	df=pd.DataFrame(csvmaker,columns=['pathway','code'].extend(comparisonStrings))
	df.to_csv(path_or_buf='pvalues.csv')

# do pathway analysis! store in one folder for each comparison
def analyze_pathways_raven(diffName, matrixName, dataName, delmited):
	pathList=findPathwayList_justRavenPathway() # identify pathways under consideration
	data, colNums= readFpkm(dataName,delmited)	# read in fpkm data
	diffs= readDiffs(diffName,delmited) # read in difference to be considered
	groups=readMatrix(matrixName,delmited, colNums) # read design matrix

	# create an index of relative activities for all comparisons
	csvmaker=[]
	RAvals=[]
	comparisonStrings=[]
	for comparison in diffs:
		comparisonStrings.append(comparison[0]+'-'+comparison[1])
		RAvals.append(makeRA(data,comparison,groups))
		if not os.path.exists(comparison[0]+'-'+comparison[1]):
			os.makedirs(comparison[0]+'-'+comparison[1])
	CVdict={}
	for keyVal in data.keys():
		CVdict[keyVal]=np.std(data[keyVal])
	# iterate over pathways and calculate scores for pathway
	for pathway in pathList:
		print(pathway[0])
		# print out graphs with importance scores, rules, and relative abundances
		for RAval, comparator in zip(RAvals, comparisonStrings):
			outputGraphs(pathway, RAval, comparator, pathway[1], pathway[4])


# calculate z score for a given pathway
def scorePathway(RAs,pathImportances, CVdict):
	score=0
	allNodes=RAs.keys()
	for node in pathImportances:
		score+=abs(RAs[node])*math.log(pathImportances[node],2)*CVdict[node]
	# print(score)
	# print('Relative abundance mean difference: '+ str(np.mean([abs(RAs[value]) for value in allNodes])))
	randomScores=[]
	for i in range(1000):
		tempscore=0
		for node in pathImportances:
			tempscore+=abs(RAs[allNodes[randint(0,len(allNodes)-1)]])*math.log(pathImportances[node],2)*CVdict[node]
		randomScores.append(tempscore)
	meaner=np.mean(randomScores)
	stdev=np.std(randomScores)
	zscore=(score-meaner)/stdev
	return zscore

if __name__ == '__main__':
	import time
	start_time = time.time()
	
	# load arguments from user
	parser = argparse.ArgumentParser(prog='BONITA') 
	parser.set_defaults(sep=',')
	parser.add_argument("-sep", "--sep", metavar="seperator", help="How are columns in datafile specified")	
	parser.add_argument("-t", action='store_const',const='\t', dest="sep",help="Tab delimited?")	
	parser.add_argument("data")

	parser.add_argument("matrix")
	parser.add_argument("diffName")
	
	results = parser.parse_args()
	matrixName= results.matrix
	# run pathway analysis
	analyze_pathways_raven(results.diffName, matrixName, results.data, results.sep)
	print("--- %s seconds ---" % (time.time() - start_time))

