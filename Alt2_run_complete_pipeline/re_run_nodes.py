import argparse as argparse
import operator
import networkx as nx
import pickle
from ctypes import *
from simulation import paramClass, modelClass, NPsync
from utils import genInitValueList, setupEmptyKOKI, writeModel
from GA import GAsearchModel, localSearch
from pathway_analysis_score_nodes import *
import subprocess
import time
import os as os

# this file uses the graph file name and iterNum passed to find all that GA iteration's local search files and re run jobs which did not complete
def runAllNodes(GAmodel, newSSS, nodeList, name):
	counter = 1 # count number of jobs started 
	for node in list(nodeList):
		# set up and run a job for any node that doesn't have results
		if not os.path.isfile(str(name)+"_"+str(node)+"_local1.pickle"):
			shellHandle=open(str(name)+"_"+str(node)+"_localSearch.sh", "w+")
			slurmCommands=str("#!/bin/sh\n#SBATCH --partition=standard\n#SBATCH -J "+name+"\n#SBATCH -o "+str(name)+"_"+str(node)+".log\n#SBATCH -t 1:00:00\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --mem=5G\nmodule load intelpython/2.7.12\nsleep 1\nmake\npython parallel_local_search.py ")+str(node)+str(" ")+str(GAmodel)+str(" ")+str(newSSS) #+str(node)+str(model)+str(indy)+str(newSSS)+str(params)+str(KOlist)+str(KIlist))
			shellHandle.write(slurmCommands)
			shellHandle.close()
			shellCommand=str(name)+"_"+str(node)+"_localSearch.sh"
			print(str(shellCommand))
			p = subprocess.Popen(['sbatch', shellCommand])
			counter= counter+1
		else: 
			print(node)
		if counter > 350: # limit total number of jobs to 350*5= 1750 to stay under BlueHive limit of 2000
			print('completed run on node: '+node)
			break

if __name__ == '__main__':
	# read in arguments from shell scripts
	parser = argparse.ArgumentParser()
	parser.add_argument("graph")
	parser.add_argument("iterNum")
	results = parser.parse_args()
	graphName=results.graph
	print(graphName)
	iterNum=int(results.iterNum)
	name=graphName[:-8]+'_'+results.iterNum
	GAmodel =pickle.load(open(name+"_GAmodel.pickle", "rb"))
	nodeList = GAmodel[0][1]

	#parallelLocalSearch(GAmodel=name+"_GAmodel.pickle", newSSS=graphName[:-8]+'_sss.pickle', nodeList=list(model1.nodeList), name=name)
	runAllNodes(GAmodel=name+"_GAmodel.pickle", newSSS=graphName[:-8]+'_sss.pickle', nodeList=list(nodeList), name=name)