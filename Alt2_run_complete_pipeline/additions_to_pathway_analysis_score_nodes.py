# import necessary modules
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

def runAllNodes(GAmodel, newSSS, nodeList, name):
    #def runAllNodes(nodeList, indy, newSSS, params, model,KOlist, KIlist):
    #all arguments except nodeList are names of pickle files
    for node in list(nodeList):
        shellHandle=open(str(name)+"_"+str(node)+"_localSearch.sh", "w+")
		#standard
        slurmCommands=str("#!/bin/sh\n#SBATCH --partition=standard\n#SBATCH -J "+name+"\n#SBATCH -o "+str(name)+"_"+str(node)+".log\n#SBATCH -t 10:00:00\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --mem=120G\nmodule load intelpython/2.7.12\nmake\npython parallel_local_search.py ")+str(node)+str(" ")+str(GAmodel)+str(" ")+str(newSSS) #+str(node)+str(model)+str(indy)+str(newSSS)+str(params)+str(KOlist)+str(KIlist))
        #debug
		#slurmCommands=str("#!/bin/sh\n#SBATCH --partition=debug\n#SBATCH -J "+name+"\n#SBATCH -o "+str(name)+"_"+str(node)+".log\n#SBATCH -t 1:00:00\n#SBATCH -n 1\n#SBATCH -c 1\n#SBATCH --mem=50G\nmodule load intelpython/2.7.12\nmake\npython parallel_local_search.py ")+str(node)+str(" ")+str(GAmodel)+str(" ")+str(newSSS) #+str(node)+str(model)+str(indy)+str(newSSS)+str(params)+str(KOlist)+str(KIlist))
        shellHandle.write(slurmCommands)
        shellHandle.close()
        shellCommand=str(name)+"_"+str(node)+"_localSearch.sh"
        print(str(shellCommand))
        p = subprocess.Popen(['sbatch', shellCommand])

def parallelLocalSearch(GAmodel, newSSS, nodeList, name):

    runAllNodes(GAmodel=GAmodel,newSSS=newSSS, nodeList=nodeList, name=name)


if __name__ == '__main__':

	start_time = time.time()

	# read in arguments from shell scripts
	parser = argparse.ArgumentParser()
	parser.add_argument("graph")
	parser.add_argument("iterNum")
	results = parser.parse_args()
	graphName=results.graph
	print(graphName)
	iterNum=int(results.iterNum)
	name=graphName[:-8]+'_'+results.iterNum
	graph = nx.read_gpickle(str(graphName[:-8])+".gpickle")

	# read in C function to run simulations
	updateBooler=cdll.LoadLibrary('./simulator.so')
	boolC=updateBooler.syncBool

	# load data
	sampleList=pickle.Unpickler(open( graphName[:-8]+'_sss.pickle', "rb" )).load()

	# set up parameters of run, model
	params=paramClass()
	model=modelClass(graph,sampleList, False)
	model.updateCpointers()

	storeModel=[(model.size), list(model.nodeList), list(model.individualParse), list(model.andNodeList) , list(model.andNodeInvertList), list(model.andLenList),	list(model.nodeList), dict(model.nodeDict), list(model.initValueList)]

	# put lack of KOs, initial values into correct format
	knockoutLists, knockinLists = setupEmptyKOKI(len(sampleList))
	newInitValueList = genInitValueList(sampleList,model)
	model.initValueList=newInitValueList

	# find rules by doing GA then local search
	model1, dev, bruteOut = GAsearchModel(model, sampleList, params, knockoutLists, knockinLists, name, boolC) # run GA

    # dump GA model
	storeModel1 = [(model1.size), list(model1.nodeList), list(model1.individualParse), list(model1.andNodeList) , list(model1.andNodeInvertList), list(model1.andLenList),	list(model1.nodeList), dict(model1.nodeDict), list(model1.initValueList)]
	pickle.dump([storeModel1, bruteOut, knockoutLists, knockinLists], open( name+"_GAmodel.pickle", "wb" ) ) # output rules

	#Run parallel local search

	GAmodel =pickle.load(open(name+"_GAmodel.pickle", "rb"))
	nodeList = GAmodel[0][1]

	parallelLocalSearch(GAmodel=name+"_GAmodel.pickle", newSSS=graphName[:-8]+'_sss.pickle', nodeList=list(nodeList), name=name)
