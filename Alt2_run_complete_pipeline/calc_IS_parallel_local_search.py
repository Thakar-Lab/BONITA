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
import subprocess
import time
import glob
from random import randint

if __name__ == '__main__':

	start_time = time.time()
	#name = "temp"

	parser = argparse.ArgumentParser()
	parser.add_argument("model")
	parser.add_argument("individuals")
	parser.add_argument("iterNum")
	results = parser.parse_args()
	iterNum=int(results.iterNum)
	storeModel1, knockoutLists, knockinLists, individual, equivs = pickle.load(open(results.model, "rb"))
	individuals = pickle.load(open(results.individuals, "rb"))
	
	individual=individuals[iterNum-1]
	model1 = modelHolder(storeModel1)

	# read in C function to run simulations
	updateBooler=cdll.LoadLibrary('./simulator.so')
	boolC=updateBooler.syncBool

	# load data
	sampleList=pickle.Unpickler(open('temp_series1_net_sss.pickle', "rb" )).load()

	# set up parameters of run, model
	params=paramClass()
	model1.modelHolder_updateCpointers()

	# # calculate importance scores and output
	scores1=calcImportance(individual,params,model1, sampleList,knockoutLists, knockinLists, boolC)
	pickle.dump(scores1, open(results.model[:-11]+'_'+str(iterNum)+"_scores1.pickle", "wb"))

	# # write rules
	# with open(name+"_rules.txt", "w") as text_file:
	# 	text_file.write(writeModel(bruteOut1, model1))
	print("--- %s seconds ---" % (time.time() - start_time))

