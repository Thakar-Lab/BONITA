# import necessary modules
import argparse as argparse
import operator
import networkx as nx
import pickle
from ctypes import *
from simulation import paramClass, modelClass, NPsync
from utils import genInitValueList, setupEmptyKOKI, writeModel
from GA import GAsearchModel, localSearch
from analysis_accuracy import modelHolder
from pathway_analysis_score_nodes import *
import subprocess
import time
import glob
from random import randint

if __name__ == '__main__':

	start_time = time.time() # start timing this file
	GAmodels=glob.glob("*_GAmodel.pickle") #open GA replicates
	
	for GAmodel in GAmodels: # loop over models from GA
		storeModel1, bruteOut1, knockoutLists, knockinLists = pickle.load(open(GAmodel, "rb")) # open GA model
		model1 = modelHolder(storeModel1) # put GA model in a modelHolder object so properties are easy to acccess
		equivs, individual, devs=[] , [], [] # initialize arrays
		modelNumber=GAmodel.split('_')[3] # identify the number for the loaded GA replicate we are assessing
		ruleLength=0 # holder to determine the length of the individual... this is a sanity check
		for node in model1.nodeList: # loop over nodes
			filename='pickles/temp_series1_net_'+str(modelNumber)+'_'+node+'_local1.pickle' # find file with output of the local search
			output=pickle.load(open(filename, "rb")) # load the output of local search
			individual.extend(output[0]) # extract individual from results
			equivs.append(output[1]) # extract ERS from results
			devs.append(output[2]) # extract final errors from results
			ruleLength+=len(equivs) # count up the length of the individual
		pickle.dump( [storeModel1, knockoutLists, knockinLists, individual, equivs ], open( 'IS_'+str(modelNumber)+"_setup.pickle", "wb" ) ) # save a combined file for each GAmodel
		# code below generates a list of random individuals drawn from the ERS for use to do IS calculation
		individualList=[]
		for j in range(50):
			newindividual=[]
			for i in range(len(equivs)):
				newindividual.extend(equivs[i][randint(0,len(equivs[i])-1)])
			individualList.append(newindividual)
			print(newindividual[0:15])
		pickle.dump(individualList, open( 'IS_individuals_'+str(modelNumber)+"_setup.pickle", "wb" ) )
		# print out the length of the individual calculated from combining and the total time it took to perform these steps
		print(ruleLength)
		print("--- %s seconds ---" % (time.time() - start_time))

