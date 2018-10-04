# import python packages needed in this file
import pickle
import networkx as nx
import operator
import argparse as argparse

from experiments import transformTest

if __name__ == '__main__':
	# use the command line to input the file that you want
	#Is this right?
	#python experiments.py graphname(to be given) binData\kMeans.bin iterNum(to be given)

	parser = argparse.ArgumentParser()
	parser.add_argument("graph")
	parser.add_argument("datafile")
	parser.add_argument("iterNum")

	results = parser.parse_args()
	graphName=results.graph
	fileName=results.datafile
	iterNum=int(results.iterNum)
	
	name=graphName[:-8]+fileName+'_'+results.iterNum ##

	# insert your code to get a ssDict here.
	
	# use fileName to read in the discretized file... you want to standardize the format
	#ssDict= SOMETHING
	
	# edit the below name to include the pertinent piece of the filename of discretized data
	#outname=graphName[:-8]+'_'+results.iterNum
	outname=graphName[:-8]+fileName+'_'+results.iterNum
	print(name) #
	graph = nx.read_gpickle(graphName)
	print(len(graph.nodes()))
	transformTest(graph,outname,fileName)
