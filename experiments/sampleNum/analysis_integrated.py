
# import python packages
import matplotlib.pyplot as plt
import numpy as numpy
import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl
import seaborn as sns
from sets import Set
import os as os
from scipy.stats import variation
import networkx as nx
# import other parts of our code
import utils as utils

# find the end of a node in the bistring
def findEnds(model, node, indiv):
	if node==len(model.nodeList)-1:
		end1=len(indiv)
	else:
		end1=model.individualParse[node+1]
	start1=model.individualParse[node]
	return start1, end1

# find the incoming edges to each 'and' connection for a given node
def findInEdges(model,node):
	inEdges=[]
	for lister in model.andNodeList[node]:
		tempTup=tuple(lister)
		inEdges.append(set(tempTup))
	return inEdges

# find the simplest form of a rule
def simplifyRule(rule, inEdges):
	for i in range(len(rule)):
		if rule[i]==1:
			for j in range(len(rule)):
				if rule[j]==1 and not i==j:
					if inEdges[i].issubset(inEdges[j]):
						rule[j]=0
	return rule


# complete and save a plot
def finishPlot(xlabeler, ylabeler, plotname):
	plt.ylabel(ylabeler, labelpad=1)
	plt.xlabel(xlabeler, labelpad=.5)
	sns.despine()
	sns.set_style("ticks")
	plt.rc('font', size=8)          # controls default text sizes
	plt.rc('axes', titlesize=8)     # fontsize of the axes title
	plt.rc('axes', labelsize=8)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=8)    # fontsize of the tick labels
	plt.rc('legend', fontsize=8)    # legend fontsize
	plt.rc('figure', titlesize=8)  # fontsize of the figure title
	plt.savefig(plotname, bbox_inches='tight', transparent=True, dpi=600)
	plt.savefig(plotname[:-3]+'svg', bbox_inches='tight', transparent=True, dpi=600)
	plt.clf()
# set up figure options and convert the dict to a pandas dataframe
def setupPlot(dictlist):
	sns.set_style("ticks")
	# sns.set_context("paper")
	sns.set_palette(sns.color_palette("Greys_r", 6))	
	return pd.DataFrame(dictlist)

# make a swarm plot
def plotSwarm(dictlist, xval, yval, title, xlabeler, ylabeler, plotname):
	df= setupPlot(dictlist)
	ax=sns.swarmplot(data=df,x=xval, y=yval, color='black')
	finishPlot(xlabeler, ylabeler, plotname)

# make a scatter plot
def plotScatter(dictlist, xval, yval, title, xlabeler, ylabeler, colorVar, plotname):
	df= setupPlot(dictlist)
	if colorVar=='none':
		ax=sns.lmplot(data=df,x=xval, y=yval, fit_reg=True, )
	else:
		ax=sns.lmplot(data=df,x=xval, y=yval,hue=hue, fit_reg=False)
	finishPlot(xlabeler, ylabeler, plotname)

# make a scatterplot with histograms on sides
def plotHistScatter(dictlist, xval, yval, title, xlabeler, ylabeler, plotname,xlims):
	df= setupPlot(dictlist)
	ax=sns.jointplot(data=df,x=xval, y=yval, kind='kde')
	finishPlot(xlabeler, ylabeler, plotname)

# make a boxplot
def plotBar(dictlist, xval, yval, title, xlabeler, ylabeler, plotname, hueBool, huer):
	# sns.axes_style("ticks", rc={'axes.edgecolor': '.8','xtick.color': '.9', 'text.color': '.85', 'xtick.color': '.85'})
	# sns.set_context("talk")
	plt.rcParams.update({'font.size': 8})
	plt.rc('font', size=8)          # controls default text sizes
	plt.rc('axes', titlesize=8)     # fontsize of the axes title
	plt.rc('axes', labelsize=8)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=8)    # fontsize of the tick labels
	plt.rc('legend', fontsize=8)    # legend fontsize
	plt.rc('figure', titlesize=8)  # fontsize of the figure title
	df= setupPlot(dictlist)
	sns.set_style("ticks")
	if not hueBool:
		fig, ax = plt.subplots(figsize=[2.6,2])
		ax=sns.barplot(data=df,ax=ax, x=xval, y=yval, palette=sns.color_palette("Greys", 6), ci='sd', capsize=.05)
	else:
		fig, ax = plt.subplots(figsize=[2.6,2])
		ax=sns.barplot(data=df,ax=ax, x=xval, y=yval, palette=sns.color_palette("Greys_r", 6), hue=huer, ci='sd', capsize=.05, errwidth=.5)
		sns.set_style("ticks")
		plt.rc('font', size=8)          # controls default text sizes
		plt.rc('axes', titlesize=8)     # fontsize of the axes title
		plt.rc('axes', labelsize=8)    # fontsize of the x and y labels
		plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
		plt.rc('ytick', labelsize=8)    # fontsize of the tick labels
		plt.rc('legend', fontsize=8)    # legend fontsize
		plt.rc('figure', titlesize=8)  # fontsize of the figure title
		plt.legend( loc="upper left", handletextpad=0.15,bbox_to_anchor=[-.05, -.2, 1, .1], bbox_transform=plt.gcf().transFigure, framealpha=0.0,edgecolor='black', ncol=4,mode="expand" )
		plt.rc('legend', fontsize=8)    # legend fontsize
		plt.setp(ax.patches, linewidth=0)
		plt.xticks(rotation=70)
	finishPlot(xlabeler, ylabeler, plotname)

# make a scatterplot with regression line
def plotScatterReg(dictlist, xval, yval, title, xlabeler, ylabeler, plotname, hueBool, huer):
	df= setupPlot(dictlist)
	if not hueBool:
		ax=sns.regplot(data=df,x=xval, y=yval, color='black', fit_reg=False)
	else:
		ax=sns.regplot(data=df,x=xval, y=yval, fit_reg=False)
		plt.legend(bbox_to_anchor=(1.1, .8), bbox_transform=plt.gcf().transFigure)
	finishPlot(xlabeler, ylabeler, plotname)

# find characteristics of upstream nodes from each particular node
def findUpstreamChars(stem, nodePPV, nodeSens, nodeRuleTruth, maxRT, model):
	graph = nx.read_gpickle('gpickles/'+stem[8:-1]+'.gpickle')
	inDegree=[] # in degree of each node
	upstreamPPV, upstreamSens, upstreamRT, upstreamMaxRT=[],[],[], [] # accuracy statistics on average across upstream nodes for each node in network
	# loop through network
	for node in model.nodeList:
		tempPPV, tempSens, tempRT, tempMaxRT=[], [], [], [] # accuracy statistics for each upstream node
		templist= graph.predecessors(node) # generate list of upstream nodes
		inDegree.append(len(templist)) # save in-degree
		# save accuracy statistics for each upstream node
		for tempNode in templist:
			if (nodePPV[model.nodeDict[tempNode]] <100):
				tempPPV.append(nodePPV[model.nodeDict[tempNode]])
				tempSens.append(nodeSens[model.nodeDict[tempNode]])
				tempRT.append(nodeRuleTruth[model.nodeDict[tempNode]])
				tempMaxRT.append(maxRT[model.nodeDict[tempNode]])
		# average upstream node accuracy statistics
		if len(tempPPV)>0:
			upstreamPPV.append(numpy.mean(tempPPV))
			upstreamSens.append(numpy.mean(tempSens))
			upstreamRT.append(numpy.mean(tempRT))
			upstreamMaxRT.append(numpy.mean(tempMaxRT))
		# return 1. if a source node
		else:
			upstreamPPV.append(1.)
			upstreamSens.append(1.)
			upstreamRT.append(1.)
			upstreamMaxRT.append(1.)
	return upstreamPPV, upstreamSens, upstreamRT, upstreamMaxRT, inDegree

def compareIndividualsNodeWise(truthList, testList, model1s, model2s,covs, equivs):

	modeler=model1s[0]
	SDs=[0. for q in truthList]
	nodeSDs=[]
	edgeSens, inDegrees, edgePPVs=[], [], []
	inCoV=[]
	TPsum, TNsum, FPsum, FNsum=0, 0,0,0
	for node in range(0,len(modeler.nodeList)):
		tempSD=0.
		FP, TP, TN, FN=0, 0,0,0
		# simplify rules at the node and find the edge-wise PPV, sens, and SDs
		inCovTemper=[]
		for k in range(len(truthList)):
			inCovtemp=[]
			# find start and end of this node in each model
			start1,end1 = findEnds(model1s[k], node, truthList[k])
			start2,end2 = findEnds(model2s[k], node, testList[k])			
			
			# find the shadow and nodes for each model
			truthInEdges=findInEdges(model1s[k],node)
			testInEdges=findInEdges(model2s[k],node)

			# find the bitstring for just this node
			truth= truthList[k][start1:end1]
			test= testList[k][start2:end2]
			
			# simplify ground truth and recovered rules
			truth= simplifyRule(truth, truthInEdges)
			test= simplifyRule(test, testInEdges)
			
			# edit overall rule list with simplified rules
			testList[k][start2:end2] = test
			truthList[k][start1:end1] = truth
		
			# find SD, PPV, etc....
			truthSet=Set([]) # edges in correct rule
			testSet=Set([]) # edges in rule found
			baseSet=Set([]) # edges possible across all rules
			
			# find edges in true rule (and edges possible), average incoming coefficient of variation
			for i in range(0,len(truth)):
				if truth[i]==1:
					for nodeToAdd in model1s[k].andNodeList[node][i]:
						truthSet.add(nodeToAdd)
						inCovtemp.append(covs[k][node])
				for nodeToAdd in model1s[k].andNodeList[node][i]:
					baseSet.add(nodeToAdd)
			# find edges in test (recovered) rule
			for i in range(0,len(test)):
				if test[i]==1:
					for nodeToAdd in model2s[k].andNodeList[node][i]:
						testSet.add(nodeToAdd)
			# find structural distance at this node. 
			SDs[k]=SDs[k]+len(truthSet.difference(testSet))+len(testSet.difference(truthSet))
			tempSD=tempSD+len(truthSet.difference(testSet))+len(testSet.difference(truthSet))
			# save edge-wise statistics for this node
			FP+=1.*len(testSet.difference(truthSet))
			TP+=1.*len(testSet.intersection(truthSet))
			FN+=1.*len(truthSet.difference(testSet))
			inCovTemper.append(numpy.mean(inCovtemp))
		# calculae and save overall edge-wise statistics
		if (TP+FN)>0:
			sensitivity=1.*TP/(TP+FN)
		else:
			sensitivity=100
		if TP+FP>0:
			PPV=1.*TP/(TP+FP)
		else:
			PPV=100
		nodeSDs.append(tempSD/len(truthList))
		edgeSens.append(sensitivity)
		edgePPVs.append(PPV)
		TPsum+=TP
		FNsum+=FN
		FPsum+=FP
		inDegrees.append(len(baseSet))
		inCoV.append(numpy.mean(inCovTemper))
	if (TPsum+FNsum)>0:
		edgeSens=1.*TPsum/(TPsum+FNsum)
	else:
		edgeSens=100
	if (FPsum+TPsum)>0:
		edgePPV= 1.*TPsum/(FPsum+TPsum)
	else:
		edgePPV=100


	nodeSens=[] # sensitivity by node
	nodePPV=[] # PPV by node 
	nodeRTs=[] # rules true by node
	nodePsens=[]
	nodepPPV=[]
	nodelister=model1s[0].nodeList # gives a node List (should be the same across all trials in a network...)
	sampleRTs=[[] for item in truthList] # Rules True for each trial
	samplePPVs=[[] for item in truthList] # PPV for each trial
	sampleSenss=[[] for item in truthList] # Sens for each trial
	equivRTsens=[[] for item in truthList] # RT sensitivity of equivalents for each trial
	equivSens=[[] for item in truthList] # sensitivity for equivalents for each trial
	equivNodeRTsens=[]
	equivNodeSens=[]

	# iterate over all nodes in the network
	for node in range(len(nodelister)):
		rtTemp=[] # stores rules true for this node across all networks
		ones=[] # stores the number of false negative and rules
		zeros=[] # stores the number of  correct and rules
		negones=[] # stores the number of false positive and rules 
		equivOnes=[] # stores the min number of false negatives across equivs 
		equivZeros=[] # stores the max correct across equivs
		equivNegOnes=[] # stores the min false positives across equivs
		sumindividual=[] # total number true positive and rules
		equivRTsensNode=[]
		equivSensNode=[]
		
		#loop over individuals provided and calculate sens, PPV, rules true
		for i in range(len(truthList)):
			
			# find start and end of this node in each model
			start1,end1 = findEnds(model1s[i], node, truthList[i])
			start2,end2 = findEnds(model2s[i], node, testList[i])			
			
			# find the values for just this node
			truth= truthList[i][start1:end1]
			test= testList[i][start2:end2]

			# set up empty lists for ands, edges, and the shadow and nodes associated with this node in each model
			truthAnds=[]
			testAnds=[]
			
			# get the set of all shadow and nodes that are actually used in each rule
			for j in range(len(model1s[i].andNodeList[node])):
				if truth[j]>0:
					truthAnds.append(tuple(model1s[i].andNodeList[node][j]))
			for j in range(len(model2s[i].andNodeList[node])):
				if test[j]>0:
					testAnds.append(tuple(model2s[i].andNodeList[node][j]))
			truthAnd=tuple(truthAnds)
			truthAnd=set(truthAnd)
			testAnd=set(tuple(testAnds))
			
			# get the set of all shadow and nodes used in each equivalent rule
			equivAnds=[]
			# print(equivs[i])
			for test1 in equivs[i][node]:
				tempEquivAnd=[]
				for j in range(len(model2s[i].andNodeList[node])):
					if test1[j]>0:
						tempEquivAnd.append(tuple(model2s[i].andNodeList[node][j]))
				testAnd1=set(tuple(tempEquivAnd))
				equivAnds.append(testAnd1)
			RTequiv=0.
			possibilityOnes=[]
			possibilityZeros=[]
			possibilityZNetones=[]
			for testAnder1 in equivAnds:	
				if (truthAnd==testAnder1):
					RTequiv=1.
				possibilityOnes.append(len(truthAnd.difference(testAnd)))
				possibilityZeros.append(len(truthAnd.intersection(testAnd)))
				possibilityZNetones.append(len(testAnd.difference(truthAnd)))
			# append results for this trial to all results
			maxpossibilityZeros=max(possibilityZeros)
			minpossiblityOnes=min(possibilityOnes)
			minpossibilityNegOnes=min(possibilityZNetones)
			equivOnes.append(minpossiblityOnes)
			equivZeros.append(maxpossibilityZeros)
			equivNegOnes.append(minpossibilityNegOnes)
			equivRTsensNode.append(RTequiv)
			equivRTsens[i].append(RTequiv)
			# calculate true positives, false positives, false negatives, and total slots for this node, trial and save
			onetemp=len(truthAnd.difference(testAnd))
			zerotemp=len(truthAnd.intersection(testAnd))
			negonetemp=len(testAnd.difference(truthAnd))
			sumindtemp=len(truthAnd)
			ones.append(onetemp)
			zeros.append(zerotemp)
			negones.append(negonetemp)
			sumindividual.append(sumindtemp)			
			# add Rules true first sample-wise then node-wise
			if len(model1s[i].andNodeList[node])>1:
				if (truthAnd==testAnd):
					sampleRTs[i].append(1.)
				else:
					sampleRTs[i].append(0.)
				if (sumindtemp-onetemp+negonetemp)>0:
					samplePPVs[i].append(1.*(sumindtemp-onetemp)/(sumindtemp-onetemp+negonetemp))
				else:
					samplePPVs[i].append(100)
				if (sumindividual[i])>0:
					sampleSenss[i].append(1.*(sumindtemp-onetemp)/(sumindtemp))
				else:
					sampleSenss[i].append(100)
			if (truthAnd==testAnd):
				rtTemp.append(1.)
			else:
				rtTemp.append(0.)

		nodeRTs.append(numpy.mean(rtTemp)) # node-wise Rules true added
		equivNodeRTsens.append(numpy.mean(equivRTsensNode))

		# calculate sensitivity for the node
		temp=[100 if sumindividual[i]==0 else 1.*(sumindividual[i]-ones[i])/(sumindividual[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			sensitivity=100
		else:
			sensitivity=(1.*numpy.sum(temp)/len(temp))

		# calculate max sensitivity for the node
		temp=[100 if sumindividual[i]==0 else 1.*(sumindividual[i]-equivOnes[i])/(sumindividual[i]) for i in range(0,len(equivOnes))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			psensitivity=100
		else:
			psensitivity=(1.*numpy.sum(temp)/len(temp))
		nodePsens.append(psensitivity)

		# calculate PPV for the node
		temp=[100 if (sumindividual[i]-ones[i]+negones[i])==0 else 1.*(sumindividual[i]-ones[i])/(sumindividual[i]-ones[i]+negones[i]) for i in range(0,len(ones))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			PPV=100
		else:
			PPV=(1.*numpy.sum(temp)/len(temp))

		# calculate PPV for the node
		temp=[100 if (sumindividual[i]-equivOnes[i]+equivNegOnes[i])==0 else 1.*(sumindividual[i]-equivOnes[i])/(sumindividual[i]-equivOnes[i]+equivNegOnes[i]) for i in range(0,len(equivOnes))]
		temp=filter(lambda a: a != 100, temp)
		if len(temp)==0:
			pPPV=100
		else:
			pPPV=(1.*numpy.sum(temp)/len(temp))
		nodepPPV.append(pPPV)

		# add to list of sensitivity and PPV by 
		nodeSens.append(sensitivity)
		nodePPV.append(PPV)
	sampleEquivRT=[1.*numpy.mean(filter(lambda a: a != 100, sampler)) for sampler in equivRTsens] # Rules True for each trial
	sampleRT=[1.*numpy.mean(filter(lambda a: a != 100, sampler)) for sampler in sampleRTs] # Rules True for each trial
	samplePPV=[1.*numpy.mean(filter(lambda a: a != 100, sampler)) for sampler in samplePPVs] # PPV for each trial
	sampleSens=[1.*numpy.mean(filter(lambda a: a != 100, sampler)) for sampler in sampleSenss] # Sens for each trial
	return sampleEquivRT, equivNodeRTsens, nodePsens,nodepPPV,sampleSens, samplePPV, nodeSens, nodePPV, sampleRT, nodeRTs, edgeSens, edgePPV, SDs, nodeSDs, len(modeler.nodeList), inDegrees, inCoV

# class to mimick model class for easy access
class modelHolder:
	def __init__(self,valueList):
		[ size, nodeList, individualParse, andNodeList, andNodeInvertList, andLenList,nodeList, nodeDict, initValueList]=valueList
		self.size=size
		self.individualParse=individualParse
		self.andNodeList=andNodeList
		self.andNodeInvertList=andNodeInvertList
		self.andLenList=andLenList
		self.nodeList=nodeList
		self.nodeDict=nodeDict
		self.initValueList=initValueList

# perform analysis for a particular graph and trial: uploads data, funs comparisons, and returns
def analyzeGraph(stem, replicates):
	model1s, model2s,truthList, testList, varLists, equivalents, devs=[], [],[], [],[], [], [] # set up variables
	# iterate over replicates
	for i in range(1,replicates+1):
		outputList=pickle.Unpickler(open( stem+str(i)+'_local1.pickle', "rb" )).load() # load in main data
		[truth, test, truthmodel, testmodel, equivs, dev]=outputList
		equivalents.append(equivs)
		model1s.append(modelHolder(truthmodel))
		model2s.append(modelHolder(testmodel))
		truthList.append(truth)
		testList.append(test)
		inputList=pickle.Unpickler(open( stem+str(i)+'_input.pickle', "rb" )).load() # load in generated data
		varLists.append([numpy.std([inputList[k][j] for k in range(len(inputList))]) for j in range(len(inputList[0])) ]) # put variances in correct order
	return model1s[0],compareIndividualsNodeWise(truthList, testList, model1s,model2s, varLists, equivalents), varLists

# analyze omics noise experiment
def analyzeOmicsNoiseExperiment(codes, end):
	noiseLabels=[1,2,3,4,5,6,7,8]		
	noiseVals=[.01,.02,.05,.1,.5,.75,1,2]
	analyzeNoiseExperiment(codes, end, noiseLabels, noiseVals)

# analyze sample number experiment
def analyzeSamplesExperiment(codes, end):
	noiseLabels=[1,2,3,4,5,6]		
	noiseVals=[2,3,4,5,10,15]
	analyzeNoiseExperiment(codes, end, noiseLabels, noiseVals)

# analyze experiment varying noise- in RPKN, omics levels etc. 
def analyzeNoiseExperiment(codes, end, noiseLabels, noiseVals):
	sensspeclist=[]
	for code in codes:
		print(code)		
		#FOR EACH RUN WITH THAT GPICKLE
		for noiseLabel, noiseVal in zip(noiseLabels,noiseVals):
			truthmodel,result, varList=analyzeGraph('pickles/'+code+str(noiseLabel)+'_',10) # analyze accuracies for this network and noise level
			# put results into a dataframe
			sampleEquivRT, equivNodeRTsens,nodePsens,nodepPPV, tempSens, tempPPV, nodeSens, nodePPV, ruleTruth, nodeRuleTruth, edgeSens, edgePPV, SD, nodeSDs, nodeNum, edgeDegree, inCoV = result
			for i in range(len(tempSens)):
				sensspeclist.append({'Percent Rules True':100.*sampleEquivRT[i],'Node_Num':nodeNum,'Sensitivity':tempSens[i],'Equiv RT Sens':sampleEquivRT[i],'PPV':tempPPV[i], 'Proportion Rules True': ruleTruth[i],'GA': 'no Adapt','Noise':noiseVal})
	
	# make plots
	df=pd.DataFrame(sensspeclist)
	groupedvalues=df.groupby('Noise')
	for k, gp in groupedvalues:
		print(k)
		print(gp.groupby('Node_Num').mean().reset_index())
	plotBar(sensspeclist,'Node_Num', 'Percent Rules True','Percent Rules True by Node Number', 'Node Number', 'Percent Rules True','percent_RT_equiv'+end+'.png', False, 'Noise')
	plotBar(sensspeclist,'Node_Num', 'Percent Rules True','Percent Rules True by Node Number', 'Node Number', 'Percent Rules True','Noise_percent_RT'+end+'.png', True, 'Noise')


	# plotBar(sensspeclist,'Node_Num', 'Equiv RT Sens','Rules True by Node Number', 'Node Number', 'Proportion Rules True','RT_equiv'+end+'.png', True, 'Noise')
	# plotBar(sensspeclist,'Node_Num', 'Proportion Rules True','Rules True by Node Number', 'Node Number', 'Proportion Rules True','Rules_True'+end+'.png', True, 'Noise')
	# plotBar(sensspeclist,'Noise', 'Equiv RT Sens','Rules True by Sample Number', 'Sample Number', 'Proportion Rules True','Noise_Equiv_RT_'+end+'.png', False, 'Noise')
	# plotBar(sensspeclist,'Noise', 'Proportion Rules True','Rules True by Sample Number', 'Sample Number', 'Proportion Rules True','Noise_RT_'+end+'.png', False, 'Noise')
# analyze accuracy at a granular level when no noise is being used 
def analyzeAccuracy(codes, end):
	degree2dictlist, degree3dictlist, SDs, nodeSDS, nodeNumbers,nodeNumlistcut, nodeNumlistuncut, dfdictlist, sensspeclist=[], [], [], [],[],[],[], [], []
	nodewiseList, nodewiseList2, nodewiseList3=[],[],[] # set up dictionaries for all nodes, in-degree 2 nodes, and in-degree 3 nodes
	sensspeclist=[] # set up dictionaries for networks as a whole
	# for each gpickle
	for code in codes:
		print(code)		
		nodeLookup={}
		tempsensitivities=[[],[],[],[]]
		tempspecificities=[[],[],[],[]]

		ruleTruthtots=[]
		#FOR EACH RUN WITH THAT GPICKLE
		truthmodel,result, varList=analyzeGraph('pickles/'+code+'_',25)
		sampleEquivRT, equivNodeRTsens,nodePsens,nodepPPV, tempSens, tempPPV, nodeSens, nodePPV, ruleTruth, nodeRuleTruth, edgeSens, edgePPV, SD, nodeSDs, nodeNum, edgeDegree, inCoV = result
		upstreamPPV, upstreamSens, upstreamRT, upstreamMaxRT, inDegree=findUpstreamChars('pickles/'+code+'_', nodePPV, nodeSens, nodeRuleTruth, equivNodeRTsens, truthmodel)
		for i in range(len(nodePPV)):
			cov=numpy.mean([varList[k][i] for k in range(25)])
			if not nodePPV[i]==100:
				nodewiseList.append({'equiv Sens': nodePsens[i],'equiv PPV': nodepPPV[i],'equiv RT':equivNodeRTsens[i], 'Node_Num':nodeNum,'Sensitivity':nodeSens[i],'PPV':nodePPV[i], 'Proportion Rules True': nodeRuleTruth[i],'CoV': cov,'in_CoV': inCoV[i], 'SD':nodeSDs[i], 'in_degree':inDegree[i], 'upstream Max RT': upstreamMaxRT[i], 'upstream RT': upstreamRT[i]})
				if edgeDegree[i]==3:
					nodewiseList3.append({'equiv Sens': nodePsens[i],'equiv PPV': nodepPPV[i],'equiv RT':equivNodeRTsens[i],'Node_Num':nodeNum,'Sensitivity':nodeSens[i],'PPV':nodePPV[i], 'Proportion Rules True': nodeRuleTruth[i],'CoV': cov,'in_CoV': inCoV[i], 'SD':nodeSDs[i], 'in_degree':inDegree[i], 'upstream Max RT': upstreamMaxRT[i], 'upstream RT': upstreamRT[i]})
				if edgeDegree[i]==2:
					nodewiseList2.append({'equiv Sens': nodePsens[i],'equiv PPV': nodepPPV[i],'equiv RT':equivNodeRTsens[i],'Node_Num':nodeNum,'Sensitivity':nodeSens[i],'PPV':nodePPV[i], 'Proportion Rules True': nodeRuleTruth[i],'CoV': cov,'in_CoV': inCoV[i], 'SD':nodeSDs[i], 'in_degree':inDegree[i], 'upstream Max RT': upstreamMaxRT[i], 'upstream RT': upstreamRT[i]})
		nodeNumbers.append(nodeNum)
		if len(ruleTruth)>0:
			ruleTruthtots.append(numpy.mean(ruleTruth))
		SDs.append(SD)
		nodeSensDict={}
		nodePPVDict={}
		nodeRTdict={}
		# assign what the var of interest is
		for i in range(len(tempSens)):
			sensspeclist.append({'Node_Num':nodeNum,'Sensitivity':tempSens[i],'Equiv RT Sens':sampleEquivRT[i],'PPV':tempPPV[i], 'Proportion Rules True': ruleTruth[i],'GA': 'no Adapt'})
		print(nodeNumbers[-1])	
	plotBar(sensspeclist,'Node_Num', 'Equiv RT Sens','Rules True by Node Number', 'Node Number', 'Proportion Rules True','Adapt_GA_Equiv_RT_sens'+end+'.png', False, 'Noise')
	plotBar(sensspeclist,'Node_Num', 'Proportion Rules True','Rules True by Node Number', 'Node Number', 'Proportion Rules True','Adapt_GA_Rules_True'+end+'.png', False, 'Noise')

	df=pd.DataFrame(sensspeclist)
	groupedvalues=df.groupby('Node_Num').mean().reset_index()
	print(groupedvalues)
	plotBar(sensspeclist,'Node_Num', 'Sensitivity','Sensitivity by Node Number', 'Node Number', 'Sensitivity','Adapt_GA_Sensitivity'+end+'.png', False, 'Noise')
	plotBar(sensspeclist,'Node_Num', 'PPV','PPV by Node Number', 'Node Number', 'PPV','Adapt_GA_PPV'+end+'.png', False, 'Noise')

if __name__ == '__main__':
	codes=[]
	for file in os.listdir("gpickles"):
		if file.endswith(".gpickle"):
			codes.append(file[:-8])
	print(codes)
	analyzeSamplesExperiment(codes,'')