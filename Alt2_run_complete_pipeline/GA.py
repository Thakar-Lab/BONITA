# import python packages
import pickle
import copy as copy
from deap import base, creator, gp, tools
from deap import algorithms as algo
from random import random, seed, shuffle, randint, sample, choice
import numpy as numpy
import operator
import math as math
from sets import Set
import gc as gc
from collections import defaultdict, deque
from itertools import chain
from operator import attrgetter, itemgetter
# import other pieces of our software
from simulation import paramClass, NP
from utils import findEnd, genRandBits, bitList, genInitValueList

# generates random bitstring with at least one value for each node
def genBits(model):
	# generate random bitlist
	startInd=list(genRandBits(model.size))
	counter=0
	# make sure bitlist isn't zero
	while numpy.sum(startInd)==0 and counter < 10000:
		startInd=list(genRandBits(model.size))
		counter+=1
	# go through nodes and make sure that there are 1-5 ones in the random list
	for node in range(0,len(model.nodeList)):
		end=findEnd(node,model)
		start=model.individualParse[node]
		if (end-start)>1:
			counter=0
			while numpy.sum(startInd[start:end])>5 and counter < 10000:
				chosen=math.floor(random()*(end-start))
				startInd[start+int(chosen)]=0
				counter+=1
			if numpy.sum(startInd[start:end])==0:
				chosen=math.floor(random()*(end-start))
				startInd[start+int(chosen)]=1
		elif (end-start)==1:
			startInd[start]=1
	return [copy.deepcopy(model),startInd]

# finds the lowest error individual in a population
def findPopBest(population):
	saveVal=-1
	minny=100000
	for i in range(len(population)):
		if numpy.sum(population[i].fitness.values)< minny:
			minny=numpy.sum(population[i].fitness.values)
			saveVal=i
	ultimate=population[saveVal]
	minvals=population[saveVal].fitness.values
	return minvals, ultimate[1], ultimate[0]

# executes two point crossover at node junctions
def cxTwoPointNode(ind1, ind2):
	# copied and modified from deap
	# needed to account for bistring only being one of two components of individual
	"""Executes a two-point crossover on the input :term:`sequence`
    individuals. The two individuals are modified in place and both keep
    their original length.
    :returns: A tuple of two individuals.
    This function uses the :func:`~random.randint` function from the Python
    base :mod:`random` module.

    Modified to cross over between rules
    """
	size = len(ind1[0].nodeList)
	cxpointer1 = randint(1, size)
	cxpointer2 = randint(1, size - 1)
	# make sure pointers are in right order
	if cxpointer2 >= cxpointer1:
		cxpointer2 += 1
	else: # Swap the two cx points
		cxpointer1, cxpointer2 = cxpointer2, cxpointer1
	cxpoint1=ind1[0].individualParse[cxpointer1]
	cxpoint2=ind1[0].individualParse[cxpointer2]
	# cross over both bitlists and the andNodeLists (as well as andNodeInvertLists)
	ind1[1][cxpoint1:cxpoint2], ind2[1][cxpoint1:cxpoint2] = ind2[1][cxpoint1:cxpoint2], ind1[1][cxpoint1:cxpoint2]
	ind1[0].andNodeList[cxpointer1:cxpointer2], ind2[0].andNodeList[cxpointer1:cxpointer2] = ind2[0].andNodeList[cxpointer1:cxpointer2], ind1[0].andNodeList[cxpointer1:cxpointer2]
	ind1[0].andNodeInvertList[cxpointer1:cxpointer2], ind2[0].andNodeInvertList[cxpointer1:cxpointer2] = ind2[0].andNodeInvertList[cxpointer1:cxpointer2], ind1[0].andNodeInvertList[cxpointer1:cxpointer2]
	# update the arrays seen by C code updateBool
	ind1[0].updateCpointers()
	ind2[0].updateCpointers()
	return ind1, ind2

# sets up GA toolbox from deap
def buildToolbox( individualLength, bitFlipProb, model, params):

	toolbox = base.Toolbox() # build baseline toolbox
	weightTup=(-1.0,) # specify weights of the errors
	for i in range(len(model.nodeList)-1):
		weightTup+=(-1.0,)
	creator.create("FitnessMin", base.Fitness, weights=weightTup) # make a fitness minimization function
	creator.create("Individual", list, fitness=creator.FitnessMin)	# create a class of individuals that are lists

	#register our bitsring generator and how to create an individual, population
	toolbox.register("genRandomBitString", genBits, model=model)
	toolbox.register("Individual", tools.initIterate, creator.Individual, toolbox.genRandomBitString)
	toolbox.register("population", tools.initRepeat, list , toolbox.Individual)
	#create statistics toolbox and give it functions
	stats = tools.Statistics(key=lambda ind: ind.fitness.values)
	stats.register("avg", numpy.mean)
	stats.register("std", numpy.std)
	stats.register("min", numpy.min)
	stats.register("max", numpy.max)

	# finish registering the toolbox functions
	toolbox.register("mate", tools.cxTwoPoint)
	toolbox.register("mutate", tools.mutFlipBit, indpb=bitFlipProb)
	toolbox.register("select", selNSGA2)
	toolbox.register("similar", numpy.array_equal)
	return toolbox, stats

# calculate  fitness for an individual
def evaluateByNode(individual, cells, model,  sss, params, KOlist, KIlist, boolC):
	boolValues=[NP(list(individual), model, cells, model.initValueList[i], params, KOlist[i], KIlist[i], boolC) for i in range(len(sss))]
	return tuple([numpy.sum([(boolValues[j][i]-sss[j][model.nodeList[i]])**2 for j in range(0,len(sss))]) for i in range(0, len(model.nodeList))])

# generates a random set of samples made up of cells by using parameteris from probInit seq
# to set up then iterating using strict Boolean modeling.
def runProbabilityBooleanSims(individual, model, sampleNum, cells, params, KOlist, KIlist, boolC):
	seeds=[]
	for i in range(0,sampleNum):
		seeds.append(random())
	samples=[sampler(individual, model, sampleNum, seeds[i], params, KOlist[i], KIlist[i], boolC) for i in range(sampleNum)]
	return samples

# generates random seed samples... i.e. generates random starting states then runs EBN
def sampler(individual, model, cells, seeder, params, KOs, KIs, boolC):
	seed(seeder)
	cellArray=[]
	sampleProbs=[]
	# generate random proportions for each node to start
	for j in range(0,len(model.nodeList)):
		sampleProbs.append(random())
	return NP(individual, model, cells, sampleProbs, params, KOs, KIs, boolC)

# generates list of offspring to be compared... decides to do crossover or mutation
def varOrAdaptive(population, toolbox, model, lambda_, cxpb, mutpb, genfrac, mutModel):
	# algorithm for generating a list of offspring... copied and pasted from DEAP with modification for adaptive mutation
	assert (cxpb + mutpb) <= 1.0, ("The sum of the crossover and mutation "
		"probabilities must be smaller or equal to 1.0.")
	offspring = []
	for _ in xrange(lambda_):
		op_choice = random()
		if op_choice < cxpb:            # Apply crossover
			ind1, ind2 = map(toolbox.clone, sample(population, 2))
			ind1, ind2 = cxTwoPointNode(ind1, ind2)
			del ind1.fitness.values
			offspring.append(ind1)
		elif op_choice < cxpb + mutpb:  # Apply mutation
			ind = toolbox.clone(choice(population))
			ind, = mutFlipBitAdapt(ind, genfrac, mutModel)
			del ind.fitness.values
			offspring.append(ind)
		else:                           # shouldn't happen... clone existing individual
			offspring.append(choice(population))
	return offspring

# select node to mutate
def selectMutNode(errors):
	normerrors=[1.*error/numpy.sum(errors) for error in errors]# normalize errors to get a probability that the node  is modified
	probs=numpy.cumsum(normerrors)
	randy=random()# randomly select a node to mutate
	return next(i for i in range(len(probs)) if probs[i]>randy)

# mutation algorithm
def mutFlipBitAdapt(indyIn, genfrac, mutModel):
	errors=list(indyIn.fitness.values) # get errors
	individual=indyIn[1]
	model=indyIn[0]
	# get rid of errors in nodes that can't be changed
	errorNodes=0
	for j in xrange(len(errors)):
		if model.andLenList[j]<2:
			errors[j]=0
		else:
			errorNodes=errorNodes+1

	if numpy.sum(errors)<.05*errorNodes or errorNodes==0:
		# condition selection on number of incoming edges + downstream edges
		pseudoerrors=[len(model.possibilityList[i]) if model.successorNums[i]==0 else len(model.possibilityList[i])*model.successorNums[i] for i in range(len(model.nodeList))]
		# zero out nodes that can't be changed
		for j in xrange(len(pseudoerrors)):
			if model.andLenList[j]<2:
				pseudoerrors[j]=0
		focusNode=selectMutNode(pseudoerrors)
	else:
		# if errors are relatively high, focus on nodes that fit the worst and have highest in-degree
		# calculate probabilities for mutating each node
		for i in range(len(errors)):
			temper=model.successorNums[i]
			if temper==0:
				errors[i]=errors[i]*len(model.possibilityList[i])
			else:
				errors[i]=errors[i]*len(model.possibilityList[i])*temper
		focusNode=selectMutNode(errors)
	# perform mutation
	if model.andLenList[focusNode]>1:
		# find ends of the node of interest in the individual
		start=model.individualParse[focusNode]
		end=findEnd(focusNode,model)
		# mutate the inputs some of the time
		if len(model.possibilityList[focusNode])>3 and random()<mutModel:
			temppermup=[]
			upstreamAdders=list(model.possibilityList[focusNode])
			rvals=list(model.rvalues[focusNode])
			while len(temppermup)<3:
				randy=random()# randomly select a node to mutate
				tempsum=sum(rvals)
				if tempsum==0:
					addNoder=int(math.floor(random()*len(upstreamAdders)))
				else:
					recalc=numpy.cumsum([1.*rval/tempsum for rval in rvals])
					addNoder=next(i for i in range(len(recalc)) if recalc[i]>randy)
				temppermup.append(upstreamAdders.pop(addNoder))
				rvals.pop(addNoder)
			model.update_upstream(focusNode,temppermup)
			model.updateCpointers()
		for i in range(start,end):
			if random()< 2/(end-start+1):
				individual[i] = 1
			else:
				individual[i] = 0
		#ensure that there is at least one shadow and node turned on
		if numpy.sum(individual[start:end])==0:
			individual[start]=1
		indyIn[0]=model
		indyIn[1]=individual
	else:
		print('did not actually check')
	return indyIn,
def selNSGA2(individuals, k):
	# NSGA2 selection taken from deap
	"""Apply NSGA-II selection operator on the *individuals*. Usually, the
	size of *individuals* will be larger than *k* because any individual
	present in *individuals* will appear in the returned list at most once.
	Having the size of *individuals* equals to *k* will have no effect other
	than sorting the population according to their front rank. The
	list returned contains references to the input *individuals*. For more
	details on the NSGA-II operator see [Deb2002]_.

	:param individuals: A list of individuals to select from.
	:param k: The number of individuals to select.
	:returns: A list of selected individuals.

	.. [Deb2002] Deb, Pratab, Agarwal, and Meyarivan, "A fast elitist
	   non-dominated sorting genetic algorithm for multi-objective
	   optimization: NSGA-II", 2002.
	"""
	pareto_fronts = sortNondominatedAdapt(individuals, k)
	for front in pareto_fronts:
		assignCrowdingDist(front)

	chosen = list(chain(*pareto_fronts[:-1]))
	k = k - len(chosen)
	if k > 0:
		sorted_front = sorted(pareto_fronts[-1], key=attrgetter("fitness.crowding_dist"), reverse=True)
		chosen.extend(sorted_front[:k])

	return chosen

# taken from deap and modified slightly to make pareto sorting less strict
def sortNondominatedAdapt(individuals, k, first_front_only=False):
	"""Sort the first *k* *individuals* into different nondomination levels
	using the "Fast Nondominated Sorting Approach" proposed by Deb et al.,
	see [Deb2002]_. This algorithm has a time complexity of :math:`O(MN^2)`,
	where :math:`M` is the number of objectives and :math:`N` the number of
	individuals.

	:param individuals: A list of individuals to select from.
	:param k: The number of individuals to select.
	:param first_front_only: If :obj:`True` sort only the first front and
							 exit.
	:returns: A list of Pareto fronts (lists), the first list includes
			  nondominated individuals.
	.. [Deb2002] Deb, Pratab, Agarwal, and Meyarivan, "A fast elitist
	   non-dominated sorting genetic algorithm for multi-objective
	   optimization: NSGA-II", 2002.
	"""
	if k == 0:
		return []

	map_fit_ind = defaultdict(list)
	for ind in individuals:
		map_fit_ind[ind.fitness].append(ind)
	fits = map_fit_ind.keys()

	current_front = []
	next_front = []
	dominating_fits = defaultdict(int)
	dominated_fits = defaultdict(list)

	# Rank first Pareto front
	for i, fit_i in enumerate(fits):
		for fit_j in fits[i+1:]:
			if dominated(fit_i, fit_j):
				dominating_fits[fit_j] += 1
				dominated_fits[fit_i].append(fit_j)
			elif dominated(fit_j, fit_i):
				dominating_fits[fit_i] += 1
				dominated_fits[fit_j].append(fit_i)
		if dominating_fits[fit_i] == 0:
			current_front.append(fit_i)

	fronts = [[]]
	for fit in current_front:
		fronts[-1].extend(map_fit_ind[fit])
	pareto_sorted = len(fronts[-1])

	# Rank the next front until all individuals are sorted or
	# the given number of individual are sorted.
	if not first_front_only:
		N = min(len(individuals), k)
		while pareto_sorted < N:
			fronts.append([])
			for fit_p in current_front:
				for fit_d in dominated_fits[fit_p]:
					dominating_fits[fit_d] -= 1
					if dominating_fits[fit_d] == 0:
						next_front.append(fit_d)
						pareto_sorted += len(map_fit_ind[fit_d])
						fronts[-1].extend(map_fit_ind[fit_d])
			current_front = next_front
			next_front = []

	return fronts
# taken from deap and modified slightly to make pareto sorting less strict
def dominated(ind1, ind2):
	"""Return true if each objective of *self* is not strictly worse than
		the corresponding objective of *other* and at least one objective is
		strictly better.
		:param obj: Slice indicating on which objectives the domination is
					tested. The default value is `slice(None)`, representing
					every objectives.
	"""
	not_equal = False
	mean1=numpy.mean(ind1.wvalues)
	mean2=numpy.mean(ind2.wvalues)
	std1=numpy.std(ind1.wvalues)
	if mean1 > mean2 :
		not_equal = True
	elif mean1 < mean2:
		return False
	return not_equal
# taken from deap
def assignCrowdingDist(individuals):
	"""Assign a crowding distance to each individual's fitness. The
	crowding distance can be retrieve via the :attr:`crowding_dist`
	attribute of each individual's fitness.
	"""
	if len(individuals) == 0:
		return

	distances = [0.0] * len(individuals)
	crowd = [(ind.fitness.values, i) for i, ind in enumerate(individuals)]

	nobj = len(individuals[0].fitness.values)

	for i in xrange(nobj):
		crowd.sort(key=lambda element: element[0][i])
		distances[crowd[0][1]] = float("inf")
		distances[crowd[-1][1]] = float("inf")
		if crowd[-1][0][i] == crowd[0][0][i]:
			continue
		norm = nobj * float(crowd[-1][0][i] - crowd[0][0][i])
		for prev, cur, next in zip(crowd[:-2], crowd[1:-1], crowd[2:]):
			distances[cur[1]] += 1.*(next[0][i] - prev[0][i]) / norm

	for i, dist in enumerate(distances):
		individuals[i].fitness.crowding_dist = dist

# master GA algorithm
def eaMuPlusLambdaAdaptive( toolbox, model, mu, lambda_, cxpb, mutpb, ngen, namer, newSSS,KOlist, KIlist, params ,boolC,stats=None, verbose=__debug__):
	# population=[[copy.deepcopy(model),genBits(model)]for i in range(params.popSize)]
	population=toolbox.population(n=params.popSize)
	mutModel=params.mutModel
	logbook = tools.Logbook()
	logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])
	lastcheck=[]
	modellist=[]
	fitnesslist=[]
	popList=[]
	# Evaluate the individuals with an invalid fitness
	invalid_ind = [ind for ind in population if not ind.fitness.valid]
	fitnesses=[evaluateByNode(indy[1], params.cells, indy[0],  newSSS, params, KOlist, KIlist, boolC) for indy in invalid_ind]

	for ind, fit in zip(invalid_ind, fitnesses):
		ind.fitness.values = fit
	fitnesslist.append([list(ind.fitness.values) for ind in population])
	popList.append([list(inder[1]) for inder in population])
	modellist.append([[(modeler[0].size), list(modeler[0].nodeList), list(modeler[0].individualParse), list(modeler[0].andNodeList) , list(modeler[0].andNodeInvertList), list(modeler[0].andLenList),	list(modeler[0].nodeList), dict(modeler[0].nodeDict), list(modeler[0].initValueList)] for modeler in population])

	record = stats.compile(population) if stats is not None else {}
	logbook.record(gen=0, nevals=len(invalid_ind), **record)
	if verbose:
		print(logbook.stream)

	breaker=True
	for j in range(len(model.andLenList)):
		if model.andLenList[j]>1:
			breaker=False
	for ind in population:
		if numpy.sum(ind.fitness.values)< .01*len(ind.fitness.values):
			breaker=True
	if breaker:
		# outputList=[fitnesslist, popList, modellist]
		# pickle.dump( outputList, open( namer+"_pops.pickle", "wb" ) )
		return population, logbook

	# Begin the generational process
	for gen in range(1, ngen+1):
		offspring = varOrAdaptive(population, toolbox, model, lambda_, .5+.5*(1.-1.*gen/ngen), (.5*gen/ngen), (1.*gen/ngen),mutModel)
		# Evaluate the individuals with an invalid fitness
		invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
		fitnesses=[evaluateByNode(indy[1], params.cells, indy[0],  newSSS, params, KOlist, KIlist, boolC) for indy in invalid_ind]
		for ind, fit in zip(invalid_ind, fitnesses):
			ind.fitness.values = fit
		# Select the next generation population
		population[:] = toolbox.select(population + offspring, mu)
		fitnesslist.append([list(ind.fitness.values) for ind in population])
		popList.append([list(inder[1]) for inder in population])
		modellist.append([[(modeler[0].size), list(modeler[0].nodeList), list(modeler[0].individualParse), list(modeler[0].andNodeList) , list(modeler[0].andNodeInvertList), list(modeler[0].andLenList),	list(modeler[0].nodeList), dict(modeler[0].nodeDict), list(modeler[0].initValueList)] for modeler in population])

		# Update the statistics with the new population
		record = stats.compile(population) if stats is not None else {}
		logbook.record(gen=gen, nevals=len(invalid_ind), **record)
		if verbose:
			print(logbook.stream)
		breaker=False
		for ind in population:
			if numpy.sum(ind.fitness.values)< .01*len(ind.fitness.values):
				breaker=True
				saveInd=ind
		if breaker:
			errorTemp=saveInd.fitness.values
			for value in errorTemp:
				if value> .1:
					breaker=False
		if breaker:
			# outputList=[fitnesslist, popList, modellist]
			# pickle.dump( outputList, open( namer+"_pops.pickle", "wb" ) )
			return population, logbook

	# outputList=[fitnesslist, popList, modellist]
	# pickle.dump( outputList, open( namer+"_pops.pickle", "wb" ) )
	return population, logbook

# wrapper for GA. eaMuPlusLambdaAdaptive does heavy lifting
def GAsearchModel(model, sampleList,params, KOlist, KIlist, namer, boolC):
	newInitValueList= genInitValueList(sampleList,model) # set up initial value list
	model.initValueList=newInitValueList # append initial value list to model
	toolbox, stats=buildToolbox(model.size,params.bitFlipProb, model, params) # set up toolbox
	# run GA, find best in population, return
	population, logbook=eaMuPlusLambdaAdaptive(toolbox, model, mu=params.mu, lambda_=params.lambd, stats=stats, cxpb=params.crossoverProb, mutpb=params.mutationProb, ngen=params.generations, namer=namer, newSSS= sampleList,KOlist=KOlist, KIlist=KIlist, params=params,  verbose=params.verbose, boolC=boolC)
	out1, out2, model  = findPopBest(population)
	return model,out1,out2

# wrapper to do local search in parrallell manner
def localSearch(model, indy, newSSS, params, KOlist, KIlist, boolC):
	outputs=[checkNodePossibilities(node, indy, newSSS, params.cells, model,params, KOlist, KIlist , boolC) for node in range(len(model.nodeList))]
	equivs=[]
	individual=[]
	devs=[]
	for output in outputs:
		individual.extend(output[0])
		equivs.append(output[1])
		devs.append(output[2])
	return individual, equivs, devs


# local search function
def checkNodePossibilities(node, indy, newSSS, cellNum, model,params, KOlist, KIlist, boolC ):
	tol=.01*len(newSSS) # set tolerance for equivalence
	end=findEnd(node,model) # find end of model for this node
	start=model.individualParse[node] # find start of model for this node
	truth=list(indy[start:end])
	equivs=[truth] 	#add if
	if (end-start)==0:
		return truth, equivs, 0.
	indOptions=[]
	indErrors=[]
	# iterate over possibilities for this node
	for i in range(1,2**(end-start)):
		# printfile=open('output_file.txt','w')
		# printfile.write(str(i)+'\n')
		tempultimate=list(indy)
		tempInd=bitList(i, len(truth))
		# printfile.write(str(tempInd)+'\n')
		# printfile.close()
		tempultimate[start:end]=tempInd # set rule to one being checked
		currentsumtemp=evaluateByNode(tempultimate, cellNum, model,  newSSS, params,KOlist, KIlist , boolC)
		currentsum=currentsumtemp[node] # save the error found
		indOptions.append(tempInd)
		indErrors.append(currentsum)
		gc.collect()
	minny= min(indErrors)
	equivs=[]
	# find the minimum error individual
	for i in range(len(indOptions)):
		if indErrors[i]< minny+tol:
			equivs.append(indOptions[i])
	truth=equivs[0]
	return truth, equivs, minny