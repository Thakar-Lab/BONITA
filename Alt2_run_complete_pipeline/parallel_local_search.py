from GA import *
from simulation import *
from utils import *
from pathway_analysis_score_nodes import *
import time as time


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

    # setup C pointers with correct lengths to pass to simulation software in C
    def modelHolder_updateCpointers(self):
        tempandnoder=[]
        tempandinverter=[]
        for currentNode in range(1000):
            tempAndNodes=[]
            tempandNodeInvertList=[]
            if currentNode<len(self.nodeList):
                tempAndNodes=[xi+[-1]*(3-len(xi)) for xi in self.andNodeList[currentNode]]
                tempandNodeInvertList=[xi+[-1]*(3-len(xi)) for xi in self.andNodeInvertList[currentNode]]
            while(len(tempAndNodes)<7):
                tempAndNodes.append([0,0,0])
                tempandNodeInvertList.append([0,0,0])
            tempandnoder.append(tempAndNodes)
            tempandinverter.append(tempandNodeInvertList)
        self.andNodeInvert=np.array(tempandinverter, dtype=np.intc, order='C')
        self.andNodes=np.array(tempandnoder, dtype=np.intc, order='C')

if __name__ == '__main__':
    start_time = time.time()
    print('loaded\n')
    # read in arguments from shell scripts
    parser = argparse.ArgumentParser()
    parser.add_argument("node") #string
    parser.add_argument("GAmodel") #pickle file
    parser.add_argument("newSSS") #pickle file
    results=parser.parse_args()
    #node=int(results.node)
    node=str(results.node)
    GAmodel=results.GAmodel
    newSSS=results.newSSS
    name= GAmodel[:-15]           #WP23_edgeList_edited.txt_net_2.rules_5_GAmodel.pickle
    #recover sample list
    newSSS=pickle.load(open(newSSS, "rb"))
    #Process GAmodel pickle to generate input for checkNodePossibilities
    GAmodel = pickle.load(open(GAmodel, "rb")) #[storeModel1, bruteOut, knockoutLists, knockinLists]
    indy=GAmodel[1]
    print("\nIndividual: "+str(indy))
    KOlist=GAmodel[2]
    print("\nKnock-out list: "+str(KOlist))
    KIlist=GAmodel[3]
    print("\nKnock-in list: "+str(KIlist))
    model = modelHolder(GAmodel[0]) #[ size, nodeList, individualParse, andNodeList, andNodeInvertList, andLenList,nodeList, nodeDict, initValueList]
    model.modelHolder_updateCpointers()
    #Load C libraries
    updateBooler=cdll.LoadLibrary('./simulator.so')
    boolC=updateBooler.syncBool
    #Initiate params class
    params=paramClass()
    #Start local search
    nodeIndex = model.nodeList.index(node)

    cellArray=[]
    simSteps= 100
    # set up knockin and knockout lists
    knockins=np.zeros(len(model.nodeList),dtype=np.intc, order='C')
    knockouts=np.zeros(len(model.nodeList),dtype=np.intc, order='C')
    for knocker in KOlist:
        knockouts[knocker]=1
    for knocker in KIlist:
        knockins[knocker]=1
    # put objects in correct format for passing to C
    nodeIndividual=np.array(indy, dtype=np.intc, order='C')
    indLen=len(indy)
    nodeNum=len(model.nodeList)
    individualParse=np.array(model.individualParse, dtype=np.intc, order='C')

    andLenList=np.array(model.andLenList, dtype=np.intc, order='C')

    # convert objects into C pointers
    print(model.nodeList)

    truth, equivs, minny = checkNodePossibilities(nodeIndex, indy, newSSS, params.cells, model, params, KOlist, KIlist, boolC)
    pickle.dump([truth, equivs, minny], open(str(name)+"_"+str(node)+"_local1.pickle", "wb"))

    end_time = time.time()
    time_taken = end_time - start_time
    print("Time taken: " + str(time_taken))
