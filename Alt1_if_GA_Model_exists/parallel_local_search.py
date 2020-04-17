from GA import *
from simulation import *
from utils import *
from pathway_analysis_score_nodes import *
from analysis_accuracy import *
import time as time

if __name__ == '__main__':
    start_time = time.time()
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
    print("Individual: "+str(indy))
    KOlist=GAmodel[2]
    print("Knock-out list: "+str(KOlist))
    KIlist=GAmodel[3]
    print("Knock-in list: "+str(KIlist))

    model = modelHolder(GAmodel[0]) #[ size, nodeList, individualParse, andNodeList, andNodeInvertList, andLenList,nodeList, nodeDict, initValueList]
    model.modelHolder_updateCpointers()

    #Load C libraries
    updateBooler=cdll.LoadLibrary('./simulator.so')
    boolC=updateBooler.syncBool

    #Initiate params class
    params=paramClass()

    #Start local search
    nodeIndex = model.nodeList.index(node)
    truth, equivs, minny = checkNodePossibilities(nodeIndex, indy, newSSS, params.cells, model, params, KOlist, KIlist, boolC)
    print(equivs)
    pickle.dump([truth, equivs, minny], open(str(name)+"_"+str(node)+"_local1.pickle", "wb"))

    end_time = time.time()
    time_taken = end_time - start_time
    print("Time taken: " + str(time_taken))
