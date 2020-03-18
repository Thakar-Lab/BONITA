from GA import *
from simulation import *
from utils import *
from pathway_analysis_score_nodes import *

if __name__ == '__main__':

    # read in arguments from shell scripts
    parser = argparse.ArgumentParser()
    parser.add_argument("node")

    #these should all be pickle file names
    parser.add_argument("model")
    parser.add_argument("indy")
    parser.add_argument("newSSS")
    parser.add_argument("params")
    parser.add_argument("KOlist")
    parser.add_argument("KIlist")

    results = parser.parse_args()
    node=results.node

    #these should all be pickles - TODO - add code in pathway_analysis_score_nodes to dump them out

    model=results.model
    indy=results.indy
    newSSS=results.newSSS
    params=results.params
    KOlist=results.KOlist
    KIlist=results.KIlist

    updateBooler=cdll.LoadLibrary('./simulator.so')
	boolC=updateBooler.syncBool

    model=pickle.load(open(model, "rb"))
    indy=pickle.load(open(indy, "rb"))
    newSSS=pickle.load(open(newSSS, "rb"))
    params=pickle.load(open(params, "rb"))
    KOlist=pickle.load(open(KOlist, "rb"))
    KIlist=pickle.load(open(KIlist, "rb"))

    truth, equivs, minny = checkNodePossibilities(node, indy, newSSS, params.cells, model, params, KOlist, KIlist, boolC)
    pickle.dump([truth, equivs, minny], )