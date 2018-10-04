#include <stdio.h>
#include <math.h>
#include <string.h>


int updateBool(int currentNode, int *oldValue,int *nodeIndividual, int andNodes[7][3], int andNodeInvertList[7][3], int nodeStart, int nodeEnd)
{
// we update node by updating shadow and nodes then combining them to update or nodes.
	//update nodes with more than one input
	// first deal with case of simple logic without need of linear regression
	int counter =0;
	int orset[300];
	int andindex;
	int indindex;
	int addnode;
	int newval;

	// go through list of possible shadow and nodes to see which ones actually contribute
	for( indindex=nodeStart; indindex< nodeEnd; indindex++)
	{
		andindex=indindex-nodeStart;
		if ( nodeIndividual[indindex])
		{
			// if a shadow and contributes, compute its value using its upstream nodes
			// calculate value of first then use and to append rest in list of predecessors
			newval=(oldValue[andNodes[andindex][0]]!=andNodeInvertList[andindex][0]);
			// printf("line 45\n");

			for( addnode=1; addnode < 3; addnode++)
			{
				// printf("line 49\n");
				if(andNodes[andindex][addnode]>(-1)){newval=(newval && (oldValue[andNodes[andindex][addnode]]!=andNodeInvertList[andindex][addnode]));}
			}
			// printf("line 53\n");
			orset[counter]=newval;
			counter++;
			
		}
	}
	//combine the shadow and nodes with or operations
	newval=orset[0];
	int q;
	for(q = 1; q<counter; q++)
	{
		newval= (newval || orset[q]);
	}
	return newval;
}

void syncBool(int *avg,int *individual,int indLen, int nodeNum, int *andLenList, int *individualParse, int andNodes[500][7][3], int andNodeInvertList[500][7][3], int simSteps, int initValues[500], int *knockouts, int *knockins){
	// do simulation. individual specifies the particular logic rules on the model. params is a generic holder for simulation parameters. 
	// set up data storage for simulation, add step 0
	int simData[500][1000];
	//iterate over number of steps necessary
	int step;
	int i;
	int nodeEnd;
	int temp;
	int nodeStart;
	int oldValue[500];
	int newValue[500];

	for(i=0; i<nodeNum; i++){newValue[i]=initValues[i];}

	for( step=0; step < simSteps; step++){
		for(i=0; i<nodeNum; i++){oldValue[i]=newValue[i];}
		for(i=0; i<nodeNum; i++){
			//find start and finish for each node to update from the individualParse list
			if(knockouts[i]){temp=0;}
			else if(knockouts[i]){temp=1;}
			else if (andLenList[i]==1){temp= (oldValue[i]!=andNodeInvertList[i][0][0]);}
			else if (andLenList[i]==0){temp=oldValue[i];}
			else{
				if (i==(nodeNum-1)){nodeEnd= indLen;}
				else{nodeEnd=individualParse[i+1];}
				nodeStart=individualParse[i];
				temp=updateBool(i, oldValue, individual, andNodes[i], andNodeInvertList[i], nodeStart, nodeEnd);
			}
			newValue[i]=temp;
			simData[i][step]=temp;
		}
	}
	for(i=0;i<nodeNum;i++){avg[i]=newValue[i]; for(step=simSteps-10; step<simSteps-1; step++){avg[i]+=simData[i][step];}}
}
