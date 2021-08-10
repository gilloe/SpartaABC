#include <iostream>

using namespace std;

#include "tree.h"

int main2(){
	tree t1("test.txt");

	cout << t1.stringTreeInPhylipTreeFormat(false) << endl;

	vector<tree::nodeP> nodeVec; // this will include the nodes of the entire tree
	t1.getAllNodes(nodeVec,t1.getRoot());
	MDOUBLE sumOfDistances = 0.0;
	for (size_t i=0 ; i<nodeVec.size(); ++i) {
		if (nodeVec[i]->isRoot()) continue;
		sumOfDistances +=nodeVec[i]->dis2father();
	}
	cout << "sum of distances:" << endl;

	cout << sumOfDistances << endl;
	// find the sum of branch lenghs of the tree✅
	size_t minNodeIndex = NULL;
	MDOUBLE minLength = SIZE_MAX;

	for (size_t i=0 ; i<nodeVec.size(); ++i) {
		if (nodeVec[i]->isRoot()) continue;
		if (nodeVec[i]->dis2father() < minLength ){
			minLength = nodeVec[i]->dis2father();
			minNodeIndex = i;
		}
	}
	cout << "minimal branch length:" << endl;
	cout << nodeVec[minNodeIndex]->dis2father() << endl;
	nodeVec[minNodeIndex]->setDisToFather(0.999);
	cout << "updated tree:" << endl;
	cout << t1.stringTreeInPhylipTreeFormat(false) << endl;
	// change the minimal branch lenth to 0.999 and print the new tree string to the screen✅


	//repeat exs 1 and 2, when the input is bifurcating tree✅
	t1.rootAt(nodeVec[minNodeIndex]->father());
	cout << "rerooted tree:" << endl;
	cout << t1.stringTreeInPhylipTreeFormat(false) << endl;

	// try to reroot the tree in a different node

	//remove the leaf "species1"

	// 
	return 0;
}
