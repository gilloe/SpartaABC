#ifndef ___MSA_DECOMPOSER_H
#define ___MSA_DECOMPOSER_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>
#include "MSA.h"
using namespace std;

class MSA_decomposer
{
public:

	MSA_decomposer(const MSA & inputMSA): _inputMSA(inputMSA)
	{
	} //constructor
	MSA_decomposer(){};
	MSA_decomposer(const MSA_decomposer& other);   //copy constructor

	MSA_decomposer& operator=(const MSA_decomposer& );
	
	void computeAllPairsMSAs(vector<MSA> & pairwiseMsasToFill);
	MSA getPairwise(const string A, const string B);
	
private:
	MSA _inputMSA;
};

#endif