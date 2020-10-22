#ifndef ___NEEDLEMAN_WUNSCH_H
#define ___NEEDLEMAN_WUNSCH_H

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
#include <climits>
#include "MSA.h"
#include "sim_matrices.h"
using namespace std;

class needleman_wunsch
{
public:

	needleman_wunsch(const vector<string> & seqArray, int similarity_mode, int match_score = 2, int mismatch_score = -1, int gap_open = 5, int gap_extend = 1): 
	  _seqs(seqArray), _similarity_mode(similarity_mode), _match_score(match_score), _mismatch_score(mismatch_score), _gap_open(gap_open), _gap_extend(gap_extend)
	{
		_simMat = sim_matrices(_similarity_mode, _match_score, _mismatch_score);
		_numberOfSequences = _seqs.size();
	} //constructor
	needleman_wunsch(){};
	needleman_wunsch(const needleman_wunsch& other);   //copy constructor

	int getNumberOfSequences() const {return _numberOfSequences;} 
	~needleman_wunsch();
	needleman_wunsch& operator=(const needleman_wunsch& );
	
	void computeAllPairsMSAs(vector<MSA> & pairwiseMsasToFill);
	MSA computeAffineNWForPair(const string A, const string B); // if gep_open = gap_extend --> same as NW without affine gap penalty
	
private:
	vector<string> _seqs;
	sim_matrices _simMat;
	int _numberOfSequences;
	int _similarity_mode;
	int _match_score;
	int _mismatch_score;
	int _gap_open;
	int _gap_extend;
};

#endif