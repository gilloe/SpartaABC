#ifndef _SUMM_STATS_WRAPPER
#define _SUMM_STATS_WRAPPER
#include <vector>
#include "MSA.h"
#include "needleman_wunsch.h"
#include "MSA_decomposer.h"
#include "SpartaABC_options.h"
#include "read_seqs.h"

using namespace std;


vector<double> getStatVec(MSA &currMSA);// This function gets an MSA as input and return a vector of summary statistics.
vector<double> getStatVec(string inputFile); // this function return the summary statistics vector from an input MSA read from a file
string NiceHeader(); // prints a nice header in the output file
vector<double> getWeightsVector(); // gets a weight vector either from an input file or using simulations.

//}

// ------------------------------- OLD VERSION
/*
class summ_stats_wrapper {
public:
	summ_stats_wrapper(	const MSA & curr_msa,
						int mode,
						int similarity_mode) { // constructor from MSA
	
		if(mode == 1) //pairwise with NW
		{
			needleman_wunsch nw_obj(curr_msa.getUnalignedSeqs(),similarity_mode);	
			nw_obj.computeAllPairsMSAs(_msas);
		}
		else if(mode == 2) //pairwise derived from MSA
		{
			MSA_decomposer msa_decomp_obj(curr_msa);
			msa_decomp_obj.computeAllPairsMSAs(_msas);
		}
		else //multiple
		{
			_msas.push_back(curr_msa);
		}
		computeAvgOnEachMsa();
		fillWeightsVec();
	}

	summ_stats_wrapper(	string inputFile,
						int mode,
						int similarity_mode) //constructor from filename - this is called only on input data (not simulated)
	{
		if((mode == 1) || (mode == 2)) //pairwise with NW - from file there is no menaing to taking from MSA
		{
			vector<string> unalignedSeqsFromFile = read_fasta_from_file(inputFile);
			needleman_wunsch nw_obj(unalignedSeqsFromFile,similarity_mode);
			nw_obj.computeAllPairsMSAs(_msas);
		}
		else //multiple, mode == 0
		{
			MSA curr_MSA(inputFile);
			_msas.push_back(curr_MSA);
		}
		computeAvgOnEachMsa();
		fillWeightsVec();
	}

	summ_stats_wrapper(){};
	summ_stats_wrapper(const summ_stats_wrapper& other);   //copy constructor
	summ_stats_wrapper& operator=(const MSA& );
	
	vector<double> getSummStatsValsVector () const {return _avgSummStatVals;}
	vector<double> getSummStatsWeightsVector () const {return _summStatWeights;}
	string getNiceHeader ();

private:
	vector<MSA> _msas;
	vector<double> _avgSummStatVals;
	vector<double> _summStatWeights;
	
	void computeAvgOnEachMsa();
	void fillWeightsVec();
	vector<double> getStatVec (MSA &currMSA);
};
*/
#endif