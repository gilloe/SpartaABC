
#include "MSA_decomposer.h"


MSA_decomposer::MSA_decomposer(const MSA_decomposer& other )//copy constructor
{
	_inputMSA = other._inputMSA;
}

void MSA_decomposer::computeAllPairsMSAs(vector<MSA> & pairwiseMsasToFill)
{
	int numberOfSequences = _inputMSA.getNumberOfSequences();
	vector<string> alignedSeqs = _inputMSA.getAlignedSeqs();

	for(int i = 0; i < numberOfSequences; i++)
	{
		for(int j = i + 1; j < numberOfSequences; j++)
		{
			MSA curr_pairwise = getPairwise(alignedSeqs[i],alignedSeqs[j]);
			pairwiseMsasToFill.push_back(curr_pairwise);
		}
	}
}

MSA MSA_decomposer::getPairwise(const string A, const string B)
{
	string charsToKeepA = "";
	string charsToKeepB = "";
	
	int originalLength = A.length(); //as A and B come from an MSA, they are of the same length

	for(int i = 0; i < originalLength; i++)
	{
		if((A[i] == '-') && (B[i] == '-'))
		{
			continue;
		}
		charsToKeepA += A[i];
		charsToKeepB += B[i];
	}

	vector<string> aligned_seqs;
	aligned_seqs.push_back(charsToKeepA);
	aligned_seqs.push_back(charsToKeepB);
	MSA curr_pairwise = MSA(aligned_seqs);
    return curr_pairwise;
}

