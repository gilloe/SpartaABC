#ifndef _MSA
#define _MSA
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <map>
#include "read_seqs.h"
#include "summary_stats_defs.h"
using namespace std;

class MSA
{
public:
	// This construction also computes all the summary statistics.
	MSA(const vector<string> & seqArray); 
	
	// This construction reads a FASTA file into the MSA and computes summary statistics
	MSA(string filename);

	MSA(){};
	MSA(const MSA& other);   //copy constructor

	int getMSAlength() const {return _alignedSeqs[0].size();}
	int getNumberOfSequences() const {return _numberOfSequences;} 
	int getTotalNumberOfIndels() const {return _totalNumberOfIndels;}
	int getTotalNumberOfUniqueIndels() const {return _totalNumberOfUniqueIndels;}
	int getNumberOfIndelsOfLengthOne() const {return _numberOfIndelsOfLengthOne;}
	int getNumberOfIndelsOfLengthTwo() const {return _numberOfIndelsOfLengthTwo;}
	int getNumberOfIndelsOfLengthThree() const {return _numberOfIndelsOfLengthThree;}
	int getNumberOfIndelsOfLengthAtLeastFour() const {return _numberOfIndelsOfLengthAtLeastFour;}
	double getAverageIndelSize() const {return _aveIndelLength;}
	double getAverageUniqueIndelSize() const {return _aveUniqueIndelLength;}
	int getMSALongestSeqLength()const {return _longestSeqLength;}
	int getMSAshortestSeqLength()const {return _shortestSeqLength;}
	int getnumberOfIndelsOfLengthOneInOnePosition()const {return _numberOfIndelsOfLengthOneInOnePosition;}
	int getnumberOfIndelsOfLengthOneInTwoPositions()const {return _numberOfIndelsOfLengthOneInTwoPositions;}
	int getnumberOfIndelsOfLengthOneInNMinus1Positions()const {return _numberOfIndelsOfLengthOneInNMinus1Positions;}
	int getnumberOfIndelsOfLengthTwoInOnePosition()const {return _numberOfIndelsOfLengthTwoInOnePosition;}
	int getnumberOfIndelsOfLengthTwoInTwoPositions()const {return _numberOfIndelsOfLengthTwoInTwoPositions;}
	int getnumberOfIndelsOfLengthTwoInNMinus1Positions()const {return _numberOfIndelsOfLengthTwoInNMinus1Positions;}
	int getnumberOfIndelsOfLengthThreeInOnePosition()const {return _numberOfIndelsOfLengthThreeInOnePosition;}
	int getnumberOfIndelsOfLengthThreeInTwoPositions()const {return _numberOfIndelsOfLengthThreeInTwoPositions;}
	int getnumberOfIndelsOfLengthThreeInNMinus1Positions()const {return _numberOfIndelsOfLengthThreeInNMinus1Positions;}
	int getnumberOfIndelsOfLengthAtLeastFourInOnePosition()const {return _numberOfIndelsOfLengthAtLeastFourInOnePosition;}
	int getnumberOfIndelsOfLengthAtLeastFourInTwoPositions()const {return _numberOfIndelsOfLengthAtLeastFourInTwoPositions;}
	int getnumberOfIndelsOfLengthAtLeastFourInNMinus1Positions()const {return _numberOfIndelsOfLengthAtLeastFourInNMinus1Positions;}
	size_t getNumberOfMSA_position_with_0_gaps() const {return _numberOfMSA_position_with_0_gaps;}
	size_t getNumberOfMSA_position_with_1_gaps() const { return _numberOfMSA_position_with_1_gaps;}
	size_t getNumberOfMSA_position_with_2_gaps() const { return _numberOfMSA_position_with_2_gaps;}
	size_t getNumberOfMSA_position_with_n_minus_1_gaps() const { return _numberOfMSA_position_with_n_minus_1_gaps;}

	double getStatValByType(stat_type statTypeToGet);
	vector<string> getUnalignedSeqs() const;
	vector<string> getAlignedSeqs () const {return _alignedSeqs;}
	void printMSA();
	
	~MSA();
	MSA& operator=(const MSA& );

private:
	int _numberOfSequences; // NUMBER OF SEQUENCES IN THE MSA
	double _aveIndelLength;
	int _totalNumberOfIndels;
	int _longestSeqLength;
	int _shortestSeqLength;
	vector<string> _alignedSeqs; //The aligned sequences
	//vector<string> _unalignedSeqs; //The unaligned sequences

	int _numberOfIndelsOfLengthOne; //counts the number of indels of size 1 over all sequences
	int _numberOfIndelsOfLengthTwo; //counts the number of indels of size 2 over all sequences
	int _numberOfIndelsOfLengthThree; //counts the number of indels of size 3 over all sequences
	int _numberOfIndelsOfLengthAtLeastFour; //counts the number of indels of size 4+ over all sequences
	
	int _numberOfIndelsOfLengthOneInOnePosition;
	int _numberOfIndelsOfLengthOneInTwoPositions;
	int _numberOfIndelsOfLengthOneInNMinus1Positions;
	int _numberOfIndelsOfLengthTwoInOnePosition;
	int _numberOfIndelsOfLengthTwoInTwoPositions;
	int	_numberOfIndelsOfLengthTwoInNMinus1Positions;
	int _numberOfIndelsOfLengthThreeInOnePosition;
	int _numberOfIndelsOfLengthThreeInTwoPositions;
	int	_numberOfIndelsOfLengthThreeInNMinus1Positions;
	int _numberOfIndelsOfLengthAtLeastFourInOnePosition;
	int _numberOfIndelsOfLengthAtLeastFourInTwoPositions;
	int	_numberOfIndelsOfLengthAtLeastFourInNMinus1Positions;

	size_t _numberOfMSA_position_with_0_gaps;
	size_t _numberOfMSA_position_with_1_gaps;
	size_t _numberOfMSA_position_with_2_gaps;
	size_t _numberOfMSA_position_with_n_minus_1_gaps;
	

	map<pair<int,int>,vector<int>> _uniqueIndelMap;
	vector<int> _indelCounter;
	// unique indels summary statistics
	double _aveUniqueIndelLength;
	int _totalNumberOfUniqueIndels;

	//methods
	void setValuesOfIndelSummStats();
	void fillUniqueGapsMap();
	//void unalignSeqs();
	void setLongestAndShortestSequenceLengths();
	void trimMSAFromAllIndelPositionAndgetSummaryStatisticsFromIndelCounter();
};
#endif