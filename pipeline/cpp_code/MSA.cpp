#include "MSA.h"

MSA::MSA(const vector<string> & seqArray) : _alignedSeqs(seqArray)
{
	_numberOfSequences = _alignedSeqs.size();
	//unalignSeqs();
	trimMSAFromAllIndelPositionAndgetSummaryStatisticsFromIndelCounter(); // also trims all indel positions
	setValuesOfIndelSummStats();
	setLongestAndShortestSequenceLengths();
	
	
} //constructor from vector string

MSA::MSA(string filename)
{
	_alignedSeqs = read_fasta_from_file(filename);
	_numberOfSequences = _alignedSeqs.size();
	trimMSAFromAllIndelPositionAndgetSummaryStatisticsFromIndelCounter();
	setValuesOfIndelSummStats();
	setLongestAndShortestSequenceLengths();
	
} //constructor from sequence file in fasta

void MSA::fillUniqueGapsMap()
{
	for(int j=0; j<_numberOfSequences; j++)	
	{
		int previousIsIndel = 0; // the previous character was not indel

		int currStartIndelPoint = -1;
		int currEndIndelPoint = -1;

		
		for(int i=0; i<getMSAlength(); i++) 
		{	
			if ((_alignedSeqs[j][i]=='-') && (previousIsIndel==0)) 
			{
				previousIsIndel = 1;
				currStartIndelPoint = i;
				currEndIndelPoint = i;

			}
			else if ((_alignedSeqs[j][i]=='-') && (previousIsIndel==1))
			{
				currEndIndelPoint++;
			}
			else // we are on a character
			{
				if(currStartIndelPoint == -1)
				{
					previousIsIndel = 0;
					continue; //this is not an indel - just a character
				}

				//we have an indel - put in map
				pair<int,int> curr_pair(currStartIndelPoint,currEndIndelPoint);
				int currLength = currEndIndelPoint - currStartIndelPoint + 1;
				if(_uniqueIndelMap.find(curr_pair) == _uniqueIndelMap.end())
				{
					// new entry
					vector<int> curr_value;
					curr_value.push_back(currLength);
					curr_value.push_back(1); // first gap of this kind
					_uniqueIndelMap[curr_pair] = curr_value;
				}
				else
				{
					// already exists - raise counter
					vector<int> curr_value = _uniqueIndelMap[curr_pair];
					curr_value[1]++;
					_uniqueIndelMap[curr_pair] = curr_value;
				}

				previousIsIndel = 0;
				currStartIndelPoint = -1;
				currEndIndelPoint = -1;
			}
		}
		if(currStartIndelPoint != -1)
		{
			//we have an indel - put in map
			pair<int,int> curr_pair(currStartIndelPoint,currEndIndelPoint);
			int currLength = currEndIndelPoint - currStartIndelPoint + 1;
			if(_uniqueIndelMap.find(curr_pair) == _uniqueIndelMap.end())
			{
				// new entry
				vector<int> curr_value;
				curr_value.push_back(currLength);
				curr_value.push_back(1); // first gap of this kind
				_uniqueIndelMap[curr_pair] = curr_value;
			}
			else
			{
				// already exists - raise counter
				vector<int> curr_value = _uniqueIndelMap[curr_pair];
				curr_value[1]++;
				_uniqueIndelMap[curr_pair] = curr_value;
			}
		}
	}
}

vector<string> MSA::getUnalignedSeqs() const {
	vector<string> unalignedSeqs;
	for(int j=0; j<_numberOfSequences; j++)	
	{
		string curr_unaligned = "";
		for(int i=0; i<getMSAlength(); i++) 
		{
			if(_alignedSeqs[j][i]!='-')
			{
				curr_unaligned += _alignedSeqs[j][i];
			}
		}
		unalignedSeqs.push_back(curr_unaligned);
	}
	return unalignedSeqs;
}


void MSA::setValuesOfIndelSummStats()
{	
	_totalNumberOfIndels = 0;
	_totalNumberOfUniqueIndels = 0;
	_numberOfIndelsOfLengthOne = 0;
	_numberOfIndelsOfLengthOneInOnePosition = 0;
	_numberOfIndelsOfLengthOneInTwoPositions = 0;
	_numberOfIndelsOfLengthOneInNMinus1Positions = 0;
	_numberOfIndelsOfLengthTwo = 0;
	_numberOfIndelsOfLengthTwoInOnePosition = 0;
	_numberOfIndelsOfLengthTwoInTwoPositions = 0;
	_numberOfIndelsOfLengthTwoInNMinus1Positions = 0;
	_numberOfIndelsOfLengthThree = 0;
	_numberOfIndelsOfLengthThreeInOnePosition = 0;
	_numberOfIndelsOfLengthThreeInTwoPositions = 0;
	_numberOfIndelsOfLengthThreeInNMinus1Positions = 0;
	_numberOfIndelsOfLengthAtLeastFour = 0;
	_numberOfIndelsOfLengthAtLeastFourInOnePosition = 0;
	_numberOfIndelsOfLengthAtLeastFourInTwoPositions = 0;
	_numberOfIndelsOfLengthAtLeastFourInNMinus1Positions = 0;
	_aveIndelLength = 0;
	_aveUniqueIndelLength = 0;

	int total_number_of_gap_chars = 0;
	int total_number_of_unique_gap_chars = 0;
	
	fillUniqueGapsMap();
	// go over the unique gap map
	typedef map<pair<int,int>, vector<int>>::iterator it_type;
	for(it_type iterator = _uniqueIndelMap.begin(); iterator != _uniqueIndelMap.end(); iterator++) 
	{
		pair<int,int> gap_indices = iterator->first;
		vector<int> length_and_count = iterator->second;

		total_number_of_gap_chars += (length_and_count[0] * length_and_count[1]);
		_totalNumberOfIndels += length_and_count[1];
		
		_totalNumberOfUniqueIndels++;
		total_number_of_unique_gap_chars += length_and_count[0];;

		if(length_and_count[0] == 1)
		{
			_numberOfIndelsOfLengthOne += length_and_count[1];
			if (length_and_count[1] == 1) _numberOfIndelsOfLengthOneInOnePosition++;
			if (length_and_count[1] == 2) _numberOfIndelsOfLengthOneInTwoPositions++;
			if (length_and_count[1] == getNumberOfSequences() - 1) _numberOfIndelsOfLengthOneInNMinus1Positions++;
			
		}
		if(length_and_count[0] == 2)
		{
			_numberOfIndelsOfLengthTwo += length_and_count[1];
			if (length_and_count[1] == 1) _numberOfIndelsOfLengthTwoInOnePosition++;
			if (length_and_count[1] == 2) _numberOfIndelsOfLengthTwoInTwoPositions++;
			if (length_and_count[1] == getNumberOfSequences() - 1) _numberOfIndelsOfLengthTwoInNMinus1Positions++;
		}
		if(length_and_count[0] == 3)
		{
			_numberOfIndelsOfLengthThree += length_and_count[1];
			if (length_and_count[1] == 1) _numberOfIndelsOfLengthThreeInOnePosition++;
			if (length_and_count[1] == 2) _numberOfIndelsOfLengthThreeInTwoPositions++;
			if (length_and_count[1] == getNumberOfSequences() - 1) _numberOfIndelsOfLengthThreeInNMinus1Positions++;
		}
		if(length_and_count[0] > 3)
		{
			_numberOfIndelsOfLengthAtLeastFour += length_and_count[1];
			if (length_and_count[1] == 1) _numberOfIndelsOfLengthAtLeastFourInOnePosition++;
			if (length_and_count[1] == 2) _numberOfIndelsOfLengthAtLeastFourInTwoPositions++;
			if (length_and_count[1] == getNumberOfSequences() - 1) _numberOfIndelsOfLengthAtLeastFourInNMinus1Positions++;
		}
	}
	if(_totalNumberOfIndels != 0)
	{
		_aveIndelLength = static_cast<double>(total_number_of_gap_chars)/static_cast<double>(_totalNumberOfIndels);
		_aveUniqueIndelLength = static_cast<double>(total_number_of_unique_gap_chars)/static_cast<double>(_totalNumberOfUniqueIndels);

	}	
}

MSA::MSA(const MSA& other )//copy constructor
{
	_alignedSeqs = other._alignedSeqs;
	//_unalignedSeqs = other._unalignedSeqs;
	_numberOfSequences = other._numberOfSequences;
	_aveIndelLength = other._aveIndelLength;
	_aveUniqueIndelLength = other._aveUniqueIndelLength;
	_totalNumberOfIndels = other._totalNumberOfIndels;
	_totalNumberOfUniqueIndels = other._totalNumberOfUniqueIndels;
	_numberOfIndelsOfLengthOne = other._numberOfIndelsOfLengthOne;
	_numberOfIndelsOfLengthTwo = other._numberOfIndelsOfLengthTwo;
	_numberOfIndelsOfLengthThree = other._numberOfIndelsOfLengthThree;
	_numberOfIndelsOfLengthAtLeastFour = other._numberOfIndelsOfLengthAtLeastFour;
	_longestSeqLength =other._longestSeqLength;
	_shortestSeqLength=other._shortestSeqLength;

}

MSA&  MSA::operator=(const MSA& other)
{
    _alignedSeqs= other._alignedSeqs;
	//_unalignedSeqs = other._unalignedSeqs;
	_numberOfSequences=other._numberOfSequences;
	_aveIndelLength = other._aveIndelLength;
	_aveUniqueIndelLength = other._aveUniqueIndelLength;
	_totalNumberOfIndels = other._totalNumberOfIndels;
	_totalNumberOfUniqueIndels = other._totalNumberOfUniqueIndels;
	_numberOfIndelsOfLengthOne = other._numberOfIndelsOfLengthOne;
	_numberOfIndelsOfLengthTwo = other._numberOfIndelsOfLengthTwo;
	_numberOfIndelsOfLengthThree = other._numberOfIndelsOfLengthThree;
	_numberOfIndelsOfLengthAtLeastFour = other._numberOfIndelsOfLengthAtLeastFour;
	_longestSeqLength = other._longestSeqLength;
	_shortestSeqLength= other._shortestSeqLength;

	return *this;
}

void MSA::setLongestAndShortestSequenceLengths()
{
	vector<int> sequencesLengths(_numberOfSequences,0);
	_longestSeqLength = 0;
	_shortestSeqLength = 0; //initialization inside loop
	for(int i=0; i<_numberOfSequences; i++)
	{
		for(size_t j=0;j<_alignedSeqs[i].size();j++)
		{
			if(_alignedSeqs[i][j]!='-')
				sequencesLengths[i]++;
		}
		if(i == 0)
		{
			_shortestSeqLength = sequencesLengths[i]; //initialize minimum search
		}
		if(sequencesLengths[i]<_shortestSeqLength)
			_shortestSeqLength=sequencesLengths[i];
		if(sequencesLengths[i]>_longestSeqLength)
			_longestSeqLength=sequencesLengths[i];
	}
}

MSA::~MSA()
{
	_alignedSeqs.clear();
}


double MSA::getStatValByType(stat_type statTypeToGet)
{
	switch(statTypeToGet)
	{
		case AVG_GAP_SIZE:
			return getAverageIndelSize();
		case MSA_LEN:
			return (double)getMSAlength();
		case MSA_MAX_LEN:
			return (double)getMSALongestSeqLength();
		case MSA_MIN_LEN:
			return (double)getMSAshortestSeqLength();
		case TOT_NUM_GAPS:
			return (double)getTotalNumberOfIndels();
		case NUM_GAPS_LEN_ONE:
			return (double)getNumberOfIndelsOfLengthOne();
		case NUM_GAPS_LEN_TWO:
			return (double)getNumberOfIndelsOfLengthTwo();
		case NUM_GAPS_LEN_THREE:
			return (double)getNumberOfIndelsOfLengthThree();
		case NUM_GAPS_LEN_AT_LEAST_FOUR:
			return (double)getNumberOfIndelsOfLengthAtLeastFour();
		case AVG_UNIQUE_GAP_SIZE:
			return getAverageUniqueIndelSize();
		case TOT_NUM_UNIQUE_GAPS:
			return (double)getTotalNumberOfUniqueIndels();
		case NUM_GAPS_LEN_ONE_POS_1_GAPS:
			return (double)getnumberOfIndelsOfLengthOneInOnePosition();
		case NUM_GAPS_LEN_ONE_POS_2_GAPS:
			return (double)getnumberOfIndelsOfLengthOneInTwoPositions();
		case NUM_GAPS_LEN_ONE_POS_N_MINUS_1_GAPS:
			return (double)getnumberOfIndelsOfLengthOneInNMinus1Positions();
		case NUM_GAPS_LEN_TWO_POS_1_GAPS:
			return (double)getnumberOfIndelsOfLengthTwoInOnePosition();
		case NUM_GAPS_LEN_TWO_POS_2_GAPS:
			return (double)getnumberOfIndelsOfLengthTwoInTwoPositions();
		case NUM_GAPS_LEN_TWO_POS_N_MINUS_1_GAPS:
			return (double)getnumberOfIndelsOfLengthTwoInNMinus1Positions();
		case NUM_GAPS_LEN_THREE_POS_1_GAPS:
			return (double)getnumberOfIndelsOfLengthThreeInOnePosition();
		case NUM_GAPS_LEN_THREE_POS_2_GAPS:
			return (double)getnumberOfIndelsOfLengthThreeInTwoPositions();
		case NUM_GAPS_LEN_THREE_POS_N_MINUS_1_GAPS:
			return (double)getnumberOfIndelsOfLengthThreeInNMinus1Positions();
		case NUM_GAPS_LEN_AT_LEAST_FOUR_POS_1_GAPS:
			return (double)getnumberOfIndelsOfLengthAtLeastFourInOnePosition();
		case NUM_GAPS_LEN_AT_LEAST_FOUR_POS_2_GAPS:
			return (double)getnumberOfIndelsOfLengthAtLeastFourInTwoPositions();
		case NUM_GAPS_LEN_AT_LEAST_FOUR_POS_N_MINUS_1_GAPS:
			return (double)getnumberOfIndelsOfLengthAtLeastFourInNMinus1Positions();
		case MSA_POSITION_WITH_0_GAPS:
			return (double)getNumberOfMSA_position_with_0_gaps();
		case MSA_POSITION_WITH_1_GAPS:
			return (double)getNumberOfMSA_position_with_1_gaps();
		case MSA_POSITION_WITH_2_GAPS:
			return (double)getNumberOfMSA_position_with_2_gaps();
		case MSA_POSITION_WITH_N_MINUS_1_GAPS:
			return (double)getNumberOfMSA_position_with_n_minus_1_gaps();
	}
}

void MSA::trimMSAFromAllIndelPositionAndgetSummaryStatisticsFromIndelCounter() {
	_numberOfMSA_position_with_0_gaps = 0;
	_numberOfMSA_position_with_1_gaps = 0;
	_numberOfMSA_position_with_2_gaps = 0;
	_numberOfMSA_position_with_n_minus_1_gaps = 0;

	//cout << "BEFORE:" << endl;
	//printMSA();
	_indelCounter.resize(getMSAlength());
	for (int i = 0; i < _numberOfSequences; i++) {
		for (int j = 0; j < getMSAlength(); j++) {
			if (_alignedSeqs[i][j] == '-') _indelCounter[j]++;
		}
	}
	// now we want to remove from the aligned sequences all positions for which the indel counter == _numberOfSequences
	for (int i = 0; i < _numberOfSequences; i++) {
		string::iterator j = _alignedSeqs[i].begin();
		while (j < _alignedSeqs[i].end()) {
			size_t placeInIndelCounter = j - _alignedSeqs[i].begin();
			if (_indelCounter[placeInIndelCounter] == _numberOfSequences) _alignedSeqs[i].erase(j, j + 1);
			if (j != _alignedSeqs[i].end()) j++; // we have to check it because this is after we erased part of the sequence
			
		}
	}

	for (size_t i = 0; i < _indelCounter.size(); ++i) {
		if (_indelCounter[i] == 0) _numberOfMSA_position_with_0_gaps++;
		else if (_indelCounter[i] == 1) _numberOfMSA_position_with_1_gaps++;
		else if (_indelCounter[i] == 2) _numberOfMSA_position_with_2_gaps++;
		else if (_indelCounter[i] == _numberOfSequences-1) _numberOfMSA_position_with_n_minus_1_gaps++;
	}
	//cout << "AFTER:" << endl;
	//printMSA();
}

void MSA::printMSA() {
	for (int i = 0; i < _numberOfSequences; i++) {
		for (int j = 0; j < getMSAlength(); j++) {
			cout<<_alignedSeqs[i][j];
		}
		cout << endl;
	}
}