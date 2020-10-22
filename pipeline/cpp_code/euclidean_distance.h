#ifndef _EUCLIDEAN_DISTANCE
#define _EUCLIDEAN_DISTANCE

#include <vector>
#include "MSA.h"
#include <math.h>

using namespace std;

class euclidean_distance
{
public:
	void setWeightsVector (vector<double>& summStatsWeights);
	vector<double> getWeightsVector();
	double computeWeightedDist(const vector<double> & realSummStats, const vector<double> & simSummStats);
private:
	vector<double> _summStatsWeights;
	double checkLengthsValidity(const vector<double> & realSummStats, const vector<double> & simSummStats);
};
#endif