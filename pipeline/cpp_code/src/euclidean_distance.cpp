#include "euclidean_distance.h"


double euclidean_distance::computeWeightedDist(const vector<double> & realSummStats, const vector<double> & simSummStats)
{
	double validityVal = checkLengthsValidity(realSummStats,simSummStats);
	if(validityVal != 0)
	{
		exit(1);
	}

	double euclideanWeightedDist = 0;
	for(int i = 0; i < realSummStats.size(); i++)
	{
		//cout << "summ statistics " << i << " real value: " << realSummStats[i] << "simulated value: " << simSummStats[i] << " weight: " << _summStatsWeights[i] << "total contribution: " << _summStatsWeights[i] * _summStatsWeights[i] * (realSummStats[i] - simSummStats[i])*(realSummStats[i] - simSummStats[i]) << endl;
		euclideanWeightedDist += _summStatsWeights[i]* _summStatsWeights[i]*(realSummStats[i]-simSummStats[i])*(realSummStats[i] - simSummStats[i]); // square
	}
	euclideanWeightedDist = pow(euclideanWeightedDist,0.5);
	return euclideanWeightedDist;
}

double euclidean_distance::checkLengthsValidity(const vector<double> & realSummStats, const vector<double> & simSummStats)
{
	if(realSummStats.size() != simSummStats.size())
	{
		cout<<"realSummStats length doesn't match simSummStats length"<<endl;
		return -1;
	}
	if(realSummStats.size() != _summStatsWeights.size())
	{
		cout<<"realSummStats length doesn't match summStatsWeights length"<<endl;
		return -1;
	}
	return 0;
}

void euclidean_distance::setWeightsVector(vector<double> & summStatsWeights)
{
	_summStatsWeights = summStatsWeights;
}

vector<double> euclidean_distance::getWeightsVector()
{
	return _summStatsWeights;
}