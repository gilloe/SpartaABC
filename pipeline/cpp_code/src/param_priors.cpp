#include "param_priors.h"

double getRandDoubleParamVal(string distributionName, double minVal, double maxVal)
{
	if(distributionName == "uniform")
	{
		double rand_unif_between_0_and_1 = ((double)rand() / (double)RAND_MAX);
		double rand_unif_in_range = minVal + rand_unif_between_0_and_1 * (maxVal - minVal);
		double rand_unif_in_range_rounded = roundParamVal(rand_unif_in_range,9);
		return rand_unif_in_range_rounded;
	}
	else
	{
		//at the moment - no other implementation
		return -1.0;
	}
}

int getRandIntParamVal(string distributionName, int minVal, int maxVal)
{
	if(distributionName == "uniform")
	{
		int rand_int_0_to_max = rand();
		int rand_int_unif_in_range = minVal + (rand_int_0_to_max % (int)(maxVal - minVal + 1));
		return rand_int_unif_in_range;
	}
	else
	{
		//at the moment - no other implementation
		return -1;
	}
}

double roundParamVal(double paramToRound, int numDigAfterDecPoint)
{
	double factor = pow(10.0,numDigAfterDecPoint);
	double roundedVal = (int)(paramToRound * (int)factor) / factor;
	return roundedVal;
}