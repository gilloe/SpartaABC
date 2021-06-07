#ifndef _PARAM_PRIORS
#define _PARAM_PRIORS

#include <string>
#include <random>
#include <cstdlib>
#include <math.h>

using namespace std;
double getRandDoubleParamVal(string distributionName, double minVal, double maxVal);
int getRandIntParamVal(string distributionName, int minVal, int maxVal);
double roundParamVal(double paramToRound, int numDigAfterDecPoint);

#endif