#pragma once
#include <vector>
#include <iostream>
#include "param_priors.h"
#include <cmath>

using namespace std;

class FastZip
{
public:
	explicit FastZip(double aParam, int max);
	~FastZip();
	int drawZip();

private:
	int numBins;
	vector<int> lowVector;
	vector<int> highVector;
	vector<double> lowHighRateVector;
};
