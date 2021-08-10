#include "geometricDist.h"
#include "param_priors.h"


void geometricDist::generateLengthDistribution() {
	MDOUBLE r = getRandDoubleParamVal("uniform", 0.0, 1.0);

	
	MDOUBLE max = 10000.0;

	int numToPush = 0;
	while (numToPush < this->maxLength) {
		int amountToPush = r*max;
		max = max - max*r;
		numToPush++;
		for (size_t i = 0; i < amountToPush; i++) {
			bins.push_back(numToPush);
		}
	}

	
}


int geometricDist::get_length() {
	// MDOUBLE r = getRandDoubleParamVal("uniform", 0.0, 1.0);

	int binChoice = getRandIntParamVal("uniform", 0, this->bins.size());
	// geometric_distribution<> d(r);
	// d(rand)
	return bins[binChoice];
}