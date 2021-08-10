#include "lengthDistribution.h"


class geometricDist: public lengthDistribution 
{
public:
	geometricDist(int max): lengthDistribution(max) {};
	void generateLengthDistribution();
	int get_length();
private:
	Vint bins;
};


