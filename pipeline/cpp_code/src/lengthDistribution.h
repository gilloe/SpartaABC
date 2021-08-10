#include <vector>
#include "definitions.h"
#include "SpartaABC_options.h"


class lengthDistribution
{
public:
	lengthDistribution(int max): maxLength(max) {};
	virtual void generateLengthDistribution() = 0;
	virtual int get_length() = 0;
protected:
	int maxLength;
};