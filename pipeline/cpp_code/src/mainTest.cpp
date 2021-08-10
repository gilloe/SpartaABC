#include <iostream>
#include "geometricDist.h"

using namespace std;

int main1(){

	geometricDist g = geometricDist(50);

	g.generateLengthDistribution();

	cout << g.get_length() << endl;
	cout << g.get_length() << endl;
	cout << g.get_length() << endl;
	cout << g.get_length() << endl;
	cout << g.get_length() << endl;
	cout << g.get_length() << endl;
	cout << g.get_length() << endl;

	return 0;
}