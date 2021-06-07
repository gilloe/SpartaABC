#include "danaRandomGenerators.h"
#include <iostream>

int seed = static_cast<int>(chrono::system_clock::now().time_since_epoch().count());
default_random_engine generator(seed);
mt19937 mt_rand(seed);

using namespace std;
/*
int drawZipfIndelible(double alpha, int n) {
	zipf_distribution<int, double> distribution(n, alpha);
	int number = distribution(mt_rand);
	return number;
}

int drawZipf(double alpha, int n) {
	drawZipfIndelible(alpha, n);
}
*/
double drawExp(double lambda) {
	
	exponential_distribution<double> distribution(lambda);
	double number = distribution(generator);
	return number;
}

double uniform() { // uniform between 0 and 1
	uniform_real_distribution<double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return number;
}

int uniform(int a, int b) { // a random number between a and b, including a and b.
	uniform_int_distribution<int> distribution(a, b);
	int number = distribution(generator);
	return number;
}

int Zipf(double q, double v);
int powerLaw(double alpha, int n)
{// very tmp code
	int tmp = Zipf(alpha,1);
	while (tmp > n) tmp = Zipf(alpha,1);
	return tmp;

	static bool first = true;      // Static first time flag
	static double c = 0;          // Normalization constant
	static double *sum_probs;     // Pre-calculated sum of probabilities
	double z;                     // Uniform random number (0 < z < 1)
	int zipf_value;               // Computed exponential value to be returned
	int    i;                     // Loop counter
	int low, high, mid;           // Binary-search bounds

								  // Compute normalization constant on first call only
	if (first == true)
	{
		for (i = 1; i <= n; i++)
			c = c + (1.0 / pow((double)i, alpha));
		c = 1.0 / c;

		sum_probs = (double*)malloc((n + 1) * sizeof(*sum_probs));
		sum_probs[0] = 0;
		for (i = 1; i <= n; i++) {
			sum_probs[i] = sum_probs[i - 1] + c / pow((double)i, alpha);
		}
		first = false;
	}

	// Pull a uniform random number (0 < z < 1)
	do
	{
		z = uniform();
	} while ((z == 0) || (z == 1));

	// Map z to the value
	low = 1, high = n, mid;
	do {
		mid = (int)floor((low + high) / 2);
		if (sum_probs[mid] >= z && sum_probs[mid - 1] < z) {
			zipf_value = mid;
			break;
		}
		else if (sum_probs[mid] >= z) {
			high = mid - 1;
		}
		else {
			low = mid + 1;
		}
	} while (low <= high);

	// Assert that zipf_value is between 1 and N
	assert((zipf_value >= 1) && (zipf_value <= n));

	return(zipf_value);
}
