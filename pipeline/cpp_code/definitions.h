#ifndef ___DEFINITIONS_H
#define ___DEFINITIONS_H

#include <vector>
#include <string>
#include <limits>
using namespace std;

// Sometimes you want to change all the double variables in your program to floats during compilation.
// In this case, you can just change the following #define.
#define MDOUBLE double
//#define MDOUBLE float


// Constants
#define PI  (3.1415926535897932384626433832795028841971693993751058)

// Instead of writing vector<vector< int>> it is easier to write VVint.
typedef vector<MDOUBLE> Vdouble;
typedef vector<int> Vint;
typedef vector<Vint> VVint;
typedef vector<VVint> VVVint;
typedef vector<char> Vchar;
typedef vector<Vdouble> VVdouble;
typedef vector<VVdouble> VVVdouble;
typedef vector<VVVdouble> VVVVdouble;
typedef vector<VVVVdouble> VVVVVdouble;
typedef vector<string> Vstring;

// Definitions of very small and very big numbers.
const MDOUBLE VERYBIG = (std::numeric_limits<MDOUBLE>::max)();
const MDOUBLE VERYSMALL = -VERYBIG;
const MDOUBLE EPSILON = numeric_limits<MDOUBLE>::epsilon(); //epsilon() returns the difference between 1 and the smallest value greater than 1 that is representable for the data type.


#endif
  
	
