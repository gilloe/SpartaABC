#ifndef ___ERROR_MSG_H
#define ___ERROR_MSG_H

#include <string>
#include <vector>
#include <iostream>
using namespace std;

// The error will always be sent to cerr.
// The error will  also be sent to the log file, with level 1, unless it is set to cerr.
// There is a way to output the error also to a file (or to a differert stream).
// The way to do it is to set the value of _errorOut.
// By defaul, it is set to NULL, but one can change it using setErrorOstream().

class errorMsg {
public:
	static void reportError(const vector<string>& textToPrint, const int exitCode=1);
	static void reportError(const string& textToPrint, const int exitCode=1);
	static void setErrorOstream(ostream* errorOut) {_errorOut = errorOut;}
private:
	static ostream* _errorOut;
};

// example of how to output to a file called error.txt
// ofstream f("error.txt");
// errorMsg::setErrorOstream(&f);
// errorMsg::reportError("cheers");
// f.close();
#endif

