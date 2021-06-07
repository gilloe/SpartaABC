#include "errorMsg.h"
#include "logFile.h"

#include <cassert>
using namespace std;

ostream *errorMsg::_errorOut= NULL;

void errorMsg::reportError(const vector<string>& textToPrint, const int exitCode) {
	for (size_t i =0 ; i < textToPrint.size() ; ++i) {
		if (&myLog::LogFile() != &cerr) {
			LOG(1, << textToPrint[i] << endl);
		}
		cerr<<textToPrint[i]<<endl;
		if (!(_errorOut) && _errorOut != &cerr)  {
			(*_errorOut)<<textToPrint[i]<<endl;
		}
	}
	assert(0); // always stop here if in DEBUG mode.
	exit(exitCode);
}

void errorMsg::reportError(const string& textToPrint, const int exitCode) {
	if (&myLog::LogFile() != &cerr) {
		LOG(1, << endl << textToPrint << endl);
	}
	cerr<<endl<<textToPrint<<endl;
	if (!(_errorOut) && _errorOut != &cerr)  {
		(*_errorOut)<<textToPrint<<endl;
	}
	assert(0); // always stop here if in DEBUG mode.
	exit(exitCode);
}


