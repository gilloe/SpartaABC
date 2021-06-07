// THE LOG CLASS IS USED TO DEBUG AND TO OUTPUT INTERMEDIATE VALUES FOR PARAMETERS DURING THE EXECUTION OF THE PROGRAM
// LAST UPDATED BY TAL PUPKO, 14 March 2017.
// A typical call using the LOG infrastructure is
// LOG(1,<<"Here is what I want to write"<<endl);
// This is reaplced using the following #define:
// #define LOG(Lev, ex) { if( Lev <= myLog::LogLevel() ) myLog::LogFile() ex; }
// to the following code:
// {
//		if (1<=myLog::LogLevel())  myLog::LogFile() <<"Here is what I want to write"<<endl;
// }
// The function LogFile() returns the _loglvl, which is by default set to 3.
// So if you send comment with values less or equal to 3 they will be sent to the log file.
// otherwise, they will be ignored.
// This mechanism allows you to select the level of details in the LOGFILE.
// The LOG mechansim allows you to direct the LOG output to either cout, cerr, or to any other output stream (ostream).
// By default, it will go to cerr

#ifndef ___LOG
#define ___LOG

#include <string>
#include <iostream>
#include <fstream>
using namespace std;
			
class myLog {
public:
	static int LogLevel() { return _loglvl;}
	static ostream& LogFile(void) {
		if (_out == NULL) return cerr;
		return *_out;
	}

	static void setLogLvl(const int newLogLvl) {_loglvl = newLogLvl;}
	static void setLogOstream(ostream* out) {_out = out;}
	static void printArgv(int loglvl, int argc, char *argv[]) ;
private:
	static ostream* _out;
	static int _loglvl;
	static bool _firstTime;
};

#endif
		

#define LOG(Lev, ex) { if( Lev <= myLog::LogLevel() ) myLog::LogFile() ex; }
#define LOGDO(Lev, ex) { if( Lev <= myLog::LogLevel() ) ex; }





