#ifndef _SIMULATION
#define _SIMULATION

#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>

#include "MSA.h"
//#include "indelible.h"

// Includes for using Dawg:
#ifdef WITH_DAWG
	#include <dawg/trick.h>
	#include <dawg/trick_parse.h>
	#include <dawg/global.h>
	#include <dawg/output.h>
	#include <dawg/ma.h>
	#include <dawg/matic.h>
#endif
// End includes for using Dawg

#ifdef _WIN32
	//define something for Windows (32-bit and 64-bit, this part is common)
	#include <windows.h>
#endif

using namespace std;

extern  int simNo;
class Simulation
{
public:
	Simulation(double indelRate, int rootLength, double indelDistributionShapeParameter): _indelRate(indelRate), _rootLength(rootLength),_indelDistributionShapeParameter(indelDistributionShapeParameter)  
	{	
		_simulationIdentified=simNo;
		_MaxdeletionLengh = 50;
		_dist = 0;
		
	};
	void simulateMSA(const int numberOfSimulations, const vector<string> & templateInstructionString);
	void generateMSA_array(int numberOfSimulations,const vector<char> & modifiedInstructionChars);
	
	double getSimulationDist(){return _dist;}
	double get_indelRate(){return _indelRate;}
	int get_rootLength(){return _rootLength;}
	double get_indelDistributionShapeParameter(){return _indelDistributionShapeParameter;}
	string get_indelOutputFileName(){return _indelOutputFileName;}
	int get_simulationIdentified(){return _simulationIdentified;};
	
	void set_indelDistributionShapeParameter(double x){_indelDistributionShapeParameter = x;}
	void setIndelRate(double x){_indelRate=x;}
	void setrootLength(int x){_rootLength=x;}
	void setSimulationProperties(double indelRate, int rootLength, double indelDistributionShapeParameter, int simulationIdentified);
	void set_simulationIdentified(int sn){_simulationIdentified = sn;}
	void setDistance(double d){_dist=d;};
	Simulation& operator = (const Simulation& );
	bool operator<( Simulation& rhs);
	vector<string> tempIndelibleStrToParamSpecIndelibleStr(int numberOfSimulations, const vector<string> & templateInstructionString);
	vector<char> convertVecStringToVecChar(const vector<string> & vecString);
	vector<MSA> msaVec;
	~Simulation();

	// Dawg related functions
#ifdef WITH_DAWG
	vector<string> tempDawgStrToParamSpecDawgStr(int numberOfSimulations, const vector<string> & templateInstructionString);
	void generateDawgMSA_array(int numberOfSimulations, const vector<char> & modifiedInstructionChars);
	void createDawgAlignments(std::vector<dawg::alignment>& alignments, const vector<char> & modifiedInstructionChars);
	void simulateDawgMSA(const int numberOfSimulations, const vector<string> & templateInstructionString);
#endif
	// End Dawg related functions
	



private: 
	double _indelRate;
	int _rootLength;
	int _simulationIdentified;
	double _indelDistributionShapeParameter;
	int _MaxdeletionLengh;
	string _indelOutputFileName;
	double _dist;
	
	
};

#endif