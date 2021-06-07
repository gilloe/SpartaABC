#ifndef _SIM_MATRICES
#define _SIM_MATRICES

#include <map>
#include <ctype.h>
#include <climits>
#include <string>
using namespace std;

class sim_matrices
{
	
public:
	sim_matrices(int mode, int default_match = (INT_MAX/2), int default_mismatch = - (INT_MAX/2)): _default_match(default_match), _default_mismatch(default_mismatch)
	{
		if(mode == 0) //simple mode
		{
			//no need to initialize
		}
		else if(mode == 1) //nuc mode
		{
			initNucMap();
		}
		else
		{
			//TO DO - PAM250
		}
		
	}
	sim_matrices(){};
	~sim_matrices(){};

	//typedef map<pair<char,char>,int> similarityMap;
	//static similarityMap _similarityMap;    

	static void initNucMap();
	int get_pair_score(char a, char b);

private:
	static map<pair<char,char>,int> _similarityMap;
	int _default_match;
	int _default_mismatch;

};




#endif