#include "sim_matrices.h"

map<pair<char,char>,int> sim_matrices::_similarityMap;

void sim_matrices::initNucMap()
{
	_similarityMap[make_pair('a','a')] = 10;
	_similarityMap[make_pair('a','c')] = -10;
	_similarityMap[make_pair('a','g')] = -5;
	_similarityMap[make_pair('a','t')] = -10;

	_similarityMap[make_pair('c','a')] = -10;
	_similarityMap[make_pair('c','c')] = 10;
	_similarityMap[make_pair('c','g')] = -10;
	_similarityMap[make_pair('c','t')] = -5;

	_similarityMap[make_pair('g','a')] = -5;
	_similarityMap[make_pair('g','c')] = -10;
	_similarityMap[make_pair('g','g')] = 10;
	_similarityMap[make_pair('g','t')] = -10;

	_similarityMap[make_pair('t','a')] = -10;
	_similarityMap[make_pair('t','c')] = -5;
	_similarityMap[make_pair('t','g')] = -10;
	_similarityMap[make_pair('t','t')] = 10;

}


int sim_matrices::get_pair_score(char a, char b)
{
	char lower_a = tolower(a);
	char lower_b = tolower(b);
	pair<char,char> curr_pair = make_pair(a,b);
	
	if(_similarityMap.find(curr_pair) == _similarityMap.end())
	{
		// not in map
		if(a == b)
		{
			return _default_match;
		}
		else
		{
			return _default_mismatch;
		}
	}
	// exists in map
	int value = _similarityMap[curr_pair];
	return value;

}