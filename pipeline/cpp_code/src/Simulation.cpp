#include "Simulation.h"

void Simulation::setSimulationProperties(double indelRate, int rootLength, double indelDistributionShapeParameter, int simulationIdentified)
{
	_indelRate = indelRate;
	_rootLength = rootLength;
	_indelDistributionShapeParameter = indelDistributionShapeParameter;
	_simulationIdentified = simulationIdentified;
	_MaxdeletionLengh = 50;
	_dist = 0;
}
void Simulation::simulateMSA(const int numberOfSimulations, const vector<string> & templateInstructionString) 
{
	vector<string> modifiedInstructionString = tempIndelibleStrToParamSpecIndelibleStr(numberOfSimulations, templateInstructionString);
	vector<char> modifiedInstructionChars = convertVecStringToVecChar(modifiedInstructionString);
	generateMSA_array(numberOfSimulations,modifiedInstructionChars);
}

vector<char> Simulation::convertVecStringToVecChar(const vector<string> & vecString)
{
	vector<char> vecChar;
	for(int i=0; i<vecString.size(); i++)
	{
		string currLine = vecString[i];
		for(int j=0; j<currLine.size(); j++)
		{
			char currChar = currLine[j];
			vecChar.push_back(currChar);
		}
		vecChar.push_back('\n');
	}
	return vecChar;
}

vector<string> Simulation::tempIndelibleStrToParamSpecIndelibleStr(int numberOfSimulations, const vector<string> & templateInstructionString)
{
	vector<string> modifiedInstructionString;
	for(int i=0; i<templateInstructionString.size(); i++) 
	{
  		string line = templateInstructionString[i];
		if(line.find("indelmodel")!=std::string::npos) 
		{
			string newLine = "  [indelmodel]  POW  "+to_string(static_cast<long double>(_indelDistributionShapeParameter))+" "+to_string(static_cast<long long>(_MaxdeletionLengh))+"\n";
			modifiedInstructionString.push_back(newLine);
		}
		else if(line.find("indelrate")!=std::string::npos )	
		{
			string newLine = "  [indelrate]   "+to_string(static_cast<long double>(_indelRate))+"\n";
			modifiedInstructionString.push_back(newLine);
		}
		else if(line.find("treename modelname")!=std::string::npos ) 
		{
			string newLine = "  [treename modelname  "+to_string(static_cast<long long>(_rootLength))+"]\n";
			modifiedInstructionString.push_back(newLine);
		}			
		else if	(line.find("EVOLVE")!=std::string::npos ) 
		{
			string newLine = "[EVOLVE] partitionname 1 "+to_string(static_cast<long long>(_simulationIdentified))+"\n";
			modifiedInstructionString.push_back(newLine);			
		}
		else
		{
			modifiedInstructionString.push_back(line);
		}
	}
	return modifiedInstructionString;
}


//input default parameter of file name 
void Simulation::generateMSA_array(int numberOfSimulations,const vector<char> & modifiedInstructionChars) 
{
	cout << "this function is no longer supported" << endl;
	exit(1);
	/*
	msaVec.clear();
	vector<string> simulSeq;
	for (int i=0; i < numberOfSimulations; ++i)
	{
		indelible(modifiedInstructionChars, simulSeq);
		MSA msa(simulSeq);
		msaVec.push_back(msa);
		simulSeq.clear();
	}*/
}

Simulation&  Simulation::operator=(const Simulation& other)
{
	_dist=other._dist;
	_indelDistributionShapeParameter = other._indelDistributionShapeParameter;
	_indelOutputFileName = other._indelOutputFileName;
	_indelRate= other._indelRate;
	_MaxdeletionLengh =other._MaxdeletionLengh;
	msaVec = other.msaVec;
	_rootLength = other._rootLength;
	_simulationIdentified = other._simulationIdentified;    
	return *this;
}

bool Simulation::operator<( Simulation& rhs)
{
	if (this->_dist< rhs.getSimulationDist())
			return true;
		return false;
}


Simulation::~Simulation()
{
	msaVec.clear();

}


// Dawg related functions:
#ifdef WITH_DAWG
vector<string> Simulation::tempDawgStrToParamSpecDawgStr(int numberOfSimulations, const vector<string> & templateInstructionString)
{
	vector<string> modifiedInstructionString;
	for (int i = 0; i<templateInstructionString.size(); i++)
	{
		string line = templateInstructionString[i];
		if (line.find("Params.Ins") != std::string::npos)
		{
			string newLine = "Params.Ins = " + to_string(static_cast<long double>(_indelDistributionShapeParameter)) + ", " + to_string(static_cast<long long>(_MaxdeletionLengh)) + "\n";
			modifiedInstructionString.push_back(newLine);
		}
		else if (line.find("Params.Del") != std::string::npos)
		{
			string newLine = "Params.Del = " + to_string(static_cast<long double>(_indelDistributionShapeParameter)) + ", " + to_string(static_cast<long long>(_MaxdeletionLengh)) + "\n";
			modifiedInstructionString.push_back(newLine);
		}
		else if (line.find("Rate.Ins") != std::string::npos)
		{
			string newLine = "Rate.Ins = " + to_string(static_cast<long double>(_indelRate)) + "\n";
			modifiedInstructionString.push_back(newLine);
		}
		else if (line.find("Rate.Del") != std::string::npos)
		{
			string newLine = "Rate.Del = " + to_string(static_cast<long double>(_indelRate)) + "\n";
			modifiedInstructionString.push_back(newLine);
		}
		else if (line.find("Length") != std::string::npos)
		{
			string newLine = "Length = " + to_string(static_cast<long long>(_rootLength)) + "\n";
			modifiedInstructionString.push_back(newLine);
		}
		else
		{
			modifiedInstructionString.push_back(line);
		}
	}
	return modifiedInstructionString;
}

void Simulation::simulateDawgMSA(const int numberOfSimulations, const vector<string> & templateInstructionString)
{
	vector<string> modifiedInstructionString = tempDawgStrToParamSpecDawgStr(numberOfSimulations, templateInstructionString);
	vector<char> modifiedInstructionChars = convertVecStringToVecChar(modifiedInstructionString);
	generateDawgMSA_array(numberOfSimulations, modifiedInstructionChars);
}

void Simulation::generateDawgMSA_array(int numberOfSimulations, const vector<char> & modifiedInstructionChars)
{

	msaVec.clear();
	// The alignments created from DAWG
    std::vector<dawg::alignment> alignments;
	for (int i = 0; i < numberOfSimulations; ++i)
	{
		createDawgAlignments(alignments, modifiedInstructionChars);
		dawg::alignment single_dawg_alignment = alignments[0];
		vector<string> dawg_aligned_seqs;
		for (int i = 0; i < single_dawg_alignment.size(); i++)
		{
			//string name = single_dawg_alignment[i].label;
			//string seq = single_dawg_alignment[i].seq;
			//cout << name << "\t" << seq << endl;

			dawg_aligned_seqs.push_back(single_dawg_alignment[i].seq);
		}

		MSA dawg_sim_MSA(dawg_aligned_seqs);
		msaVec.push_back(dawg_sim_MSA);
		alignments.clear();
	}
}

void Simulation::createDawgAlignments(std::vector<dawg::alignment>& alignments, const vector<char> & modifiedInstructionChars)
{
	dawg::trick input;

	input.parse(modifiedInstructionChars.begin(), modifiedInstructionChars.end());

	// process aliases
	input.read_aliases();

	dawg::global_options glopts;
	glopts.read_section(input.data.front());

	dawg::output write_aln;

	// Since no output file has been specified,
	// the output will go to std::cout
	if (!write_aln.open(/*glopts.output_file.c_str()*/ "aln:-",
		glopts.sim_reps - 1,
		glopts.output_split, // false
		glopts.output_append, // false
		glopts.output_label)) // false
	{
		DAWG_ERROR("bad configuration");
		return;
	}
	write_aln.set_blocks(glopts.output_block_head.c_str(),
		glopts.output_block_between.c_str(),
		glopts.output_block_tail.c_str(),
		glopts.output_block_before.c_str(),
		glopts.output_block_after.c_str()
	);

	std::vector<dawg::ma> configs;
	if (!dawg::ma::from_trick(input, configs)) {
		DAWG_ERROR("bad configuration");
		return;
	}

	// Create the object that will do all the simulation
	// work for us.  Configure its sections.
	dawg::matic kimura;

	if (!kimura.configure(configs.begin(), configs.end())) {
		DAWG_ERROR("bad configuration");
		return;
	}

	// create sets of aligned sequences;
	dawg::alignment aln;
	kimura.pre_walk(aln);
	for (unsigned int i = 0; i<glopts.sim_reps; ++i) {
		kimura.walk(aln);
		alignments.insert(alignments.end(), aln);
		//write_aln(aln); // this would print the aln data out to std::cout or a file
	}
}
#endif
//End Dawg related functions


