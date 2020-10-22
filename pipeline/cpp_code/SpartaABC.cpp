//
//  main.cpp
//  SpartaABC
//	Integrate ABC concepts to the SPARTA algorithm
//	
#include "simulator.h"
#include "SpartaABC.h"
#include "FastZip.h"


int simNo=0;
size_t totalNumberOfSimulations = 0;

vector<double> simulateSequencesAndReturnSummaryStatistics(size_t randRootLength,
	double randAInsertionParam,
	double randADeletionParam,
	double randInsertRatio,
	double randDeletionRatio) {
	Simulator sim(randRootLength, randAInsertionParam, randADeletionParam, randInsertRatio, randDeletionRatio);
	vector<string> simulatedSeqs = sim.simulateBasedOnTree(SpartaABC_options::_inputTreeFileName);// first of all- notice the change- now we get back vector<vector<string>>
	if (simulatedSeqs[0] == "too long")
	{
		vector<double> too_long;
		too_long.push_back(-1);
		return too_long;
	}
	MSA simulated_MSA(simulatedSeqs);// upgraded vector - my MSA
	//summ_stats_wrapper simulatedSummStats(simulated_MSA, SpartaABC_options::_alignmentMode, SpartaABC_options::_similarityMode);
	//vector<double> summStatsSim = simulatedSummStats.getSummStatsValsVector();
	//cout << "is equal: " << is_equal(simulated_MSA_original.getAlignedSeqs(), simulated_MSA.getAlignedSeqs()) << endl;
	
	vector<double> summStatsSim = getStatVec(simulated_MSA);// upgraded summary sttatistics.
	return summStatsSim;



}

int main(int argc, const char * argv[])
{
	//begining of the main
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();// time generator
	srand((unsigned)(time(0))); // for random number generation. Do not delete.

	if (argc < 2) {// argc is 1 (the name of the program) plus the number of arguments
		cout<<"SpartaABC requires a command file as an argument and no such file was provided."<<endl;
		exit(1);
	}
	
	string spartaABCConfigFile = argv[1];
	SpartaABC_options::initOptions(spartaABCConfigFile);
	vector<double> summStatsReal = getStatVec(SpartaABC_options::_inputRealMSAFile);
	
	// in case summary statics file is inserted as argument
	ifstream ss_file;
	if (argc > 2)
	{
		cout << " SpartaABC got an SumStats file" << endl;
		string ss_path = argv[2];
		
		ss_file.open(ss_path);
		if (ss_file.fail()||!ss_file)
		{
			cout << "can not open file:" << ss_path << endl;
			exit(1);
		}
		if (ss_file.good())
		{
			
			vector<double> new_sumstats;
			double ss_double;
			string ss;
			for (size_t i = 0; i < summStatsReal.size(); i++)
			{
				if (getline(ss_file, ss))
				{
					ss_double = stod(ss);
					new_sumstats.push_back(ss_double);
				}
			}
			summStatsReal = new_sumstats;
		}
	
		

	}
	SpartaABC_options::setRLPrior(static_cast<int>(summStatsReal[3]), static_cast<int>(summStatsReal[2]));


	//open result file where the parameter combinations that yielded small enough distance are kept
	ofstream resFile;
	resFile.open(SpartaABC_options::_outputGoodParamsFile);
	if (!resFile) 
	{
       cout<<"couldn't open output file: "<<SpartaABC_options::_outputGoodParamsFile<<endl;
       exit(1);
	}
	//read template control of Indelible/Dawg into vector<string>
	vector<string> templateCtrl;
	if (SpartaABC_options::_dawgSimulator)
	{
		templateCtrl = readIndelibleTemplateControlFile(SpartaABC_options::_dawgTemplateControlFile);
	}
	else
	{
		templateCtrl = readIndelibleTemplateControlFile(SpartaABC_options::_indelibleTemplateControlFile);
	}

	int numberOfCollectedSamples = 0;
	int numberSimulations = 0;
	//added my vector of res //OI 29.3
	vector<vector<double>> allSimulationsResults;

	

	//in order to use weights
	//cout << "computing weights" << endl;
	vector<double>  weightStats = getWeightsVector();
	euclidean_distance euclideanDistObj;
	euclideanDistObj.setWeightsVector(weightStats);

	string niceHeader = NiceHeader();
	resFile << niceHeader << endl;

	//print statistics for input MSA
	resFile << "input_msa\t?\t?\t?\t?\t?";
	for (size_t i = 0; i < summStatsReal.size(); i++)
	{
		resFile << "\t" << summStatsReal[i];
	}

	resFile << endl;
	resFile << "weights\t\t\t\t\t";
	for (size_t i = 0; i < weightStats.size(); i++)
	{
		resFile<<"\t" << weightStats[i];
	}
	resFile << endl;


	double currentCutOffValue = 1000.0; // very large
	size_t totalNumberOfSpartaTests = 100000;//temp ,suppose to be 100,000
	while(totalNumberOfSimulations < totalNumberOfSpartaTests)
	{
		totalNumberOfSimulations++;
		int randRootLength = getRandIntParamVal(SpartaABC_options::_priorDistTypeRL, SpartaABC_options::_minRLVal, SpartaABC_options::_maxRLVal);
		double randAInsertionParam = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeA, SpartaABC_options::_minAVal, SpartaABC_options::_maxAVal);
		double randADeltetionParam = randAInsertionParam;//getRandDoubleParamVal(SpartaABC_options::_priorDistTypeA, SpartaABC_options::_minAVal, SpartaABC_options::_maxAVal);
		double randInsertRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
		double randDeletionRatio = randInsertRatio;//getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
		while (randInsertRatio == 0/* && randDeletionRatio == 0*/) {
			 randInsertRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
			 randDeletionRatio = randInsertRatio;// getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
		}
		if (numberSimulations % 1000 == 0) {
			cout << "Simulation number:\t" << numberSimulations << endl;// "\tKept thus far:\t" << numberOfCollectedSamples << "\tParams are:\tRL " << randRootLength << "\tA " << randAParam << "\tInsR " << randInsertRatio << "\tDelR " << randDeletionRatio << endl;
		}
		//MSA simulated_MSA = simulateSingleMSA(randRootLength, randAParam, randIndelRatio, templateCtrl);

		// our own version to get an MSA
		// 1. simulate along the tree
		//SIMULATION!!!
		vector<double> summStatsSim=simulateSequencesAndReturnSummaryStatistics(randRootLength, randAInsertionParam, randADeltetionParam,
			randInsertRatio, randDeletionRatio);
			
		double euclideanDist;
		//initilizing very large parameter if it exploded
		if (summStatsSim[0] == -1)
		{
			summStatsSim[0] = 10000000;
			euclideanDist = 100000; //very big value
		}
		else
		{
			//adjusting the results vector.
			euclideanDist = euclideanDistObj.computeWeightedDist(summStatsReal, summStatsSim);
			if (euclideanDist < currentCutOffValue)
			{
				vector<double> tmp;
				tmp.push_back(randRootLength);
				tmp.push_back(randAInsertionParam);
				tmp.push_back(randADeltetionParam);
				tmp.push_back(randInsertRatio);
				tmp.push_back(randDeletionRatio);
				for (size_t j = 0; j < summStatsSim.size(); ++j) {
					tmp.push_back(summStatsSim[j]);
				}
				tmp.push_back(euclideanDist); // not a parameter
				size_t k = 0;
				for (; k < allSimulationsResults.size(); ++k) {
					if (allSimulationsResults[k][tmp.size() - 1] > euclideanDist) break;
				}
				allSimulationsResults.insert(allSimulationsResults.begin() + k, tmp);
				if (allSimulationsResults.size() > SpartaABC_options::_numberOfSamplesToKeep) {
					allSimulationsResults.pop_back();
					currentCutOffValue = allSimulationsResults[allSimulationsResults.size() - 1][tmp.size() - 1];
				}
			}
			numberSimulations++;
		}
		// for Gil :
		resFile <<euclideanDist<<"\t" <<randRootLength << '\t' << randAInsertionParam << '\t' << randADeltetionParam << '\t' << randInsertRatio << '\t' << randDeletionRatio;
		for (size_t l = 0; l < summStatsSim.size(); ++l) {
			resFile << "\t" << summStatsSim[l];
		}
		resFile << "\n";
		
	}	
	

	
	//calculating the posterior params.
	double posteriorRL = 0.0;
	double posteriorInsR = 0.0;
	double posteriorDelR = 0.0;
	double posteriorAInsertion = 0.0;
	double posteriorADeletion = 0.0;
	//resFile << "\n\n";
	for (size_t i = 0; i < allSimulationsResults.size(); ++i) {
	/*	resFile << allSimulationsResults[i][allSimulationsResults[i].size()-1] << "\t" << allSimulationsResults[i][0] << "\t" << allSimulationsResults[i][1]
			<< "\t" << allSimulationsResults[i][2] << "\t" << allSimulationsResults[i][3] << "\t" << allSimulationsResults[i][4];
		for (size_t l = 5; l < allSimulationsResults[i].size()-1; ++l) {
			resFile << "\t" << allSimulationsResults[i][l];
		}*/
		//resFile << "\n";
		posteriorRL += allSimulationsResults[i][0];
		posteriorAInsertion += allSimulationsResults[i][1];
		posteriorADeletion += allSimulationsResults[i][2];
		posteriorInsR += allSimulationsResults[i][3];
		posteriorDelR += allSimulationsResults[i][4];

	}
	posteriorRL /= allSimulationsResults.size();
	posteriorInsR /= allSimulationsResults.size();
	posteriorDelR /= allSimulationsResults.size();
	posteriorAInsertion /= allSimulationsResults.size();
	posteriorADeletion /= allSimulationsResults.size();
	resFile << endl;
	resFile << "The posterior expectation for the root length (RL) is: " << posteriorRL << endl;
	resFile << "The posterior expectation for the A Insertion parameter is: " << posteriorAInsertion << endl;
	resFile << "The posterior expectation for the A Deletion parameter is: " << posteriorADeletion << endl;
	resFile << "The posterior expectation for the indel rate (InsR) is: " << posteriorInsR << endl;
	resFile << "The posterior expectation for the indel rate (DelR) is: " << posteriorDelR << endl;
	
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	/*resFile << "The simulation running time is: "<< std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;*/
	resFile << "Total number of simulations: " << numberSimulations << endl;
	resFile.close();
	ss_file.close();
    return 0;
}

vector<string> readIndelibleTemplateControlFile(string indelibleTemplateCtrlFile)
{
	vector<string> templateInstructionString;
	ifstream InFile;
	InFile.open(indelibleTemplateCtrlFile); 
	if(!InFile.is_open()) 
	{
		cout<<"can't open control file "<<indelibleTemplateCtrlFile<<endl;
		exit(1);
	}
	string line;
	while(! InFile.eof()) 
	{
  		getline(InFile,line);
		templateInstructionString.push_back(line);
	}
	return templateInstructionString;
}


MSA simulateSingleMSA(int rootLength, double a_param, double indel_rate_ratio, vector<string> templateCtrl)
{

	Simulation sim(indel_rate_ratio, rootLength, a_param);
	// if compilation is with Dawg (still user can choose between INDELible and Dawg)
	#ifdef WITH_DAWG
	if (SpartaABC_options::_dawgSimulator)
	{
		sim.simulateDawgMSA(1, templateCtrl);
	}
	else
	{
		sim.simulateMSA(1, templateCtrl);
	}
	#elif WITH_INDELIBLE

	// compilation is without Dawg
	sim.simulateMSA(1, templateCtrl);
	#else // USING OUR OWN SIMULATOR
	sim.simulateMSA(1, templateCtrl);
	#endif

	MSA simulatedMSA = sim.msaVec[0];
	return simulatedMSA;
}

