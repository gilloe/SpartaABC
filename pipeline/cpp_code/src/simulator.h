#ifndef _SIMULATOR_H
#define _SIMULATOR_H

#include<vector>
#include <string>
#include<chrono>
#include "tree.h"
#include "FastZip.h"
#include "SpartaABC_options.h"
#include "Block.h"
#include "Event.h"
#include "Node.h"
#include "LinkedList.h"

using namespace std;

typedef unsigned short nucID;
const int MaxIndelSize = 50; // maximum length of insertion or deletion

class Simulator {
	public:
		explicit Simulator(size_t rootLength, double A_Insertion_param,double A_Deletion_param, double IR, double DR) : _rootLength(rootLength), _A_Insertion_param(A_Insertion_param), _A_Deletion_param(A_Deletion_param), _IR(IR), _DR(DR), _fastZInsertions(A_Insertion_param, MaxIndelSize), _fastZDeletions(A_Deletion_param, MaxIndelSize) {};
		vector<string> simulateBasedOnTree(string treeFileName);
		vector<string> leafNames;
	private:
		//returning bool since 30.3 //OI //temp return void- will ba back soon
		bool simulateAlongTree(tree::nodeP t,
			 const vector<Node*>& father_seq,//change- added "const vector<Node*> father_seq" //17.3 OI
			 vector<vector<Node*>>& upgeaded_simulated_sequences);//another change- added "vector<vector<Node*>>& upgeaded_simulated_sequences" //OI 17.3
		//new 7.4
		bool simualteWithIndelsAlongAspecificBranch_while_creating_events(double branchlength,
			const vector<Node*>& ancestralsequence, vector<Node*>& output_sequence);

		vector<Node*> generateRootSequence(); //new OI 17.3

		void back_to_numbers(const vector<vector<Node*>>& upgraded_vector, vector<vector<nucID>>& numbered_vector);//new OI 17.3
		void simulationToMSA_with_SuperSequence(const vector<vector<Node*> >& simulatedLeavesSequences, vector<string>& msa);//new OI 17.3
		void simulationToMSA_new_with_changes(const vector<vector<Node*> >& simulatedLeavesSequences, vector<string>& msa);//new 11.4
		void printSequences(const vector<vector<nucID> >& v);


		void printSequenceold(const vector<nucID>& v);
		void printSequencenew(const vector<Node*>& v);
		void printblockvec(vector<Block>& vb);
		void outputInFastaFormat(const vector<Node*>& leafSequence);
		void outputInFastaFormat(const vector<vector<Node*> >& simulatedLeavesSequences);

		/***
			 OLD METHODS
						***/
		bool isequal(const vector<nucID>& v, const vector<Node*>& up_v);
		bool isequal(const vector<string>& v, const vector<string>& up_v);
		void simualteWithIndelsAlongAspecificBranch(const vector<nucID>& ancestralSequence,
			double branchLength, vector<nucID>& outputSeq);
		void simualteWithIndelsAlongAspecificBranch2(const vector<Event>& eventvec, const vector<Node*>& ancestralsequence
			, vector<Node*>& output_sequence);//new OI 15.3
		bool create_event_vector(double branchlength, int ancestralSequencelength, vector<Event>& event_vec);//new OI 15.3
		vector<string> change_to_new_msa(vector<string>& msa);
		//from 30.3 retunring boolean
		void simualteWithIndelsAlongAspecificBranch_with_events_not_rec(
			const vector<nucID>& ancestralSequence,
			double branchLength,
			vector<nucID>& outputSeq, vector<Event>& event_vec);
	


		//VARIBLES
		size_t _rootLength;
		double _A_Insertion_param; // this parameter dictates the Zippfian distribution of insertion sizes
		double _A_Deletion_param; // this parameter dictates the Zippfian distribution of deletion sizes
		double _IR; // note that this is the insertion rate per position
		double _DR; // note that this is the deletion rate per position

		FastZip _fastZInsertions; // generates random number according to a zippfian distribution
		FastZip _fastZDeletions;
};
// The simulator simulates sequences along a tree without substitutions
// The output is a list of strings, each corresponding to a single sequence
// The name of each sequence is not given in this output.
// For example, if we simulate along the tree ((s1:0.3; s2: 0.4), s3:0.5);
// We will get as output a vector of string, say, V, so that:
// v[0] = "ACGAAG-"
// v[1] = "A---AG-"
// v[2] = "-CGAAGA"
// In truth, we do not simulate substitutions, so that vectors will actually look like:
// v[0] = "AAAAAA-"
// v[1] = "A---AA-"
// v[2] = "-AAAAAA"

#endif