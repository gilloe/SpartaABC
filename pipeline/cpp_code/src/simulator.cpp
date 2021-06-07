#include <iostream>
#include "tree.h"
#include "treeIt.h"
#include "simulator.h"
#include "danaRandomGenerators.h"

using namespace std;

//const int MaxIndelSize = 50; // maximum length of insertion or deletion
int tmparray[MaxIndelSize];
//updated vraibles //OI 15.3
LinkedList superSequence;
vector<nucID> numbered_superSequnce;
nucID currIdToInsert;
bool is_boomed;




//void simulationToMSA(const vector<vector<nucID> > & simulatedLeavesSequences, vector<string> & msa);


//the simulation of the given tree which reprsented by "treeFileName"
vector<string> Simulator::simulateBasedOnTree(string treeFileName) {
	// In the beginning, each sequence is stored as a vector 
	// of integeres, i.e., vector<nucID> where nucID is unsigned short for example
	// a root sequence of length six would be V = {0,1,2,3,4,5}
	//remark- for each defention I will add one of my own //OI 17.3
	//another note- when I write "upgraded"- the meaning is the new code, or my code //OI  17.3
	
	//cout << "simulate on tree" << endl;
	//creating to important things- the simulated leeave vector and the current id to insert...
	vector<vector<Node*>> simulated_leaves_sequences;
	currIdToInsert = _rootLength;
	// Initating the root sequence
	//ancestral sequnce
	vector<Node*> ancestral_sequence = generateRootSequence();
	// After this code, we have a root sequence that looks like this
	// ancestralSequence = {0,1,2,3,...,N}
	
	tree t(treeFileName);
	bool length_flag=simulateAlongTree(t.getRoot(),  ancestral_sequence
		,simulated_leaves_sequences);//the updated method //0I 17.3
	
	//extreme case- when the gene length is to long.
	if (!length_flag)
	{
		vector<string> too_long_vec;
		too_long_vec.push_back("too long");
		return too_long_vec;
	}
	vector<string> resultingAlignment;


	//new alignment- the upgraded one- creating the alignment
	vector<string> resulting_alignment;
	simulationToMSA_with_SuperSequence(simulated_leaves_sequences, resulting_alignment);
	

	//outputInFastaFormat(simulatedLeavesSequences);

	//releasing memory
	superSequence.delete_LL();
	return resulting_alignment;
}

//vector<nucID> Simulator::generateRootSequence() {
//	vector<nucID> ancestralSequence;
//	for (nucID i = 0; i < _rootLength; ++i) {
//		ancestralSequence.push_back(i);
//	}
//	return ancestralSequence;
//}



void Simulator::printSequences(const vector<vector<nucID> >& v) {
	for (size_t i = 0; i < v.size(); ++i) {
		for (size_t j = 0; j < v[i].size(); ++j) {
			cout << v[i][j] << " ";
		}
		cout << endl;
	}
}

//temp funcs- added 24.3 , again in order to check equality //OI
void Simulator::printSequenceold(const vector<nucID>& v) {
	for (size_t i = 0; i < v.size(); ++i) {
			cout << v[i]<< " ";
	}
	cout << endl;
}

//temp funcs- added 24.3 , again in order to check equality //OI
void Simulator::printSequencenew(const vector<Node*>& v) {
	for (size_t i = 0; i < v.size(); ++i) {
		cout << (*v[i]).getData() << " ";
	}
	cout << endl;
}


//temp funcs- added 24.3 , again in order to check equality //OI
void Simulator::printblockvec(vector<Block>& vb)
{
	for (size_t i = 0; i < vb.size(); i++)
	{
		cout << "[" << vb[i].getIndex() << ", " << vb[i].getLength() << " ] ,";
	}
	cout << endl;
}




// In each step, this function simulate the sons of node t and if it is a leaf
// it pushes the results to the simSequences vector
bool Simulator::simulateAlongTree(tree::nodeP t,
						const vector<Node*>& father_seq,
						vector<vector<Node*>>& simulated_sequences ) {
	vector<nucID> outSequence;
	//cout << "simulate along branch" << endl;
	if (t->isLeaf()) {
		//here we added the upgraded leaves vector... //OI 17.3
		simulated_sequences.push_back(father_seq);
		leafNames.push_back(t->name());
	}
	else {// internal node}
		for (size_t i = 0; i < t->getNumberOfSons(); ++i) {
			bool flag = true;
			// first we simulate the sequence of son i,
			// the results is stored in "outSequence"


			vector <Node*> output_sequence;
			//small change- now returns boolean
			flag=simualteWithIndelsAlongAspecificBranch_while_creating_events(t->getSon(i)->dis2father(), father_seq, output_sequence);
			if (!flag)
			{
				return false;
			}
			
			flag=simulateAlongTree(t->getSon(i), output_sequence,simulated_sequences);
			//recursion!!!!
			if (!flag)
			{
				return false;
			}
		}
		
	}
	return true;
}

//the Simulator
bool Simulator::simualteWithIndelsAlongAspecificBranch_while_creating_events(double branchlength,
	const vector<Node*>& ancestralsequence, vector<Node*>& output_sequence)
{
	/*
	//strarted 15.3 OI

	simulating the events of insertions/ deletions along a specific branch in the tree.
	In the first half, we literly simulating the branch with our block vector.
	For example: ([0,100])->([0,50],[-1,5],[51,100])->([0,40],[45,50],[-1,5],[51,100]) is insertion of (51,5) and then deletion of (41,5)

	In the second half, we trnasfer the form from the block version to the node version, and we update the super sequence
	*/
	//initilizing the first block- with the size of the sequence.... //OI 15.3
	Block allsequence(0, ancestralsequence.size());
	vector<Block> block_sequence;
	block_sequence.push_back(allsequence);

	size_t sequence_length = ancestralsequence.size();
	//running over the branch...
	double sequenceWiseInsertionRate = 1.0 * _IR * (sequence_length + 1);
	double sequenceWiseDeletionRate = 1.0 * _DR * sequence_length;
	double waitingTime = drawExp(sequenceWiseInsertionRate + sequenceWiseDeletionRate);
	while (waitingTime < branchlength)
	{//running through all the events
		//cout << current_event.getIndel_Flag() << " , " << current_event.getIndex() << " , " << current_event.getLength() << endl;//in order to print events... OI 24.3	
		if (sequence_length > SpartaABC_options::_maxRLVal * 25)
		{
			// if the length is to large.
			return false;
		}
		bool isDeletion = true;
		if (uniform() < sequenceWiseInsertionRate / (sequenceWiseInsertionRate + sequenceWiseDeletionRate)) isDeletion = false;
		if (!isDeletion)
		{//Flag==TRUE---->> insertion //OI 15.3
			//initilizing the insert params //OI 15.3
			int insertion_index = uniform(0, sequence_length + 1);
			size_t insertion_length = _fastZInsertions.drawZip();
			sequence_length += insertion_length;
			//initilizing an iterator to run all over the block //comment added OI 16.3
			vector<Block>::iterator block_sequence_iterator = block_sequence.begin();
			while (block_sequence_iterator != block_sequence.end() && insertion_index >= (*block_sequence_iterator).getLength())
			{// iterating through the block sequence in order to find the right spot to insert //OI 15.3
				insertion_index -= (*block_sequence_iterator).getLength();
				++block_sequence_iterator;

			}
			Block insertion_block(-1, insertion_length);//generate new block- an insertion bloc //OI 15.3
			if (block_sequence_iterator == block_sequence.end())
			{
				//extreme condition- where the insertion is at the end //OI 15.3
				block_sequence.push_back(insertion_block);
			}
			else
			{
				if ((*block_sequence_iterator).getIndex() == -1)
				{
					// if it an insertion block- we just need to append the length of the block //OI 15.3
					(*block_sequence_iterator).setLength((*block_sequence_iterator).getLength() + insertion_length);
				}
				else
				{
					if (insertion_index == 0)
					{
						//easy condition- where we dont need to create to blocks
						block_sequence.insert(block_sequence_iterator, insertion_block);
					}
					else
					{
						// here, we generate the block we should create **after** the insertion block //OI 15.3
						Block after_insertion((*block_sequence_iterator).getIndex() + insertion_index, (*block_sequence_iterator).getLength() - insertion_index);
						//then we change the length of the original block //OI 15.3
						(*block_sequence_iterator).setLength(insertion_index);
						//now we have to different states
						if ((block_sequence_iterator != block_sequence.end()) && (next(block_sequence_iterator) == block_sequence.end()))// change of //16.3-->>18.3 // i did a wrong condition //OI
						{
							//1. this is the last bloc in the vector
							//thus we will **push back** the blocks //OI 15.3
							block_sequence.push_back(insertion_block);
							block_sequence.push_back(after_insertion);
						}
						else
						{
							//2. the block isnt the last block 
							// thus we will insert after the "block sequence iterator" //OI 15.3
							vector<Block>::iterator after_current_iterator = next(block_sequence_iterator);
							//note: change:
							//now we insert the after insertion block first, and then we point the iterator on him, and just after that we insert the insertion block... //OI 18.3
							after_current_iterator = block_sequence.insert(after_current_iterator, after_insertion);
							after_current_iterator = block_sequence.insert(after_current_iterator, insertion_block);


						}
					}
				}
			}
		}//ending insertion //OI 16.3
		else
		{//Flag==FALSE---->> deletion //OI 16.3
			//initilizing the deletion params ///OI 16.3
			int deletion_index = uniform(0, sequence_length - 1);
			size_t deletion_size = _fastZDeletions.drawZip();
			if (deletion_index + deletion_size > sequence_length)
			{// extreme cond- when we deletion length passes the length of the sequence //OI 15.3
				deletion_size = sequence_length - deletion_index;
			}
			sequence_length -= deletion_size;
			//initilizing an iterator to run all over the block //OI 16.3
			vector<Block>::iterator block_sequence_iterator = block_sequence.begin();
			while (deletion_size > 0 && block_sequence_iterator != block_sequence.end())
			{
				//now we are running on the block sequence and deketing the necessary //OI 16.3
				if ((*block_sequence_iterator).getLength() <= deletion_index)
				{
					//first step- if the current block isn't the block we need to delete //OI 16.3
					deletion_index -= (*block_sequence_iterator).getLength();//adjusting the index of the deletion
					++block_sequence_iterator;
				}
				// from here we start to delete //OI 16.3
				else
				{
					if (deletion_size >= (*block_sequence_iterator).getLength() - deletion_index)
					{

						//temp for debuging OI //28.3
						if ((*block_sequence_iterator).getIndex() == -1) {
							//indel_nums += (*block_sequence_iterator).getLength() - deletion_index;// for debuging
						}
						//first condition- if the deletion length is longer then the length of the block //OI 16.3
						deletion_size = deletion_size - (*block_sequence_iterator).getLength() + deletion_index; //adjusting the deletion size
						if (deletion_index == 0)
						{
							//extreme condition- if the deletion starts fron the 0-index (we need to delete the block)// OI 16.3

							block_sequence_iterator = block_sequence.erase(block_sequence_iterator);
						}
						else
						{
							// else, we	will delete the rest of the block [until del_index|after del_index]--->>[until del_index|] //OI 16.3
							//adjusting the length of the block
							(*block_sequence_iterator).setLength(deletion_index);
							++block_sequence_iterator;
						}
						deletion_index = 0;//adjusting the relative index of deletion
					}
					else
					{
						// now, deletion_size < (*block_sequence_iterator).getLength() - deletion_index; 
						//thus we need to seperate the block to 2 different parts- before deletion and after deletion
						//[before delteion|deletion area|after deletion]--->> [before deletion],[after deletion]
						// OI 16.3
						if ((*block_sequence_iterator).getIndex() == -1)// extreme case- the deletion is in insertion block //added 24.3 OI
						{
							(*block_sequence_iterator).setLength((*block_sequence_iterator).getLength() - deletion_size);
							//indel_nums += deletion_size;//for debuging
							deletion_size = 0;


						}
						else
						{

							if (deletion_index == 0)
							{
								//extreme condition- if the deletion is in the start
								//here, we will change the start_index of the current block.// OI 16.3
								(*block_sequence_iterator).setIndex((*block_sequence_iterator).getIndex() + deletion_size);
								(*block_sequence_iterator).setLength((*block_sequence_iterator).getLength() - deletion_size);//addition 24.3 , we need to change the length of the block //OI
								deletion_size = 0;//adjusting the size of deletion to be 0 //OI 16.3
							}
							else
							{
								//avarege case- if the deletion is in the middle of the block
								//here, we will create another block and push him after the current block //OI 16.3
								// creation of the new block //OI 16.3
								Block insert_after_deletion((*block_sequence_iterator).getIndex() + deletion_index + deletion_size, (*block_sequence_iterator).getLength() - deletion_index - deletion_size);
								(*block_sequence_iterator).setLength(deletion_index);
								if (next(block_sequence_iterator) == block_sequence.end())
								{
									block_sequence.push_back(insert_after_deletion);
								}
								else
								{
									vector<Block>::iterator after_current_block = next(block_sequence_iterator);
									block_sequence.insert(after_current_block, insert_after_deletion);
								}
								deletion_size = 0;//adjusting the size of deletion to be 0 //OI 16.3
							}
						}


					}
				}
			}
		}// end of deletion //16.3 OI

		//printing temp //OI 29.3
		//printblockvec(block_sequence);
		//updating branch length and the other parameters. 
		branchlength = branchlength - waitingTime;
		sequenceWiseInsertionRate = 1.0 * _IR * (sequence_length + 1);
		sequenceWiseDeletionRate = 1.0 * _DR * sequence_length;
		waitingTime = drawExp(sequenceWiseInsertionRate + sequenceWiseDeletionRate);
	}
	//printing temp //OI 29.3
	//printblockvec(block_sequence);
	//now we finished the first part, and we are moving to the second... //OI 16.3

	//vector<Node*> output_sequence;// this is our output
	for (size_t i = 0; i < block_sequence.size(); i++)
	{
		Block current_block = block_sequence[i];
		if (current_block.getIndex() == -1)
		{
			//in this part, we will create the part of the new sequence, that was made by **insertion blocks** //OI 16.3
			//first, we need to create new LL
			LinkedList inserstion_LL_temp;//initilizing the temp linkedlist of insertion
			//now, before adding the node to the new sequence, we need to understand if it is empty //OI 17.3
			bool empty_flag = output_sequence.empty();

			//addition, in order to let me insert all the node to the sequence before I update the super sequence// OI 18.3 
			Node* insertafter = NULL;
			if (!empty_flag)
			{
				insertafter = output_sequence.back();
			}
			//adding the first node to our LinkedList
			Node* prev = new Node(currIdToInsert);
			inserstion_LL_temp.addnode_to_list(prev, NULL);
			//adding the first node of insertion to the output vector //OI 16.3
			output_sequence.push_back(prev);

			for (size_t i = 1; i < current_block.getLength(); i++)
			{
				// now we inserting the rest of the insertion blocks to the temp LinkedList, and inserting them to the vector
				Node* temp = new Node(currIdToInsert + i);
				inserstion_LL_temp.addnode_to_list(temp, prev);
				//adding rest of the nodes of insertion to the output vector //OI 16.3
				output_sequence.push_back(temp);
				prev = temp;//updating the node we want to be our "insert after" node (look at LinkedList method-> "addnode_to_list")
			}
			//updating the current Id 
			currIdToInsert += current_block.getLength();
			if (empty_flag)// small change- we now use the boolean we initilized earlier //OI 17.3
			{
				//this part is dealing with the case when the insertion block is at the start. //Oi 16.3
				superSequence.add_list_to_list(inserstion_LL_temp, NULL);
			}
			else
			{
				//regular case- when the insertion block **is not** at the start.
				superSequence.add_list_to_list(inserstion_LL_temp, insertafter);
			}
		}//ended the case of insertion block // OI 16.3
		else
		{
			//here, it is an original block //OI 16.3
			for (int i = current_block.getIndex(); i < current_block.getIndex() + current_block.getLength(); i++)
			{
				//In this part, we will create the part of the new sequence, that was made by original blocks(blocks from the ancestral sequence) //OI 16.3
				output_sequence.push_back(ancestralsequence[i]);//here, we pushing in the node pointer of the rightposition.
			}
		}//ended the original block case //16.3 OI
	}//ended the second half- updating the Super Sequence and creating the output sequence //OI 16.3
	//output_sequence.shrink_to_fit();
	return true;
}

//alternative way to initilize the MSA
void Simulator::simulationToMSA_new_with_changes(const vector<vector<Node*> >& simulatedLeavesSequences, vector<string>& msa) {
	vector<size_t> counter_SS(superSequence.getLength(), -1);
	for (size_t i = 0; i < simulatedLeavesSequences.size(); i++)
	{
		for (size_t j = 0; j < simulatedLeavesSequences[i].size(); j++)
		{
			if (counter_SS[(*simulatedLeavesSequences[i][j]).getData()] == -1)
			{
				counter_SS[(*simulatedLeavesSequences[i][j]).getData()] += 2;
			}
			else
			{
				counter_SS[(*simulatedLeavesSequences[i][j]).getData()]++;
			}

		}
	}

	for (size_t l = 0; l < simulatedLeavesSequences.size(); ++l) {//go over all the leaves
		string tmp;												  // reset string for the current sequence
		size_t j = 0;											  // reset index for the current sequence
		Node* pos = superSequence.getHead();
		for (size_t i = 0; i < superSequence.getLength(); ++i) {		  // go over the superSequence
			if (counter_SS[(*pos).getData()] != -1)
			{
				if (j == simulatedLeavesSequences[l].size() || counter_SS[(*pos).getData()] == 0) {		  // if currenr sequence ended, add "-"
					tmp.append("-");
				}
				else if ((*pos).getData() == (*simulatedLeavesSequences[l][j]).getData()) { // if current sequence = superSequence
					tmp.append("A");											// add arbitrary character
					++j;
					counter_SS[(*pos).getData()]--;
				}
				else {														  // else add "-"
					tmp.append("-");
				}
			}
			pos = (*pos).getNext();
		}
		msa.push_back(tmp);												// add current string to MSA
	}
}

//small changes from Dana's method(in the super sequnce name) //OI 17.3
// changed it to work with the LL class ... //OI 5.4
void Simulator::simulationToMSA_with_SuperSequence(const vector<vector<Node*> >& simulatedLeavesSequences, vector<string>& msa) {
	vector<size_t> counter_SS(superSequence.getLength());
	for (size_t i = 0; i < simulatedLeavesSequences.size(); i++)
	{
		for (size_t j = 0; j < simulatedLeavesSequences[i].size(); j++)
		{

			counter_SS[(*simulatedLeavesSequences[i][j]).getData()]++;

		}

	}

	for (size_t i = 0; i < counter_SS.size(); i++)
	{
		if (counter_SS[i] == 0) { counter_SS[i] = -1; }
	}
	for (size_t l = 0; l < simulatedLeavesSequences.size(); ++l) {//go over all the leaves
		string tmp;												  // reset string for the current sequence
		size_t j = 0;											  // reset index for the current sequence
		Node* pos = superSequence.getHead();
		for (size_t i = 0; i < superSequence.getLength(); ++i) {		  // go over the superSequence
			if (counter_SS[(*pos).getData()] != -1)
			{
				if (j == simulatedLeavesSequences[l].size() || counter_SS[(*pos).getData()] == 0) {		  // if currenr sequence ended, add "-"
					tmp.append("-");
				}
				else if ((*pos).getData() == (*simulatedLeavesSequences[l][j]).getData()) { // if current sequence = superSequence
					tmp.append("A");											// add arbitrary character
					++j;
					counter_SS[(*pos).getData()]--;
				}
				else {														  // else add "-"
					tmp.append("-");
				}
			}
			pos = (*pos).getNext();
		}
		msa.push_back(tmp);												// add current string to MSA
	}
}


//printing made by Dana.
void Simulator::outputInFastaFormat(const vector<Node*> & leafSequence) {
	size_t j = 0;
	Node* pos = superSequence.getHead();
	for (size_t i = 0; i < superSequence.getLength(); ++i) {
		if (j == leafSequence.size()) cout << "- ";
		else if (pos->getData() == leafSequence[j]->getData()) {
			cout << pos->getData()<<" ";
			++j;
		}
		else {
			cout << "- ";
		}
		pos = pos->getNext();
	}
}

void Simulator::outputInFastaFormat(const vector<vector<Node*> > & simulatedLeavesSequences) {
	for (size_t i = 0; i < simulatedLeavesSequences.size(); ++i) {
		cout << ">S" << i << endl;
		outputInFastaFormat(simulatedLeavesSequences[i]);
		cout << endl;
	}
}




//MY CODE //OI 15.3
bool Simulator::create_event_vector(double branchlength, int ancestralsequencelength, vector<Event>& event_vec)
{


	/*
	started and finished 15.3 by OI
	creating the event vector of each brnach, 15.3 OI

	here, we will run through the brnach, and will create the vector of events.
	logicly, it will be FiFo.
	*/
	//initialize starting parameters.

	// temp bool- to check "boomers"
	double sequenceWiseInsertionRate = 1.0 * _IR * (ancestralsequencelength + 1);
	double sequenceWiseDeletionRate = 1.0 * _DR * ancestralsequencelength;
	double waitingTime = drawExp(sequenceWiseInsertionRate + sequenceWiseDeletionRate);
	while (waitingTime < branchlength)
	{//running over the branch
		if (ancestralsequencelength > SpartaABC_options::_maxRLVal * 25)
		{
			//SpartaABC_options::boomed_times += 1;
			////cout << SpartaABC_options::boomed_times << endl;
			//is_boomed = false;
			return false;
		}

		bool isDeletion = true;
		if (uniform() < sequenceWiseInsertionRate / (sequenceWiseInsertionRate + sequenceWiseDeletionRate)) isDeletion = false;
		if (isDeletion) {

			//Deletion
			//initializing the start point and length of deletion //15.3 OI
			int theStartingPoint = uniform(0, ancestralsequencelength - 1);
			size_t deletionLength = _fastZDeletions.drawZip();
			if (theStartingPoint == 0 && deletionLength >= ancestralsequencelength)
			{//extreme cond- when we try to delete the entire sequence //OI 15.3
				deletionLength = ancestralsequencelength - 1;
			}
			if (theStartingPoint + deletionLength > ancestralsequencelength)
			{// extreme cond- when we deletion length passes the length of the sequence //OI 15.3
				deletionLength = ancestralsequencelength - theStartingPoint;
			}
			// inserting the event to the vevent vector and updating the length of the sequence.
			Event del(false, theStartingPoint, deletionLength);
			event_vec.push_back(del);
			ancestralsequencelength -= deletionLength;
		}
		else
		{//insertion
			//initializing the start point and length of insertion// 15.3 OI
			size_t theStartingPoint = uniform(0, ancestralsequencelength + 1);
			size_t insertionLength = _fastZInsertions.drawZip();
			// inserting the event to the vevent vector and updating the length of the sequence//OI  15.3
			Event ins(true, theStartingPoint, insertionLength);
			event_vec.push_back(ins);
			ancestralsequencelength += insertionLength;
		}
		//updating branch length and the other parameters. 
		branchlength = branchlength - waitingTime;
		sequenceWiseInsertionRate = 1.0 * _IR * (ancestralsequencelength + 1);
		sequenceWiseDeletionRate = 1.0 * _DR * ancestralsequencelength;
		waitingTime = drawExp(sequenceWiseInsertionRate + sequenceWiseDeletionRate);
	}
	return true;
}


//new funcs- genereated 15.3 //OI
void Simulator::simualteWithIndelsAlongAspecificBranch2(const vector<Event>& eventvec, const vector<Node*>& ancestralsequence,vector<Node*>& output_sequence)
{
	/*
	//strarted 15.3 OI

	simulating the events of insertions/ deletions along a specific branch in the tree.
	In the first half, we literly simulating the branch with our block vector.
	For example: ([0,100])->([0,50],[-1,5],[51,100])->([0,40],[45,50],[-1,5],[51,100]) is insertion of (51,5) and then deletion of (41,5)

	In the second half, we trnasfer the form from the block version to the node version, and we update the super sequence
	*/
	//initilizing the first block- with the size of the sequence.... //OI 15.3
	Block allsequence(0, ancestralsequence.size());
	vector<Block> block_sequence;
	block_sequence.push_back(allsequence);
	
	for (Event current_event:eventvec)
	{//running through all the events
		//cout << current_event.getIndel_Flag() << " , " << current_event.getIndex() << " , " << current_event.getLength() << endl;//in order to print events... OI 24.3	
		if (current_event.getIndel_Flag())
		{//Flag==TRUE---->> insertion //OI 15.3
			//initilizing the insert params //OI 15.3
			int insertion_index = current_event.getIndex();
			size_t insertion_length = current_event.getLength();
			//initilizing an iterator to run all over the block //comment added OI 16.3
			vector<Block>::iterator block_sequence_iterator = block_sequence.begin();
			while ( block_sequence_iterator != block_sequence.end()&&insertion_index>=(*block_sequence_iterator).getLength())
			{// iterating through the block sequence in order to find the right spot to insert //OI 15.3
				insertion_index -= (*block_sequence_iterator).getLength();
				++block_sequence_iterator;

			}
			Block insertion_block(-1, insertion_length);//generate new block- an insertion bloc //OI 15.3
			if (block_sequence_iterator == block_sequence.end())
			{
				//extreme condition- where the insertion is at the end //OI 15.3
				block_sequence.push_back(insertion_block);
			}
			else
			{ 
				if ((*block_sequence_iterator).getIndex() == -1)
				{
				// if it an insertion block- we just need to append the length of the block //OI 15.3
					(*block_sequence_iterator).setLength((*block_sequence_iterator).getLength() + insertion_length);
				}
				else
				{
					if (insertion_index == 0)
					{
						//easy condition- where we dont need to create to blocks
						block_sequence.insert(block_sequence_iterator, insertion_block);
					}
					else
					{
						// here, we generate the block we should create **after** the insertion block //OI 15.3
						Block after_insertion((*block_sequence_iterator).getIndex() + insertion_index, (*block_sequence_iterator).getLength() - insertion_index);
						//then we change the length of the original block //OI 15.3
						(*block_sequence_iterator).setLength(insertion_index);
						//now we have to different states
						if ((block_sequence_iterator != block_sequence.end()) && (next(block_sequence_iterator) == block_sequence.end()))// change of //16.3-->>18.3 // i did a wrong condition //OI
						{
							//1. this is the last bloc in the vector
							//thus we will **push back** the blocks //OI 15.3
							block_sequence.push_back(insertion_block);
							block_sequence.push_back(after_insertion);
						}
						else
						{
							//2. the block isnt the last block 
							// thus we will insert after the "block sequence iterator" //OI 15.3
							vector<Block>::iterator after_current_iterator = next(block_sequence_iterator);
							//note: change:
							//now we insert the after insertion block first, and then we point the iterator on him, and just after that we insert the insertion block... //OI 18.3
							after_current_iterator = block_sequence.insert(after_current_iterator, after_insertion);
							after_current_iterator = block_sequence.insert(after_current_iterator, insertion_block);


						}
					}
				}
			}
		}//ending insertion //OI 16.3
		else
		{//Flag==FALSE---->> deletion //OI 16.3
			//initilizing the deletion params ///OI 16.3
			size_t deletion_size = current_event.getLength();
			int deletion_index=current_event.getIndex();
			//initilizing an iterator to run all over the block //OI 16.3
			vector<Block>::iterator block_sequence_iterator = block_sequence.begin();
			while (deletion_size > 0 && block_sequence_iterator!=block_sequence.end())
			{
				//now we are running on the block sequence and deketing the necessary //OI 16.3
				if ((*block_sequence_iterator).getLength() <= deletion_index)
				{
					//first step- if the current block isn't the block we need to delete //OI 16.3
					deletion_index -= (*block_sequence_iterator).getLength();//adjusting the index of the deletion
					++block_sequence_iterator;
				}
				// from here we start to delete //OI 16.3
				else
				{
					if (deletion_size >= (*block_sequence_iterator).getLength() - deletion_index)
					{

						//temp for debuging OI //28.3
						if ((*block_sequence_iterator).getIndex() == -1) {
							//indel_nums += (*block_sequence_iterator).getLength() - deletion_index;// for debuging
						}
						//first condition- if the deletion length is longer then the length of the block //OI 16.3
						deletion_size = deletion_size - (*block_sequence_iterator).getLength() + deletion_index; //adjusting the deletion size
						if (deletion_index == 0)
						{
							//extreme condition- if the deletion starts fron the 0-index (we need to delete the block)// OI 16.3
							
							block_sequence_iterator = block_sequence.erase(block_sequence_iterator);
						}
						else
						{
							// else, we	will delete the rest of the block [until del_index|after del_index]--->>[until del_index|] //OI 16.3
							//adjusting the length of the block
							(*block_sequence_iterator).setLength(deletion_index);
							++block_sequence_iterator;
						}
						deletion_index = 0;//adjusting the relative index of deletion
					}
					else
					{
						// now, deletion_size < (*block_sequence_iterator).getLength() - deletion_index; 
						//thus we need to seperate the block to 2 different parts- before deletion and after deletion
						//[before delteion|deletion area|after deletion]--->> [before deletion],[after deletion]
						// OI 16.3
						if ((*block_sequence_iterator).getIndex()==-1)// extreme case- the deletion is in insertion block //added 24.3 OI
						{
							(*block_sequence_iterator).setLength((*block_sequence_iterator).getLength() - deletion_size);
							//indel_nums += deletion_size;//for debuging
							deletion_size = 0;
							
						
						}
						else
						{

							if (deletion_index == 0)
							{
								//extreme condition- if the deletion is in the start
								//here, we will change the start_index of the current block.// OI 16.3
								(*block_sequence_iterator).setIndex((*block_sequence_iterator).getIndex() + deletion_size);
								(*block_sequence_iterator).setLength((*block_sequence_iterator).getLength() - deletion_size);//addition 24.3 , we need to change the length of the block //OI
								deletion_size = 0;//adjusting the size of deletion to be 0 //OI 16.3
							}
							else
							{
								//avarege case- if the deletion is in the middle of the block
								//here, we will create another block and push him after the current block //OI 16.3
								// creation of the new block //OI 16.3
								Block insert_after_deletion((*block_sequence_iterator).getIndex() + deletion_index + deletion_size, (*block_sequence_iterator).getLength() - deletion_index - deletion_size);
								(*block_sequence_iterator).setLength(deletion_index);
								if (next(block_sequence_iterator) == block_sequence.end() )
								{
									block_sequence.push_back(insert_after_deletion);
								}
								else
								{
									vector<Block>::iterator after_current_block = next(block_sequence_iterator);
									block_sequence.insert(after_current_block, insert_after_deletion);
								}
								deletion_size = 0;//adjusting the size of deletion to be 0 //OI 16.3
							}
						}
						

					}
				}					
			}
		}// end of deletion //16.3 OI

		//printing temp //OI 29.3
		//printblockvec(block_sequence);
	}
	//printing temp //OI 29.3
	//printblockvec(block_sequence);
	//now we finished the first part, and we are moving to the second... //OI 16.3
	
	//vector<Node*> output_sequence;// this is our output
	for (size_t i=0;i<block_sequence.size();i++)
	{
		Block current_block = block_sequence[i];
		if (current_block.getIndex() == -1)
		{
			//in this part, we will create the part of the new sequence, that was made by **insertion blocks** //OI 16.3
			//first, we need to create new LL
			LinkedList inserstion_LL_temp;//initilizing the temp linkedlist of insertion
			//now, before adding the node to the new sequence, we need to understand if it is empty //OI 17.3
			bool empty_flag = output_sequence.empty();

			//addition, in order to let me insert all the node to the sequence before I update the super sequence// OI 18.3 
			Node* insertafter = NULL;
			if (!empty_flag)
			{
				insertafter=output_sequence.back();
			}
			//adding the first node to our LinkedList
			Node* prev=new Node(currIdToInsert);
			inserstion_LL_temp.addnode_to_list(prev, NULL);
			//adding the first node of insertion to the output vector //OI 16.3
			output_sequence.push_back(prev);
			
			for (size_t i = 1; i < current_block.getLength(); i++)
			{
				// now we inserting the rest of the insertion blocks to the temp LinkedList, and inserting them to the vector
				Node* temp = new Node(currIdToInsert + i);
				inserstion_LL_temp.addnode_to_list(temp, prev);
				//adding rest of the nodes of insertion to the output vector //OI 16.3
				output_sequence.push_back(temp);
				prev = temp;//updating the node we want to be our "insert after" node (look at LinkedList method-> "addnode_to_list")
			}
			//updating the current Id 
			currIdToInsert += current_block.getLength();
			if (empty_flag)// small change- we now use the boolean we initilized earlier //OI 17.3
			{
				//this part is dealing with the case when the insertion block is at the start. //Oi 16.3
				superSequence.add_list_to_list(inserstion_LL_temp,NULL);
			}
			else
			{
				//regular case- when the insertion block **is not** at the start.
				superSequence.add_list_to_list(inserstion_LL_temp, insertafter);
			}
		}//ended the case of insertion block // OI 16.3
		else
		{
			//here, it is an original block //OI 16.3
			for (int i = current_block.getIndex(); i < current_block.getIndex() + current_block.getLength(); i++)
			{
				//In this part, we will create the part of the new sequence, that was made by original blocks(blocks from the ancestral sequence) //OI 16.3
				output_sequence.push_back(ancestralsequence[i]);//here, we pushing in the node pointer of the rightposition.
			}
		}//ended the original block case //16.3 OI
	}//ended the second half- updating the Super Sequence and creating the output sequence //OI 16.3
	//output_sequence.shrink_to_fit();
	return;
}

/***
	OLD STUFF
			***/

//new generation of sequence //OI 17.3
vector<Node*> Simulator::generateRootSequence()
{
	//int this method, we will generate a vector<Node*> ancestral sequence, while we initilize our Super Sequnce...
	//OI 17.3
	vector<Node*> ancestral_sequnce;
	//creating the first Node of the sequnce
	Node* prev = new Node(0);
	ancestral_sequnce.push_back(prev);
	//creating temp Super Sequnce	//OI 17.3
	LinkedList superSequence_temp;
	superSequence_temp.addnode_to_list(prev, NULL);
	for (size_t i = 1; i < _rootLength; i++)
	{
		//now we will run and create all our nodes....
		Node* pos = new Node(i);
		ancestral_sequnce.push_back(pos);
		//inserting the nodes to the supersequnce //OI 17.3
		superSequence_temp.addnode_to_list(pos, prev);
		prev = pos;
	}
	//updating the **real** super sequence //OI 17.3
	superSequence = superSequence_temp;
	return ancestral_sequnce;

}

//changing the vectors of simulated leaves from nodes to numbers
void Simulator::back_to_numbers(const vector<vector<Node*>>& upgraded_vector, vector<vector<nucID>>& numbered_vector)
{
	//again. the method is seperated to 2 parts
	//the first- changing the simulated vectors
	//the second- chaning the upgraded super sequnce
	//OI 17.3

	for (vector<Node*> simulated_vector : upgraded_vector)
	{//first part- as I  mentioned above
		vector<nucID>* temp = new vector<nucID>;
		numbered_vector.push_back(*temp);
		for (size_t i = 0; i < simulated_vector.size(); i++)
		{
			numbered_vector.back().push_back((*simulated_vector[i]).getData());//tranforming the simulated leaves --->>> to the numbered vectors...
		}
	}

	//the second part as I mentioned above
	vector<nucID> temp_super_sequence;
	Node* pos = superSequence.getHead();
	while (pos != NULL)
	{
		//here, we run through the linked list ...
		temp_super_sequence.push_back(pos->getData());// transforming the super sequence to a numbered one. 
		pos = pos->getNext();
	}
	numbered_superSequnce = temp_super_sequence;//assigning the numbered_upgraded super sequence //OI 17.3
}


/*****
	Even Older Code
					*****/


//Dana's Simulator
//void Simulator::simualteWithIndelsAlongAspecificBranch(
//	const vector<nucID>& ancestralSequence,
//	double branchLength,
//	vector<nucID>& outputSeq){ // the output
//	//cout << "entering simulate along a branch" << endl;
//	// The son sequence is initiated to be the father, and then, later, indels are introduced
//	outputSeq = ancestralSequence;
//
//	// The waiting time for an event depends on the insertion rate and the deletion rate
//	// The waiting time is for the entire sequence and not for a specific position
//	// The _IR param is for position, so to get the sequence-leven insertion rate
//	// one has to multiply the _IR parameter with the current length of the sequence.
//	// Same for the DR.
//	// However, assume a sequence of size 3 xyz
//	// deletions can start before x, before y and before z -> 3 places.
//	// insertions can start before x, before y, before z or after z -> 4 places
//
//	double sequenceWiseInsertionRate = 1.0*_IR*(ancestralSequence.size() + 1);
//	double sequenceWiseDeletionRate = 1.0*_DR *ancestralSequence.size();
//	double waitingTime = drawExp(sequenceWiseInsertionRate + sequenceWiseDeletionRate);
//	if (waitingTime > branchLength) {
//		return;
//	}
//	// generate one events, i.e., an insertion or a deletion
//	bool isDeletion = true;
//	
//	if (uniform() < sequenceWiseInsertionRate / (sequenceWiseInsertionRate + sequenceWiseDeletionRate)) isDeletion = false;
//	if (isDeletion) {
//		//cout << "entering isDeletion ";
//		int theStartingPoint = uniform(0, ancestralSequence.size() - 1);
//		size_t deletionLength = _fastZDeletions.drawZip();
//		if (theStartingPoint == 0 && deletionLength >= ancestralSequence.size()) {
//			deletionLength = ancestralSequence.size()-1; // we do not want the entire sequence to be deleted...
//		}
//	//	cout << "deletion length=" << deletionLength<< endl;
//		//int deletionLength = powerLaw(powerLawParam, N);
//		size_t positionInWhichDeletionEnds = theStartingPoint + deletionLength;
//		if (positionInWhichDeletionEnds > outputSeq.size()) positionInWhichDeletionEnds = outputSeq.size();
//		outputSeq.erase(outputSeq.begin() + theStartingPoint, outputSeq.begin() + positionInWhichDeletionEnds);
//		simualteWithIndelsAlongAspecificBranch(outputSeq, branchLength - waitingTime, outputSeq);
//	//	cout << "exiting deletion" << endl;
//		return;
//	}
//	else { // insertions
//		//cout << "entering insertion ";
//		size_t theStartingPoint = uniform(0, ancestralSequence.size()+1);
//		//size_t insertionLength = drawZipf(_A_Insertion_param, MaxIndelSize);
//		size_t insertionLength = _fastZInsertions.drawZip();
//		//cout << "insertion length=" << insertionLength << endl;
//		//int insertionLength = powerLaw(powerLawParam, N);
//		vector<nucID>::iterator it;
//		if (theStartingPoint > ancestralSequence.size()) {
//			it = outputSeq.end();
//		}
//		else {
//			it = outputSeq.begin()+theStartingPoint;
//		}
//		//Initiating the array of new IDs to insert
//		for (size_t m = 0; m < insertionLength; m++) {
//			tmparray[m] = currIdToInsert;
//			currIdToInsert++;
//		}
//
//		// this part updates the super sequence
//		int numberBeforeInsertion = -1; // this indicates an insertion before 0.
//		if (theStartingPoint > outputSeq.size()) {
//			numberBeforeInsertion = outputSeq[outputSeq.size() - 1];
//		}
//		else {
//			if (theStartingPoint >= 1) numberBeforeInsertion = outputSeq[theStartingPoint - 1];
//		}
//		int j = 0; 
//		if (numberBeforeInsertion == -1) j = -1;
//		else {
//			for (; j < superSequence.size(); ++j) {
//				if (superSequence[j] == numberBeforeInsertion) break;
//			}
//		}
//		vector<nucID>::iterator it2 = superSequence.begin()+(j +1);
//		superSequence.insert(it2, tmparray, tmparray+ insertionLength);
//
//		// updating the sequence itself
//		//it2 = it2 + insertionLength;
//		outputSeq.insert(it, tmparray, tmparray + insertionLength);
//
//		simualteWithIndelsAlongAspecificBranch(outputSeq, branchLength - waitingTime, outputSeq);
//
//	}
//}
//For testing validity? Im not sure...
//void Simulator::simualteWithIndelsAlongAspecificBranch_with_events(
//	const vector<nucID>& ancestralSequence,
//	double branchLength,
//	vector<nucID>& outputSeq,vector<Event>& event_vec) { // the output
//	// WE DONT USE THIS METHOD ANYMORE!! 
//	//Addition-29.3- we change to be iterative look at "simualteWithIndelsAlongAspecificBranch_with_events_not_rec" the itertiave way.
//	// The recursive way causes stack overflow
//
//	// this method is used in order to merge the events of Dana's code with the new code events....
//
//
//	//cout << "entering simulate along a branch" << endl;
//	// The son sequence is initiated to be the father, and then, later, indels are introduced
//	outputSeq = ancestralSequence;
//	
//	// The waiting time for an event depends on the insertion rate and the deletion rate
//	// The waiting time is for the entire sequence and not for a specific position
//	// The _IR param is for position, so to get the sequence-leven insertion rate
//	// one has to multiply the _IR parameter with the current length of the sequence.
//	// Same for the DR.
//	// However, assume a sequence of size 3 xyz
//	// deletions can start before x, before y and before z -> 3 places.
//	// insertions can start before x, before y, before z or after z -> 4 places
//
//	double sequenceWiseInsertionRate = 1.0 * _IR * (ancestralSequence.size() + 1);
//	double sequenceWiseDeletionRate = 1.0 * _DR * ancestralSequence.size();
//	double waitingTime = drawExp(sequenceWiseInsertionRate + sequenceWiseDeletionRate);
//	if (waitingTime > branchLength) {
//		// the to next lines are our stop condition.... //OI 17.3
//		//cout << "finished" << endl;
//		return;
//	}
//	// generate one events, i.e., an insertion or a deletion
//	bool isDeletion = true;
//	if (uniform() < sequenceWiseInsertionRate / (sequenceWiseInsertionRate + sequenceWiseDeletionRate)) isDeletion = false;
//	if (isDeletion) {
//		//cout << "entering isDeletion ";
//		int theStartingPoint = uniform(0, ancestralSequence.size() - 1);
//		size_t deletionLength = _fastZDeletions.drawZip();
//		if (theStartingPoint == 0 && deletionLength >= ancestralSequence.size()) {
//			deletionLength = ancestralSequence.size() - 1; // we do not want the entire sequence to be deleted...
//		}
//		//	cout << "deletion length=" << deletionLength<< endl;
//			//int deletionLength = powerLaw(powerLawParam, N);
//		size_t positionInWhichDeletionEnds = theStartingPoint + deletionLength;
//		if (positionInWhichDeletionEnds > outputSeq.size()) positionInWhichDeletionEnds = outputSeq.size();
//		outputSeq.erase(outputSeq.begin() + theStartingPoint, outputSeq.begin() + positionInWhichDeletionEnds);
//		//cout <<"old: "<< false << ", " << theStartingPoint << ", " << deletionLength << endl;// added 24.3 OI,   printed in order to check equality
//		Event deletion_event(false, theStartingPoint, deletionLength);
//		event_vec.push_back(deletion_event);
//		simualteWithIndelsAlongAspecificBranch_with_events(outputSeq, branchLength - waitingTime, outputSeq,event_vec);//small change- added event vec
//		//**THE INTRESTING PART OF THE CHANGE!!!!**
//		//here, we add the event of the indel (in this time - deletion)
//		//and we add it in the start : [in, del, del ,in]--->>>[del,in,del,del,in]
//		//OI 17.3
//		
//		//	cout << "exiting deletion" << endl;
//		return ;
//	}
//	else { // insertions
//		//cout << "entering insertion ";
//		size_t theStartingPoint = uniform(0, ancestralSequence.size() + 1);
//		//size_t insertionLength = drawZipf(_A_Insertion_param, MaxIndelSize);
//		size_t insertionLength = _fastZInsertions.drawZip();
//		
//		//cout << "insertion length=" << insertionLength << endl;
//		//int insertionLength = powerLaw(powerLawParam, N);
//		vector<nucID>::iterator it;
//		if (theStartingPoint > ancestralSequence.size()) {
//			it = outputSeq.end();
//		}
//		else {
//			it = outputSeq.begin() + theStartingPoint;
//		}
//		//Initiating the array of new IDs to insert
//		for (size_t m = 0; m < insertionLength; m++) {
//			tmparray[m] = currIdToInsert;
//			currIdToInsert++;
//		}
//
//		// this part updates the super sequence
//		int numberBeforeInsertion = -1; // this indicates an insertion before 0.
//		if (theStartingPoint > outputSeq.size()) {
//			numberBeforeInsertion = outputSeq[outputSeq.size() - 1];
//		}
//		else {
//			if (theStartingPoint >= 1) numberBeforeInsertion = outputSeq[theStartingPoint - 1];
//		}
//		int j = 0;
//		if (numberBeforeInsertion == -1) j = -1;
//		else {
//			for (; j < superSequence.size(); ++j) {
//				if (superSequence[j] == numberBeforeInsertion) break;
//			}
//		}
//		vector<nucID>::iterator it2 = superSequence.begin() + (j + 1);
//		superSequence.insert(it2, tmparray, tmparray + insertionLength);
//
//		// updating the sequence itself
//		//it2 = it2 + insertionLength;
//		outputSeq.insert(it, tmparray, tmparray + insertionLength);
//		//from here - new! //OI 17.3
//		//**THE INTRESTING PART OF THE CHANGE!!!!**
//		//here, we add the event of the indel (in this time - insertion)
//		//and we add it in the start : [in, del, del ,in]--->>>[in,in,del,del,in]
//		//OI 17.3
//		//cout <<"old: "<< true << ", " << theStartingPoint << "," << insertionLength << "," << endl;// added 24.3 OI,   printed in order to check equality
//		Event insertion_event(true, theStartingPoint, insertionLength);
//		event_vec.push_back(insertion_event);
//		simualteWithIndelsAlongAspecificBranch_with_events(outputSeq, branchLength - waitingTime, outputSeq,event_vec);
//		
//		return ;
//	}
//}


//Another old Simulator
//void Simulator::simualteWithIndelsAlongAspecificBranch_with_events_not_rec(
//	const vector<nucID>& ancestralSequence,
//	double branchLength,
//	vector<nucID>& outputSeq, vector<Event>& event_vec) { // the output
//
//	// this method is used in order to merge the events of Dana's code with the new code events....
//	//Addition-29.3- we change to be iterative.
//	// The recursive way causes stack overflow
//
//	//cout << "entering simulate along a branch" << endl;
//	// The son sequence is initiated to be the father, and then, later, indels are introduced
//	outputSeq = ancestralSequence;
//
//	// The waiting time for an event depends on the insertion rate and the deletion rate
//	// The waiting time is for the entire sequence and not for a specific position
//	// The _IR param is for position, so to get the sequence-leven insertion rate
//	// one has to multiply the _IR parameter with the current length of the sequence.
//	// Same for the DR.
//	// However, assume a sequence of size 3 xyz
//	// deletions can start before x, before y and before z -> 3 places.
//	// insertions can start before x, before y, before z or after z -> 4 places
//
//	double sequenceWiseInsertionRate = 1.0 * _IR * (outputSeq.size() + 1);
//	double sequenceWiseDeletionRate = 1.0 * _DR * outputSeq.size();
//	double waitingTime = drawExp(sequenceWiseInsertionRate + sequenceWiseDeletionRate);
//	while (waitingTime < branchLength)
//	{
//		// generate one events, i.e., an insertion or a deletion
//		bool isDeletion = true;
//		if (uniform() < sequenceWiseInsertionRate / (sequenceWiseInsertionRate + sequenceWiseDeletionRate)) isDeletion = false;
//		if (isDeletion) {
//			//cout << "entering isDeletion ";
//			int theStartingPoint = uniform(0, outputSeq.size() - 1);
//			size_t deletionLength = _fastZDeletions.drawZip();
//			if (theStartingPoint == 0 && deletionLength >= outputSeq.size()) {
//				deletionLength = outputSeq.size() - 1; // we do not want the entire sequence to be deleted...
//			}
//			//	cout << "deletion length=" << deletionLength<< endl;
//				//int deletionLength = powerLaw(powerLawParam, N);
//			size_t positionInWhichDeletionEnds = theStartingPoint + deletionLength;
//			if (positionInWhichDeletionEnds > outputSeq.size()) positionInWhichDeletionEnds = outputSeq.size();
//			outputSeq.erase(outputSeq.begin() + theStartingPoint, outputSeq.begin() + positionInWhichDeletionEnds);
//			//cout <<"old: "<< false << ", " << theStartingPoint << ", " << deletionLength << endl;// added 24.3 OI,   printed in order to check equality
//			Event deletion_event(false, theStartingPoint, deletionLength);
//			event_vec.push_back(deletion_event);
//			
//			//**THE INTRESTING PART OF THE CHANGE!!!!**
//			//here, we add the event of the indel (in this time - deletion)
//			//and we add it in the start : [in, del, del ,in]--->>>[del,in,del,del,in]
//			//OI 17.3
//
//			//	cout << "exiting deletion" << endl;
//		}
//		else { // insertions
//			//cout << "entering insertion ";
//			size_t theStartingPoint = uniform(0, outputSeq.size() + 1);
//			//size_t insertionLength = drawZipf(_A_Insertion_param, MaxIndelSize);
//			size_t insertionLength = _fastZInsertions.drawZip();
//
//			//cout << "insertion length=" << insertionLength << endl;
//			//int insertionLength = powerLaw(powerLawParam, N);
//			vector<nucID>::iterator it;
//			if (theStartingPoint > outputSeq.size()) {
//				it = outputSeq.end();
//			}
//			else {
//				it = outputSeq.begin() + theStartingPoint;
//			}
//			//Initiating the array of new IDs to insert
//			for (size_t m = 0; m < insertionLength; m++) {
//				tmparray[m] = currIdToInsert;
//				currIdToInsert++;
//			}
//
//			// this part updates the super sequence
//			int numberBeforeInsertion = -1; // this indicates an insertion before 0.
//			if (theStartingPoint > outputSeq.size()) {
//				numberBeforeInsertion = outputSeq[outputSeq.size() - 1];
//			}
//			else {
//				if (theStartingPoint >= 1) numberBeforeInsertion = outputSeq[theStartingPoint - 1];
//			}
//			int j = 0;
//			if (numberBeforeInsertion == -1) j = -1;
//			else {
//				for (; j < superSequence.size(); ++j) {
//					if (superSequence[j] == numberBeforeInsertion) break;
//				}
//			}
//			vector<nucID>::iterator it2 = superSequence.begin() + (j + 1);
//			superSequence.insert(it2, tmparray, tmparray + insertionLength);
//
//			// updating the sequence itself
//			//it2 = it2 + insertionLength;
//			outputSeq.insert(it, tmparray, tmparray + insertionLength);
//			//from here - new! //OI 17.3
//			//**THE INTRESTING PART OF THE CHANGE!!!!**
//			//here, we add the event of the indel (in this time - insertion)
//			//and we add it in the start : [in, del, del ,in]--->>>[in,in,del,del,in]
//			//OI 17.3
//			//cout <<"old: "<< true << ", " << theStartingPoint << "," << insertionLength << "," << endl;// added 24.3 OI,   printed in order to check equality
//			Event insertion_event(true, theStartingPoint, insertionLength);
//			event_vec.push_back(insertion_event);
//		}
//		//here we change the paramters in order to keep runing overthe branch //OI 29.3
//		branchLength = branchLength - waitingTime;
//		sequenceWiseInsertionRate = 1.0 * _IR * (outputSeq.size() + 1);
//		sequenceWiseDeletionRate = 1.0 * _DR * outputSeq.size();
//		waitingTime = drawExp(sequenceWiseInsertionRate + sequenceWiseDeletionRate);
//	}
//	return;
//}

//old MSA initilizer
//void simulationToMSA(const vector<vector<nucID> >& simulatedLeavesSequences, vector<string>& msa) {
//	vector<size_t> counter_SS(superSequence.size());
//	for (size_t i = 0; i < simulatedLeavesSequences.size(); i++)
//	{
//		for (size_t j = 0; j < simulatedLeavesSequences[i].size(); j++)
//		{
//			counter_SS[simulatedLeavesSequences[i][j]]++;
//		}
//
//	}
//	for (size_t i = 0; i < counter_SS.size(); i++)
//	{
//		if (counter_SS[i] == 0) { counter_SS[i] = -1; }
//	}
//	for (size_t l = 0; l < simulatedLeavesSequences.size(); ++l) {//go over all the leaves
//		string tmp;												  // reset string for the current sequence
//		size_t j = 0;											  // reset index for the current sequence
//		for (size_t i = 0; i < superSequence.size(); ++i) {		  // go over the superSequence
//			if (counter_SS[superSequence[i]] != -1)
//			{
//				if (j == simulatedLeavesSequences[l].size() || counter_SS[superSequence[i]] == 0) {		  // if currenr sequence ended, add "-"
//					tmp.append("-");
//				}
//				else if (superSequence[i] == simulatedLeavesSequences[l][j]) { // if current sequence = superSequence
//					tmp.append("A");											// add arbitrary character
//					++j;
//					counter_SS[superSequence[i]]--;
//				}
//				else {														  // else add "-"
//					tmp.append("-");
//				}
//			}
//		}
//		msa.push_back(tmp);												// add current string to MSA
//	}
//}



//vector<string> Simulator::mainSimualteAlongAtree()
//{
//	int length = 3000;
//	vector<vector<nucID> > simulatedLeavesSequences;
//	vector<nucID> ancestralSequence;
//	int currIdToInsert = length;
//	for (int i = 0; i < length; ++i) {
//		ancestralSequence.push_back(i);
//	}
//	superSequence = ancestralSequence;
//	string treeFile = "C:\\Users\\Dana Rapoport\\Downloads\\SpartaABC_with_INDELible_20170320_bundle\\SpartaABC_with_INDELible_20170320_bundle\\run_example\\small_tree.txt";
//	tree t(treeFile);
//	simulateAlongTree(t.getRoot(), ancestralSequence, simulatedLeavesSequences, 1.3, 0.02, 0.02);
//	vector<string> msa;
//	simulationToMSA(simulatedLeavesSequences, msa);
//	return msa;
//	//for (int k = 0; k < msa.size(); ++k) { cout << msa[k] << endl; }
//	/*printSequences(simulatedLeavesSequences);
//	// printing the super sequence:
//	cout << "\nthe super sequence is: " << endl;
//	for (size_t k = 0; k < superSequence.size(); ++k) {
//		cout << superSequence[k] << " ";
//	}
//	cout << endl;
//	cout << " the final leaf sequences in Fasta format\n";
//	outputInFastaFormat(simulatedLeavesSequences);
//	*/
//	//return 0;
//}
