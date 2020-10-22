#ifndef ___ALPHABET_H
#define ___ALPHABET_H

#include <string>
#include <vector>
using namespace std;

// Explanations added by Tal Pupko on 1/6/2017.
//Alphabet is a pure virtual class. 
//It means that every time you want to use an alphabet, you must derive a specific alphabet subclass from this class.
//An example of an alphabet is that describing amino acids, and is defined in the file "amino.h". 
//So in this file, we set the "roles", what each class of alphabet should contain.
//Each alphabet must have a size. For example, the size of the amino alphabet is 20.
//The alphabet is actually a map between characters and integers. Each character can be converted to an integer and back.
//Thus, for example, in class amino, the amino acid "E" is coded as 6.
//The idea is that a sequence would then be stored as a vector of integers.
//Note that when a codon alphabet is read, each three characters will be converted to a single integer. 
//This also explains the motivation to use integers rather than characters.
//This explain the second function, "fromString". It takes an input a string and returns a vector of integers.
//The reverse of this function is the function "fromInt". It gets the integer and returns the string associated with it.
//In class amino, fromInt(6) will return E. In class codon, fromInt(7) will return a string of length 3, which is the codon encoded as the number 7.
//Each alphabet must also have an integer that represents a gap. In class amino for example, a gap is represented by the number -1.
//When reading a sequence, the character '-' is encoded as -1.
//The function gap() returns the integer associated with the gap character (in amino it will return -1).
//Similarly, some sequences contain the character "?" or "x". This usually indicates that this part of the sequence is unknown.
//This is also coded as an integer (in class amino, it is given the integer -2).
//The function stringSize() returns the number of characters that are encoded to each integer. 
// For codons it would be 3, and for nucleotides it will be 1.
//The function "from char(string, pos) gets as input a string (e.g., AAAGGG) and a position (e.g., 0) and returns the 
//integer that is associated with that position. So in case of codons, it will read 3 characters and will return the 
//integer associated with AAA. In case of nucleotide, it will read only the first "A" and return the integer associated with "A".
//In each alphabet some characters represent genuine characters and some represent gaps and ambigiuous characters.
//For example, ?, - and X are not real nucleotides. N is not real too, as it represents either A, C, G, or T.
//The function isSpecific will return true only if it gets as input A,C,G,T and will return false eitherwise.

// Explanation about relations:
// Sometimes the sequences contain letters like R which means, in the case of nucleodies, G or A.
// When calculating the likelihood of such sequences we have to take this into acount.
// For example, the tree :		                            A
//				 				                           / \
//													   t1 /   \ t2
//												    	 /     \
//													    R  	    A
//
//		L = P(A)*P(A->A)(t1)*P(A->A)(t2) + P(A)*P(A->G)(t1)*P(A->A)(t2)
//		= P(A)*P(A->A)(t2)* [ P(A->A)(t1) + P(A->G)(t1) ]
//
//		Note that we don't divide it by 2.
//		A table (VVint _relation) keeps this information :

//		A	C	G	T
//A		1	0	0	0
//C		0	1	0	0
//G		0	0	1	0
//T		0	0	0	1
//U		0	0	0	1
//R		1	0	1	0
//Y		0	1	0	1
// The function relations(const int charInSeq, const int charToCheck) will return
// _relation[charToCheck][charInSeq];
//In the above example, assume R is coded as 5 and A is coded as 0, relation (5,0) will return 1, because R includes A, and relation (5,1)
// will return 0, because R is either A or G, but not C (which is encoded as 1).




class alphabet {
public:
	virtual size_t size() const = 0;
	//size of the alphabet. e.g. the size of the amino alphabet is 20.
	virtual vector<size_t> fromString(const string& str) const = 0;
	//fromString- get a sequence string and transform it to an integer vector 
	//(performed in the relevant subclass)
	virtual string fromInt(const size_t in_id) const = 0;
	//fromString- get an integer vector that represents a sequence
	//and transform it to a sequence string 
	//(performed in the relevant subclass)
	virtual size_t gap() const = 0;
	//The function gap() returns the integer associated with the gap character(in amino it will return -1).
	virtual size_t unknown() const = 0;
	virtual size_t stringSize() const = 0;
	//The function stringSize() returns the number of characters that are encoded to each integer. 
	// For codons it would be 3, and for nucleotides it will be 1.
	virtual size_t fromChar(const string& seq, const size_t pos) const = 0;
	//The function "from char(string, pos) gets as input a string (e.g., AAAGGG) and a position (e.g., 0) and returns the 
	//integer that is associated with that position. So in case of codons, it will read 3 chars, and in case of nucleotides only 1 char
	virtual bool isSpecific(const size_t size_t) const = 0; // "specific" here is not unknown, nor ambiguity, nor gap (for example, for nucleotides it will true for A,C,G, or T).

	virtual size_t relations(const size_t charInSeq, const size_t charToCheck) const =0;
	
	
	
	virtual ~alphabet()=0;
	virtual alphabet* clone() const = 0;
	
};

#endif

