#include "read_seqs.h"

vector<string> read_fasta_from_file(string filename)
{
	vector<string> seqs;
	ifstream myFile;
	myFile.open(filename.c_str());	
	if(!myFile) 
	{
      cout<<"can not open file:" <<filename<<endl;
	  exit(1);
    }
	
	string name;
	string seq;
	string temp;
	
	getline(myFile,temp);
	while (! (myFile.peek() == EOF)) 
	{
		
		// line may be empty - ignore blank lines
		if (temp.empty())
		{
			getline(myFile,temp);
			continue;
		}
			
		if (temp[0] == '>') 
		{
			name = temp.substr(1);
			//remove white spaces from end of 'name':
			name.erase(name.find_last_not_of(" \n\r\t")+1);

			seq.clear();
			getline(myFile,temp);
			while (! (myFile.peek() == EOF))
			{
				//remove white spaces from end of 'temp':
				temp.erase(temp.find_last_not_of(" \n\r\t")+1);

				// line may be empty - ignore blank lines
				if (temp.empty())
				{
					getline(myFile,temp);
					continue;
				}

				// reached next header - push seq of previous to vector
				if(temp[0] == '>')
				{
					seqs.push_back(seq);
					seq.clear();
					break;
				}

				// concatenate sequence line to seq and read the next line
				seq += temp;
				getline(myFile,temp);
			}
			
		}

		// the last sequence in the file is not pushed in the inner while, so we push it here
		if ((myFile.peek() == EOF)&&((seq.length() > 0)||(temp.length()>0)))
		{
			temp.erase(temp.find_last_not_of(" \n\r\t") + 1);
			seq += temp;
			seqs.push_back(seq);
			seq.clear();
			break;
		}

	}

	myFile.close();
	return seqs;
}