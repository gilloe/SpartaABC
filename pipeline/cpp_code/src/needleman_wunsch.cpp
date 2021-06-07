#include "needleman_wunsch.h"

needleman_wunsch::~needleman_wunsch()
{
	_seqs.clear();
}

needleman_wunsch::needleman_wunsch(const needleman_wunsch& other )//copy constructor
{
	_seqs = other._seqs;
	_numberOfSequences = other._numberOfSequences;
	_match_score = other._match_score;
	_mismatch_score = other._mismatch_score;
	_gap_open = other._gap_open;
	_gap_extend = other._gap_extend;
	_simMat = other._simMat;
	_similarity_mode = other._similarity_mode;
}

void needleman_wunsch::computeAllPairsMSAs(vector<MSA> & pairwiseMsasToFill)
{
	for(int i = 0; i < _numberOfSequences; i++)
	{
		for(int j = i + 1; j < _numberOfSequences; j++)
		{
			MSA curr_pairwise = computeAffineNWForPair(_seqs[i],_seqs[j]); //linear can be achieved if gap_open=gap_extend
			pairwiseMsasToFill.push_back(curr_pairwise);
		}
	}
}


MSA needleman_wunsch::computeAffineNWForPair(const string A, const string B)
{
	int m = A.length();
	int n = B.length();

	// deal with edge cases where one or more of the sequences is of size zero:
	if((m == 0) && (n == 0))
	{
		vector<string> aligned_seqs;
		aligned_seqs.push_back("");
		aligned_seqs.push_back("");
		MSA curr_pairwise = MSA(aligned_seqs);
		return curr_pairwise;
	}
	if(m == 0)
	{
		string retB = B;
		string retA = "";
		for(int j=0; j<n; j++)
		{
			retA += '-';
		}
		vector<string> aligned_seqs;
		aligned_seqs.push_back(retA);
		aligned_seqs.push_back(retB);
		MSA curr_pairwise = MSA(aligned_seqs);
		return curr_pairwise;
	}
	if(n == 0)
	{
		string retA = A;
		string retB = "";
		for(int i=0; i<m; i++)
		{
			retB += '-';
		}
		vector<string> aligned_seqs;
		aligned_seqs.push_back(retA);
		aligned_seqs.push_back(retB);
		MSA curr_pairwise = MSA(aligned_seqs);
		return curr_pairwise;
	}
	// end deal with edge cases where one or mor of the sequences is of size zero

	vector<vector<int>> H_mat(m+1, vector<int>(n+1));
	vector<vector<int>> M_mat(m+1, vector<int>(n+1));
	vector<vector<int>> I_mat(m+1, vector<int>(n+1));
	vector<vector<int>> J_mat(m+1, vector<int>(n+1));

	//matrices initialization
	for(int j=0; j<=n; j++)
	{
		H_mat[0][j] = - _gap_open - _gap_extend *(j - 1);
		M_mat[0][j] = - (INT_MAX /2);
		I_mat[0][j] = - (INT_MAX /2);
		I_mat[1][j] = - _gap_open;
		J_mat[0][j] = - _gap_open - _gap_extend *(j - 1);
	}

	for(int i=0; i<=m; i++)
	{
		H_mat[i][0] = - _gap_open - _gap_extend *(i - 1);
		M_mat[i][0] = - (INT_MAX /2);
		I_mat[i][0] = - _gap_open - _gap_extend *(i - 1);
		J_mat[i][0] = - (INT_MAX /2);
		J_mat[i][1] = - _gap_open;
	}

	H_mat[0][0] = 0;
	M_mat[0][0] = 0;
	I_mat[0][0] = - (INT_MAX /2);
	J_mat[0][0] = - (INT_MAX /2);
	// end initialization

	//score computation - dynamic programming
	for(int i=1; i<=m; i++)
	{
		for(int j=1; j<=n; j++)
		{
			//int S = (A[i-1] == B[j-1]) ? _match_score : -_mismatch_score; //this should be changed to use the blosum
			int S = _simMat.get_pair_score(A[i-1],B[j-1]);
			int M_diag = M_mat[i-1][j-1] + S;
			int I_diag = I_mat[i-1][j-1] + S;
			int J_diag = J_mat[i-1][j-1] + S;
			M_mat[i][j] = max(M_diag, max(I_diag, J_diag));
			
			int M_up = M_mat[i-1][j];
			int I_up = I_mat[i-1][j];
			I_mat[i][j] = max((M_up - _gap_open),(I_up - _gap_extend));

			int M_left = M_mat[i][j-1];
			int J_left = J_mat[i][j-1];
			J_mat[i][j] = max((M_left - _gap_open),(J_left - _gap_extend));

			H_mat[i][j] = max(M_mat[i][j], max(I_mat[i][j], J_mat[i][j]));
		}
	}
	// end score computation - dynamic programming

	// traceback
	string retA, retB;
    stack<char> SA, SB;

    int ii = m;
	int jj = n;

    while (ii != 0 || jj != 0)
    {
        if (ii == 0)
        {
            SA.push('-');
            SB.push(B[jj-1]);
            jj--;
        }
        else if (jj == 0)
        {
            SA.push(A[ii-1]);
            SB.push('-');
            ii--;
        }
        else
        {
			int curr_H_val = H_mat[ii][jj];
			if(curr_H_val == M_mat[ii][jj])
			{
				//go diag
				SA.push(A[ii-1]);
                SB.push(B[jj-1]);
                ii--; jj--;

			}
			else if(curr_H_val == I_mat[ii][jj])
			{
				//take from A
				SA.push(A[ii-1]);
                SB.push('-');
                ii--;
			}
			else
			{
				//take from B
				SA.push('-');
                SB.push(B[jj-1]);
                jj--;
			}
        }
    }

	while (!SA.empty())
    {
        retA += SA.top();
        retB += SB.top();
        SA.pop();
        SB.pop();
    }
	
	vector<string> aligned_seqs;
	aligned_seqs.push_back(retA);
	aligned_seqs.push_back(retB);
	MSA curr_pairwise = MSA(aligned_seqs);
    return curr_pairwise;

}