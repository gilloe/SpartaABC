#ifndef ___PAML_H
#define ___PAML_H

#include <string>
#include <vector>
using namespace std;


string newrandomtree(int ntaxa, double birth, double death, double sample, double mut, int randomseed, int option);
int DiscreteGamma (vector<double> &cumfreqK, vector<double> &rK, double pinv,
	double alfa, double beta, int K, int median);
int DiscreteNSsites(double par[], int ngamcat, int model, vector<double> &output);
double rndgamma (double s);

vector<vector<double> > matexp (vector<vector<double> > Qvec, vector<double> &basefreqs,  double t);
vector<vector<double> > PMatQRev(vector<vector<double> > Qvec, vector<double> &basefreqs, double t);



#endif
