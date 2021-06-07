/* 
   INDELible V1.03
    "A comprehensive and flexible simulator of molecular sequence evolution"
    Copyright (C) 2009 William Fletcher

    If using this software please cite the following reference:

    Fletcher, W. and Yang, Z. 2009. 
	"INDELible: a flexible simulator of biological sequence evolution." 
	Mol. Biol. and Evol. (in press). 
 
    If you need to contact me with any queries, questions, problems or (god-forbid) bugs
    then please go to http://abacus.gene.ucl.ac.uk/software/indelible/bug.php

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    Please visit http://www.gnu.org/licenses/ for a full copy of the license.
*/

             

#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>  //???
//#include "paml.h"
#include "models.h"

using namespace std;

#pragma warning(disable:4786)	


vector<vector<double> > totalusermodels;



// functions below are all of form int double vector<double>&
// this is so generic pointer to these functions works

int choosenewzipf(int M, double a, vector<double> &z) {int u; do {u=newZipf(a);} while(u>M || u<1); return u;} //do {u=newZipf(a);} while(u>M); return u;}
int chooseoldzipf(int M, double a, vector<double> &z) {int u; do {u=oldZipf(a);} while(u>M || u<1); return u;} //do {u=oldZipf(a);} while(u>M); return u;}

	int choosenewNB(int r, double q, vector<double> &z) {return randnegbin(r,q);}
	int chooseoldNB(int r, double q, vector<double> &z) {return oldrandnegbin(r,q);}
		
	int userrandomsize(int x, double y, vector<double> &usermodel)
	{
		double myrand=mtrand1();

		for(int i=0; i<usermodel.size(); i++) if(myrand<usermodel.at(i)) return i+1;

		return usermodel.size();
	}

int f(){return 1;}
int f2(double x, int y) {return int(x);}

 int type;	// 1 for NUCLEOTIDE, 2 for AMINOACID, 4 for CODON


extern bool controldebug;



void controlerrorprint2(string blocktype, string blockname, string commandname, string instring, string myline)
{
	// this functions provides a framework for a standard output of an error to the screen and log file
	// mostly it is just formatting the border of the error box and white space etc

//	cout<<"\r                                                                      \n ERROR occurred in "<<blocktype<<" block "<<blockname<<".\n See below for details or consult indelibleLOG.txt                       ";
	cout<<"\r                                                                      \n";

	vector<string> toprint;
	
	string startline="ERROR in "+blocktype+" block ";
	if(blockname!="")   startline+=blockname;
	if(commandname!="") startline+=" in command "+commandname;

//	stringstream fd1; fd1<<linecount; string fd1s=fd1.str();
//	startline+=fd1s; startline+=" of control file:";
	
	toprint.push_back(startline);

	string tempstring;
	char c;
	int themaxsize=startline.size();
	
	for(int j0=0; j0<instring.size(); j0++) 
	{
		c=instring[j0]; 
		if(c=='\n') 
		{
			toprint.push_back(tempstring);
			if(themaxsize<tempstring.size()) themaxsize=tempstring.size();
			tempstring="";
		} 
		else tempstring+=c;
	} 
	toprint.push_back(tempstring);if(themaxsize<tempstring.size()) themaxsize=tempstring.size();
	
	if(myline!="") {
		string endline="Last Input read was: ";
		endline+=myline;
		if(themaxsize<endline.size()) themaxsize=endline.size();
		toprint.push_back(endline);
	}

	cout<<endl<<endl; 
	//(*indelibleLOG)<<endl;

	for(int i0=0; i0<toprint.size(); i0++)
	{
		string tempstring2=toprint.at(i0);
		for(int h1=tempstring2.size(); h1<themaxsize; h1++) tempstring2+=" ";
		toprint.at(i0)=tempstring2;
	}

	cout<<endl<<" +";
	//(*indelibleLOG)<<endl<<"+";  
	for(int i1=0; i1<themaxsize+2; i1++) {
		//cout<<"-";
		//(*indelibleLOG)<<"-";
	}
	cout<<"+"<<endl;	
	//(*indelibleLOG)<<"+"<<endl;
	
	for(int i2=0; i2<toprint.size(); i2++) {
		cout<<" | "<<toprint.at(i2)<<" |"<<endl;
		//(*indelibleLOG)<<"| "<<toprint.at(i2)<<" |"<<endl;
	}
	
	cout<<" +";
	//(*indelibleLOG)<<"+"; 
	for(int i3=0; i3<themaxsize+2; i3++) {
		//cout<<"-";
		//(*indelibleLOG)<<"-";
	}
	cout<<"+"<<endl<<endl;
//	(*indelibleLOG)<<"+"<<endl<<endl;

	
}


vector<int> allowedcodes(int gencode)
{
	vector<int> allowedlist;  

	if(gencode==2 || gencode==6 || gencode ==14 || gencode==22 || gencode==23) allowedlist.push_back(gencode);

	else if(gencode==15 || gencode==16) {allowedlist.push_back(15); allowedlist.push_back(16);}

	else if(gencode==1 || gencode==11 || gencode==12) {allowedlist.push_back(1);allowedlist.push_back(11);allowedlist.push_back(12);}

	else {allowedlist.push_back(3);allowedlist.push_back(4);allowedlist.push_back(5);allowedlist.push_back(9);
	allowedlist.push_back(10);allowedlist.push_back(13);allowedlist.push_back(21);}

	return allowedlist;
}

vector<int> getstops(int geneticcode)
{
	// finds the stops in genetic codes listed above

	vector<int> stops; 

	for(int i=0; i<64; i++) if(GeneticCodeTable[geneticcode][i]=='*') stops.push_back(i);

	return stops;
}

void enforcestops(int geneticcode, vector<double> &basefreqs)
{
	// makes base frequencies of stop codons equal to zero

	vector<int> stops=getstops(geneticcode);

	for(int i=0; i<stops.size(); i++) basefreqs.at(stops.at(i))=0;
}

int querystops(int geneticcode, vector<double> &basefreqs)
{
	// returns non zero base frequencies at stop codons

	vector<int> stops=getstops(geneticcode);

	for(int i=0; i<stops.size(); i++) if(basefreqs.at(stops.at(i))!=0) return stops.at(i);

	return -1;
}
//////////////////////////////////////////////////////////////////////////////////////////


model::model(int mymodelpos, int &mytype, string &myname, int &mymodelnumber, int &mygeneticcode,
		bool &mycopiedmodel, double &myinsertrate, double &mydeleterate, double &myalpha, double &mypinv, 
		int &myngamcat, double mycodonrates[], vector<double> &mybasefreqs, vector<double> &myrootbasefreqs, 
		vector<double> &myinsertfreqs, vector<double> &myparams, vector<double> &aamodel, indelmodel &insertmodel,
		indelmodel &deletemodel)
	{


		// set deletion model
		delmeansize=deletemodel.meansize;

		if(deletemodel.type == 0 || deletemodel.type == 3 || deletemodel.type == 12 ) 
		{
			delV=deletemodel.usermodel;
			delrandomsize=&userrandomsize;
		}
		else if(deletemodel.type == 2 || deletemodel.type == 13) 
		{
			//Zipfian model

			delI=deletemodel.M;   
			delD=deletemodel.a;

			if(deletemodel.type==2) 
			{
				double v=1, q=delD;

				q1=1-q; q2=1/q1;

				Hx0 = H(0.5,q1,q2,v) - exp( log(v) * (-q) );

				s = 1 - H1( H(1.5,q1,q2,v) -exp( log(v+1) * (-q) )  ,q1,q2,v );

				Himax = H( imax + 0.5 ,q1,q2,v);

				delrandomsize=&choosenewzipf;
			}
			else					delrandomsize=&chooseoldzipf;
		}
		else if(deletemodel.type == 1 || deletemodel.type == 11) 
		{
			// Negative Binomial or Geometric

			delI=deletemodel.r; 
			delD=deletemodel.q;

			if(deletemodel.type==1) delrandomsize=&choosenewNB;
			else					delrandomsize=&chooseoldNB;
		}


		// set insertion model
		if(insertmodel.type == 0 || insertmodel.type == 3 || insertmodel.type == 12 ) 
		{
			insV=insertmodel.usermodel;

			insrandomsize=&userrandomsize;
		}
		else if(insertmodel.type == 2 || insertmodel.type == 13) 
		{
			//Zipfian model

			insI=insertmodel.M; 
			insD=insertmodel.a;

			if(insertmodel.type==2) 
			{
				double v=1, q=insD;

				q1=1-q; q2=1/q1;

				Hx0 = H(0.5,q1,q2,v) - exp( log(v) * (-q) );

				s = 1 - H1( H(1.5,q1,q2,v) -exp( log(v+1) * (-q) )  ,q1,q2,v );

				Himax = H( imax + 0.5 ,q1,q2,v);

				insrandomsize=&choosenewzipf;
			}
			else					insrandomsize=&chooseoldzipf;
		}
		else if(insertmodel.type == 1 || insertmodel.type == 11) 
		{
			// Negative Binomial or Geometric

			insI=insertmodel.r; 
			insD=insertmodel.q;

			if(insertmodel.type==1) insrandomsize=&choosenewNB;
			else					insrandomsize=&chooseoldNB;
		}





		continuousgamma=false;
		numberofsiteclasses=1;

		insertrate	=myinsertrate;
		deleterate	=mydeleterate;
		subrate=1;
		//indelrate=insertrate+deleterate;
		//subrate		=1-insertrate-deleterate;

		modelpos=mymodelpos;

		error=1;
		type=mytype;
		name=myname;
		modelnumber=mymodelnumber;
		geneticcode=mygeneticcode;
		copiedmodel=mycopiedmodel;

		alpha=myalpha;
		pinv=mypinv;
		ngamcat=myngamcat;

		medianORmean=0;

		if(type!=3)
		{
			// set up discrete gamma and/or pinv --> actual rates are picked in SetSiteRates in main skeleton file.
			if(alpha>0) 
			{
				if(ngamcat==0) {Rrates.push_back(1); cumfreqs.push_back(1);  continuousgamma=true;  }
				else
				{
					DiscreteGamma(cumfreqs,Rrates,pinv, alpha,alpha,ngamcat,medianORmean);
				}
			}
			else
			{
				if(pinv>0) {Rrates.push_back(0); Rrates.push_back(1/(1-pinv));  cumfreqs.push_back(pinv); cumfreqs.push_back(1); }
				else {Rrates.push_back(1); cumfreqs.push_back(1); }
			} 

			numberofsiteclasses=Rrates.size();
		}

		// set up codon position specific relative substitution rates.
		// has no effect if gamma model used - that supercedes this setting.

		if( alpha==0 && (mycodonrates[0]!=1 || mycodonrates[1]!=1 || mycodonrates[2]!=1) ) codonratestrue=true; else codonratestrue=false;
		codonrates[0]=mycodonrates[0];
		codonrates[1]=mycodonrates[1];
		codonrates[2]=mycodonrates[2];

		if(type==1 && !copiedmodel)
		{
			// for DNA
			if(mybasefreqs.empty()) 
			{
				// if base frequencies are empty
				if(modelnumber%2==1 && modelnumber<16 ) //&& !copiedmodel)  
				{
					//when they shouldn't be give a warning
					stringstream med; med<<modelnumber; string modelnumberX=med.str();
					controlerrorprint2("[MODEL]",name,"basefreqs","A model with unequal base frequencies was chosen: model "+modelnumberX+"\nbut you have not specified base frequencies.  They will be equal.",""); 
				}
				// make equal frequencies either way
				makeequalfreqs(type,basefreqs);
			}
			else
			{
				// if base frequencies are given
				if(modelnumber%2==0 && modelnumber<16 ) // && !copiedmodel)  
				{
					// on a model that wants equal frequencies
					// force frequencies to be equal and give warning
					makeequalfreqs(type,basefreqs);
					stringstream med; med<<modelnumber; string modelnumberX=med.str();
					controlerrorprint2("[MODEL]",name,"basefreqs","A model with equal base frequencies was chosen: model "+modelnumberX+"\nbut you have specified base frequencies.  They will be set equal.",""); 
				}
				// otherwise set frequencies as given
				else basefreqs=mybasefreqs;  
			}


		}

		if(type==2 && !mybasefreqs.empty()) basefreqs=mybasefreqs;  // this will force +F models		


		if(type==3) 
		{
			if(mybasefreqs.empty()) makeequalfreqs(type,basefreqs);			// fequal frequencies

			else if(mybasefreqs.size()==4) basefreqs=fx(mybasefreqs,1);		// f1x4 frequencies

			else if(mybasefreqs.size()==12) basefreqs=fx(mybasefreqs,3);	// f3x4 frequencies

			else if(mybasefreqs.size()==64) 								// fcodon frequencies
			{
				basefreqs=mybasefreqs;			

				testmyfreqs(basefreqs,"[basefreq]");
			}
			else cout<<"INTERNAL ERROR 463"<<endl;

		}




		// these make the correct Q matrix for a given type and model number

		// make Qvec and Jvec for nucleotide models
		if(type==1) 
		{

			//	for(int y=0; y<Rrates.size(); y++)
			//	{ 
			Qvec=getDNA(name,myparams,basefreqs, modelnumber);   Qvecs.push_back(Qvec);
			getJvec(0, /*Rrates.at(y)*/name,myrates,Qvec,Jvec,basefreqs); 

			if(Rrates.size()!=0) Jvecs.assign(Rrates.size(),Jvec); else Jvecs.push_back(Jvec);

			//	}
		}

		// make Qvec and Jvec for amino acid models
		if(type==2) 
		{
			//	for(int y=0; y<Rrates.size(); y++)
			//	{ 
			Qvec=getAA( name,myparams,basefreqs, modelnumber, aamodel); Qvecs.push_back(Qvec);

			getJvec(0,name,myrates,Qvec,Jvec,basefreqs); 

			if(Rrates.size()!=0) Jvecs.assign(Rrates.size(),Jvec); else Jvecs.push_back(Jvec);

			//	}
		}

		if(type==3) 
		{

			/*
			(*) Codon models for variable dN/dS ratios among sites
			(com.nrate includes kappa & omega) (see also CDFdN_dS)

			NSsites          npara

			0  one-ratio     0:    one ratio for all sites
			1  neutral       1:    p0 (w0=0, w1=1)
			2  selection     3:    p0, p1, w2 (w0=0, w1=1)
			3  discrete      2K-1: p0,p1,..., and w0,w1,...
			4  freqs         K:    p's (w's are fixed)
			5  gamma         2:    alpha, beta
			6  2gamma        4:    p0, alpha1,beta1, alpha2=beta2
			7  beta          2:    p_beta, q_beta
			8  beta&w        4:    p0, p_beta, q_beta, w estimated
			9  beta&gamma    5:    p0, p_beta, q_beta, alpha, beta
			10  beta&1+gamma  5:    p0, p_beta, q_beta, alpha, beta (1+gamma used)
			11  beta&1>normal 5:    p0, p_beta, q_beta, mu, s    (normal truncated w>1)
			12  0&2normal     5:    p0, p1, mu2, s1, s2
			13  3normal       6:    p0, p1, mu2, s0, s1, s2
			14  M8a:beta&w=1  3:    p0, p_beta, q_beta, w=1 fixed
			15  M8a:beta&w>=1 4:    p0, p_beta, q_beta, w>=1 estimated


			*/		

			// CODON MODELS - numbered 0 to 15

			// only 3 (M3) and 14, 15 are in proper use.  M0-M13 can all be expressed as M3

			// 0 is the codon model of Goldman and Yang 1994 ; Goldman, N., and Z. Yang. 1994. A codon-based model of nucleotide substitution for protein-coding DNA sequences. Molecular Biology and Evolution 11:725-736.

			// 1 to 13 are the codon sites-models of :Yang, Z., Nielsen, R., Goldman, N. and Pedersen, A-M. K. (2000) Codon-Substitution Models for Heterogeneous Selection Pressure at Amino Acid Sites. Genetics 155: 431-439 (May 2000)

			// 14 and 15 are the empirical codon models ECM : Kosiol, C., Holmes, I. and Goldman, N. (2007) An Empirical Codon Model for Protein Sequence Evolution.  Molecular Biology and Evolution 24(7): 1464-1479.


			if(modelnumber==14) 
			{
				Qvec=getECMr();
				d(Qvec,scalefactors.at(0)); 
				Qvecs.push_back(Qvec); 
				getJvec(0,name,myrates,Qvec,Jvec,basefreqs); 
				Jvecs.push_back(Jvec); 
				cumfreqs.push_back(1);
			}
			else if(modelnumber==15) 
			{
				Qvec=getECMu();
				d(Qvec,scalefactors.at(0));
				Qvecs.push_back(Qvec); 
				getJvec(0,name,myrates,Qvec,Jvec,basefreqs); 
				Jvecs.push_back(Jvec); 
				cumfreqs.push_back(1); 
			}
			else
			{
				// all other models  follow same pattern of generation.
				// for each site class make the Qvec for that omega, and calculate the scale factor for that Qvec
				// then calculate overall scale factor from the other scale factors and the proportions
				// divide every Qvec for every site class by the overall scale factor.using d()
				// then calculate Jvecs for each site class' Qvec

				if(myparams.size()==0)
				{
					myparams.push_back(1); myparams.push_back(1);
					controlerrorprint2("[MODEL]",name,"","No kappa/omega have been defined so they have been set equal to 1.","");

				}
				double kappa=myparams.at(0);
				double omega;

				if(modelnumber==0) 
				{
					// (Goldman and Yang, 1994)  

					cumfreqs.push_back(1);
					omega=myparams.at(1);  myomegas.push_back(omega);
					Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega);
					d(Qvec,scalefactors.at(0));
					Qvecs.push_back(Qvec);

					getJvec(0,name,myrates,Qvec,Jvec,basefreqs);
					Jvecs.push_back(Jvec); 
				}

				else if(modelnumber==1)
				{
					//cout<<myparams.at(1)<<" 1 1 "<<endl;
					double p0=myparams.at(1), p1=1-p0;

					cumfreqs.push_back(p0);
					omega=myparams.at(2);   myomegas.push_back(omega);
					Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);					

					//cout<<1-myparams.at(1)<<" 1 2 "<<endl;
					cumfreqs.push_back(p1);
					omega=1;   myomegas.push_back(omega);
					Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);

					//double S=(p0*scalefactors.at(0))+(p1*scalefactors.at(1));

					double S=0; for(int y1=0; y1<cumfreqs.size(); y1++) 
					{
						//ccout<<endl<<scalefactors.at(y1)<<"  "<<cumfreqs.at(y1)<<"  "<<scalefactors.at(y1)*cumfreqs.at(y1)<<endl;

						S+=(scalefactors.at(y1)*cumfreqs.at(y1));
					}

					d(Qvecs.at(0),S); getJvec(S,name,myrates,Qvecs.at(0),Jvec,basefreqs); Jvecs.push_back(Jvec); 					
					d(Qvecs.at(1),S); getJvec(S,name,myrates,Qvecs.at(1),Jvec,basefreqs); Jvecs.push_back(Jvec); 					

				}
				else if(modelnumber==2)
				{
					//	cout<<myparams.at(1)<<" 1 1 "<<endl;
					//	cout<<myparams.at(2)<<" 1 2 "<<endl;
					//	cout<<myparams.at(3)<<" 1 3 "<<endl;
					//	cout<<1-myparams.at(1)-myparams.at(2)<<" 1 4 "<<endl;


					double p0=myparams.at(1), p1=myparams.at(2), p2=1-p0-p1;
					cumfreqs.push_back(p0);
					omega=myparams.at(3); //cout<<p0<<"  "<<omega<<endl; 
					myomegas.push_back(omega);
					Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);					

					cumfreqs.push_back(p1);
					omega=1;  //cout<<p1<<"  "<<omega<<endl;
					myomegas.push_back(omega);
					Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);					

					cumfreqs.push_back(p2);
					omega=myparams.at(4); // cout<<p2<<"  "<<omega<<endl;
					myomegas.push_back(omega);
					Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);

					double S=0; for(int y1=0; y1<cumfreqs.size(); y1++) S+=(scalefactors.at(y1)*cumfreqs.at(y1));

					//for(int t1=0; t1<Qvecs.size(); t1++) {getJvec(S,name,myrates,Qvecs.at(t1),Jvec,basefreqs); Jvecs.push_back(Jvec); 	}	
					d(Qvecs.at(0),S); getJvec(S,name,myrates,Qvecs.at(0),Jvec,basefreqs); Jvecs.push_back(Jvec); 
					d(Qvecs.at(1),S); getJvec(S,name,myrates,Qvecs.at(1),Jvec,basefreqs); Jvecs.push_back(Jvec); 
					d(Qvecs.at(2),S); getJvec(S,name,myrates,Qvecs.at(2),Jvec,basefreqs); Jvecs.push_back(Jvec); 

				}
				else if(modelnumber==3)
				{

					int mybit=myparams.size()/2; 
					double sum=0;
					for(int yf=1; yf<mybit; yf++)   
					{
						//cout<<yf<<" "<<"1"<<endl;
						omega=myparams.at(yf+mybit-1);   myomegas.push_back(omega);
						//cout<<yf<<" "<<"2 "<<myparams.at(yf+mybit-1)<<endl;
						Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);


						//cout<<yf<<" "<<"3"<<endl;
						sum+=myparams.at(yf);
						//cout<<yf<<" "<<"4"<<endl;
						cumfreqs.push_back(myparams.at(yf));
						//cout<<yf<<" "<<"5 "<<myparams.at(yf)<<endl;
					}

					if(sum<=1)
					{
						//cout<<"BLAH"<<endl;
						omega=myparams.at(2*mybit-1);   myomegas.push_back(omega);
						Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);


						cumfreqs.push_back(1-sum);
						//cout<<"BLAH "<<1-sum<<" "<<omega<<endl;
					}
					else cout<<"Error in sum of category frequencies in codon model 3"<<endl;

					double S=0; for(int y1=0; y1<cumfreqs.size(); y1++) S+=(scalefactors.at(y1)*cumfreqs.at(y1));
					for(int yf0=0; yf0<Qvecs.size(); yf0++)   {d(Qvecs.at(yf0),S); getJvec(S,name,myrates,Qvecs.at(yf0),Jvec,basefreqs); Jvecs.push_back(Jvec); }
				}




				// model numbers 4 to 13 are not used any mor
				// N.B.  THERE IS NO VECTOR OF DOUBLES CALLED myomegas IN THIS SECTION _ JUST IN CASE OF CRASHES
				else if(modelnumber==4)
				{

					//as K-1 for 4, 2K-1 for 3, but is K and 2K because of kappa, difference of K
					// for modelnumber 4, K=5 means size is 5 (kappa, p0, p1, p2, p3) with p4=1-p1-p2-p3 etc
					double sum=0, mysize=myparams.size()-2;  // mysize is 3
					for(int yf=1; yf<mysize+2; yf++)   //from 1 in the list (p0 after kappa) to mysize+2=5 goes up to p3
					{
						//cout<<yf<<" "<<"1"<<endl;
						omega=(yf-1)/mysize;  //omega is 0, 1/3, 2/3. 1
						//cout<<yf<<" "<<"2 "<<endl;
						Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);


						//cout<<yf<<" "<<"3"<<endl;
						sum+=myparams.at(yf);
						//cout<<yf<<" "<<"4"<<endl;
						cumfreqs.push_back(myparams.at(yf));
						//cout<<yf<<" "<<"5 "<<myparams.at(yf)<<"  "<<omega<<endl;

					}
					if(sum<=1)
					{
						//cout<<"BLAH"<<endl;
						omega=mysize; //omega is 3
						Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);


						cumfreqs.push_back(1-sum);
						//cout<<"BLAH "<<1-sum<<" "<<omega<<endl;

					}
					else cout<<"Error in sum of category frequencies in codon model 4"<<endl;

					double S=0; for(int y1=0; y1<cumfreqs.size(); y1++) S+=(scalefactors.at(y1)*cumfreqs.at(y1));
					for(int yf0=0; yf0<Qvecs.size(); yf0++)   {d(Qvecs.at(yf0),S); getJvec(S,name,myrates,Qvecs.at(yf0),Jvec,basefreqs); Jvecs.push_back(Jvec); }

				}
				else 
				{
					if(modelnumber==12 || modelnumber==8)
					{
						if(modelnumber==12)
						{
							omega=0;
							Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);


							cumfreqs.push_back(myparams.at(1));
							for(int hfd=0; hfd<ngamcat; hfd++) cumfreqs.push_back((1-myparams.at(1))/double(ngamcat));
						}
						else
						{
							omega=myparams.at(4);
							Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);


							cumfreqs.push_back(1-myparams.at(1));
							for(int hfd=0; hfd<ngamcat; hfd++) cumfreqs.push_back(myparams.at(1)/double(ngamcat));
						}
					}
					else cumfreqs.assign(ngamcat,1/double(ngamcat));

					//for(int hfd=0; hfd<ngamcat; hfd++) cumfreqs.push_back(1/ngamcat);

					vector<double> output;

					//double *mypar; mypar=new double[myparams.size()];

					double mypar[10]={0,0,0,0,0,0,0,0,0,0};

					//mypar[0]=0;
					for(int iu=1; iu<myparams.size(); iu++) {mypar[iu-1]=myparams.at(iu); }  

					// this function of Ziheng's calculates the discrete rates for different site classes from the model parameters

					DiscreteNSsites(mypar, ngamcat, modelnumber, output);

					for(int i=0; i<ngamcat; i++)
					{
						omega=output.at(i);
						Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);

					}

					double S=0; for(int y1=0; y1<cumfreqs.size(); y1++) S+=(scalefactors.at(y1)*cumfreqs.at(y1));
					for(int yf0=0; yf0<Qvecs.size(); yf0++)   {d(Qvecs.at(yf0),S);  getJvec(S,name,myrates,Qvecs.at(yf0),Jvec,basefreqs); Jvecs.push_back(Jvec); }

				}//end of model>4 bracket

			} //end of solitary else bracket


			//make cumfreqs cumulative;
			double mysum=0;
			vector<double> blah=cumfreqs;

			numberofsiteclasses=cumfreqs.size();

			for(int gfv=1; gfv<numberofsiteclasses; gfv++) cumfreqs.at(gfv)+=cumfreqs.at(gfv-1);

			for(int f=0; f<numberofsiteclasses; f++) Rrates.push_back(1);

		}//end of type==3 bracket

		// the Qvec/Jvec must be generated before the root/insert freqs are set in case the base frequencies come from 
		// the model like in empirical sub models (codon/protein) or when using the UNREST model for DNA.

		if(myrootbasefreqs.empty()) rootbasefreqs=basefreqs; else {rootbasefreqs=myrootbasefreqs; if(type==3) testmyfreqs(rootbasefreqs,"[rootbasefreq]"); }

		if(myinsertfreqs.empty()) insertfreqs=basefreqs; else {insertfreqs=myinsertfreqs; if(type==3) testmyfreqs(insertfreqs,"[ibasefreq]"); }


	}


///////////////////
/*
double main(int argc, char* argv[])
{	

	return 0;
}


*/

