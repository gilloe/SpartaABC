#ifndef ___RANDOMS_H
#define ___RANDOMS_H


#include "MersenneTwister.h"

extern MTRand mtrand1;
extern MTRand mtrand2;

extern double myrand;
extern int idum;
int randnegbin(int r, double q);
int Zipf(double q, double v);
#define newZipf(a) Zipf(a,1)
int oldZipf(double a);
int oldrandnegbin(int r, double q);
#define H(x,q1,q2,v) exp(q1*log(v+x))*q2
#define H1(x,q1,q2,v) -v+exp(q2*log((1-q)*x))
extern double Zq1,  Zq2,  ZHx0,  Zs,  ZHimax;

double expmyrand(double myrate);





#endif


