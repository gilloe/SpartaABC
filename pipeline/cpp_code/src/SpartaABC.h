#ifndef _SPARTAABC
#define _SPARTAABC
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <chrono>

#include "MSA.h"
#include "Simulation.h"
#include "SpartaABC_options.h"
#include "summary_stats_defs.h"
#include "euclidean_distance.h"
#include "param_priors.h"
#include "summ_stats_wrapper.h"

using namespace std;

MSA simulateSingleMSA(int rootLength, double a_param, double indel_rate_ratio, vector<string> templateIndelibleCtrl);
vector<double> getStatVec (MSA & currMSA);
vector<double> getWeightsVec();
vector<string> readIndelibleTemplateControlFile(string indelibleTemplateCtrlFile);

#endif