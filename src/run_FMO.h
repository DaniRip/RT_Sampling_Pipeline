// Existing headers
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <deque>
#include <algorithm>

// Merging properties
#include <fstream>
#include <ctime>
#include <math.h>
#include <new>
#include <vector>

// cplex header files
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilocsvreader.h>
#include <ilconcert/iloexpression.h>

// 
#include "MatlabDataArray.hpp"

using namespace std;  // use with caution in headers!

void FMO(float &objVal, vector<float> &intensityVals, string &status, float &runtime, bool *objectives, bool *cons, matlab::data::TypedArray<double> &d_target, matlab::data::TypedArray<double> &d_OAR, int num_target_voxels, int num_OAR_voxels, int num_beamlets, double target_dose);
//void readIn(float &D_Matrix[][], matlab::data::TypedArray &downsampled, int downsampled_beamlets, int downsampled_voxels)
static void assignValues(IloEnv env, IloNumArray vDose, IloNumArray2 D_tumor, IloNumArray2 D_OAR, matlab::data::TypedArray<double> &d_target, matlab::data::TypedArray<double> &d_OAR, int num_target_voxels, int num_OAR_voxels, int num_beamlets, double target_dose);
static void minimizeDose(IloInt numVoxels, IloNumArray2 D, IloNumVarArray w, IloExpr dose, string name);
static void enforceMinTumorDose(IloNumArray2 D, IloNumVarArray w, IloNumArray vDose, IloEnv env, IloModel model);
static void minimizeSqDose(IloInt numVoxels, IloNumArray2 D, IloNumVarArray w, IloExpr& totalPenalty, IloEnv env, IloNumArray vDose, string name);
static void minimizeAbsDose(IloInt numVoxels, IloNumArray2 D, IloNumVarArray w, IloExpr dose, IloEnv env, IloModel model, IloNumArray vDose, string name);

#pragma once
