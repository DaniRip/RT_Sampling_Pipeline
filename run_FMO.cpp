#include "mex.hpp"
#include "mexAdapter.hpp"
#include "run_FMO.h"

// to compile call the line below in matlab
// mex('-IC:\Program Files\IBM\ILOG\CPLEX_Studio201\cplex\include','-IC:\Program Files\IBM\ILOG\CPLEX_Studio201\concert\include','-LC:\Program Files\IBM\ILOG\CPLEX_Studio201\concert\lib\x64_windows_msvc14\stat_mda','-LC:\Program Files\IBM\ILOG\CPLEX_Studio201\cplex\lib\x64_windows_msvc14\stat_mda','-lcplex2010.lib','-lilocplex.lib','-lconcert.lib','run_FMO.cpp')
class MexFunction : public matlab::mex::Function {
public:
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
	{
		// This function runs when run_FMO is called in MATLAB
		matlab::data::ArrayFactory factory;

		matlab::data::CharArray charVector = inputs[0];
		string samplingType = charVector.toAscii();

		string y = "0";
      charVector = inputs[1];
		y = charVector.toAscii();
      double voxels = stod(y);

		charVector = inputs[2];
		y = charVector.toAscii();
      double layers = stod(y);

      matlab::data::TypedArray<bool> obj_array = move(inputs[3]);
      bool objectives[3];
      int count = 0;
      for (auto& val : obj_array)
      {
         objectives[count] = val;
         count++;
      }

      matlab::data::TypedArray<double> d_target = move(inputs[4]); //2d array, not in column major order
      matlab::data::TypedArray<double> d_OAR = move(inputs[5]);

      matlab::data::Array doubleArray= move(inputs[6]);
      int num_target_voxels = doubleArray[0];
      doubleArray= move(inputs[7]);
      int num_OAR_voxels = doubleArray[0];
      doubleArray= move(inputs[8]);
      int num_beamlets = doubleArray[0];
      doubleArray= move(inputs[9]);
      double target_dose = doubleArray[0];

		float objVal, runtime;
      vector<float> intensityVals;
		string status;
		FMO(objVal, intensityVals, status, runtime, objectives, d_target, d_OAR, num_target_voxels, num_OAR_voxels, num_beamlets, target_dose);

      string s = "Status = " + status;
      outputs[0] = factory.createCharArray(s);

      s = "Objective = " + to_string(objVal);
		outputs[1] = factory.createCharArray(s);

      s = "Runtime = " + to_string(runtime) + " sec";
		outputs[2] = factory.createCharArray(s);

      //s = "Voxels(input) = " + to_string(voxels);
      s = "num_target_voxels = " + to_string(num_target_voxels);
		outputs[3] = factory.createCharArray(s);

      //s = "Layers(input) = " + to_string(layers);
      /*for (int x=0;x<5;x++){
         s = s + to_string(d_target[x]);
      }*/
      //s = "" + to_string(d_target[0][0]);
		//outputs[4] = factory.createCharArray(s);

      //Array a = factory.createArray<float>({1,intensityVals.size()}, {intensityVals});
      //outputs[2] = a;
	}
};

typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloArray<IloNumVarArray> > IloNumVarArray3;
typedef IloArray<IloBoolVarArray> IloBoolVarArray2;
typedef IloArray<IloArray<IloBoolVarArray> > IloBoolVarArray3;

/* Stripped down version of Dani's code for c++ template for William! */

IloInt NUMBEAMLETS, NUMTUMORVOXELS, NUMOARVOXELS; // Globals to assign programatically

// Body of the code
void FMO(float &objVal, vector<float> &intensityVals, string &status, float &runtime, bool *objectives, matlab::data::TypedArray<double> &d_target, matlab::data::TypedArray<double> &d_OAR, int num_target_voxels, int num_OAR_voxels, int num_beamlets, double target_dose) {

   IloEnv env; // Initialize CPLEX model environment
   try{
      IloNumArray2 D_OAR(env), D_tumor(env);    // indexes of D: [numvox][numbeamlets]
      IloNumArray vDose(env);
      // You need to assign the Dij matrices into these guys!
      IloInt b, v;

      //float d_[downsampled_beamlets][downsampled_voxels];
      //readIn(D_Matrix, downsampled, downsampled_beamlets, downsampled_voxels);
      //readIn();
      assignValues(env, vDose, D_tumor, D_OAR, d_target, d_OAR, num_target_voxels, num_OAR_voxels, num_beamlets, target_dose);

      // Initialize CPLEX model
      IloModel model(env);

      // Declare variables
      IloNumVarArray w(env); // The intensity (decision variable)
      for (b=0; b < NUMBEAMLETS; b++){
         w.add(IloNumVar(env,0,IloInfinity,ILOFLOAT));
      }

      // Add basic FMO constraint:
      enforceMinTumorDose(D_tumor, w, vDose, env, model);

      IloExpr ObjDose(env);
      IloExtractable OriginalObj;
      // Uncomment at least one obj. function.
      if (objectives[0] == 1 && objectives[1] == 1)
      {
         minimizeDose(NUMOARVOXELS, D_OAR, w, ObjDose, "OAR"); //1
         minimizeDose(NUMTUMORVOXELS, D_tumor, w, ObjDose, "tumor"); //2
      }
      else if (objectives[1] == 1)
      {
         minimizeDose(NUMTUMORVOXELS, D_tumor, w, ObjDose, "tumor"); //2
      }
      else
      {
         minimizeDose(NUMOARVOXELS, D_OAR, w, ObjDose, "OAR"); //1
      }

      model.add(OriginalObj = IloMinimize(env, ObjDose));
      ObjDose.end();

      // build model
      IloCplex cplex(model);
      cplex.exportModel("basicFMO.lp");
      cplex.setParam(IloCplex::EpGap, 1e-3);
      //cplex.setParam(IloCplex::Param::TimeLimit, 100);

      float* start {new float (cplex.getTime())};
      cout << "Solving, here we go!" << "\n";
      cplex.solve();
      runtime = cplex.getTime()-*start;
      env.out() << endl << "Runtime = " << endl ;
      delete start, runtime;

      env.out() << endl << "Status = " << cplex.getCplexStatus() << endl ;
      status = to_string(cplex.getCplexStatus());

      env.out() << endl << "Objective Value = " << cplex.getObjValue() << endl ;
      objVal = cplex.getObjValue();

      env.out() << endl << "Intensity Values: " << endl ;
      for (b = 0; b < NUMBEAMLETS; b++) {
         cout << cplex.getValue(w[b]) << endl;
         intensityVals.push_back(cplex.getValue(w[b]));
      }
   }
   catch (IloException& ex) {
      cerr << "Error: " << ex << endl;
   }
   catch (...) {
      cerr << "Error" << endl;
   }

   env.end();

   return;
}
/*
void readIn(float &D_Matrix[][], matlab::data::TypedArray &downsampled, int downsampled_beamlets, int downsampled_voxels)
{
   int countrows = 0, countcols = 0;
   for (auto it = downsampled.begin(); it != downsampled.end(); it++)
   {
      if(countrows == downsampled_beamlets)
      {
         countcols++;
         countrows = 0;
      }

      if (countrows < downsampled_voxels)
      {
         D_Matrix[countrows][countcols] = *it;
      }
      countrows++;
   }
}*/

/* Data readin and assignment function */
static void assignValues(IloEnv env, IloNumArray vDose, IloNumArray2 D_tumor, IloNumArray2 D_OAR, matlab::data::TypedArray<double> &d_target, matlab::data::TypedArray<double> &d_OAR, int num_target_voxels, int num_OAR_voxels, int num_beamlets, double target_dose)
{
   // Get size of D matrix
   NUMBEAMLETS = num_beamlets;
   NUMTUMORVOXELS = num_target_voxels;
   NUMOARVOXELS = num_OAR_voxels;

   IloInt i;

   // Add target doses
   IloNumArray input = IloNumArray(env);
   for (int x = 0; x<NUMTUMORVOXELS; x++) {
      input.add(target_dose);
   }
   vDose.add(input);

   for (i=0; i < NUMTUMORVOXELS; i++)
      D_tumor.add(IloNumArray(env));

   for (i=0; i < NUMOARVOXELS; i++)
      D_OAR.add(IloNumArray(env));

   // Assign values from Matlab matrix to CPLEX matrix
   for (int x=0;x<num_target_voxels;x++){
      for (int y=0;y<num_beamlets;y++){
         D_tumor[x].add(IloNum(d_target[x][y]));
      }
   }

   for (int x=0;x<num_OAR_voxels;x++){
      for (int y=0;y<num_beamlets;y++){
         D_OAR[x].add(IloNum(d_OAR[x][y]));
      }
   }

   /*print
   IloNumArray test = IloNumArray(env, 5);
   for (int x=0;x<num_OAR_voxels;x++){
      cout << "\n";
      for (int y=0;y<num_beamlets;y++){
         cout << " " << d_OAR[x][y];
         test.add(d_OAR[x][y]);
      }
   }

   cout << "size = " << test.getSize() << " ";
   for (int x=0;x<test.getSize();x++){
      cout << "\n";
      cout << " " << test[x];
   }*/

   /*D_tumor[0] = IloNumArray(env, NUMBEAMLETS, 0.1, 0.0, 0.07, 0.03, 0.02);
   D_tumor[1] = IloNumArray(env, NUMBEAMLETS, 0.2, 0.15, 0.04, 0.0, 0.0);
   D_tumor[2] = IloNumArray(env, NUMBEAMLETS, 0.0, 0.02, 0.05, 0.0, 0.08);

   D_OAR[0] = IloNumArray(env, NUMBEAMLETS, 0.06, 0.03, 0.01, 0.04, 0.09);*/
}

/* Constraints for restricting the minimum dose on all tumor voxels*/
static void enforceMinTumorDose(IloNumArray2 D, IloNumVarArray w, IloNumArray vDose, IloEnv env, IloModel model)
{
   IloInt b, v;

   cout << "Adding tumor dose-minimum constraint" << endl;

   for (v = 0; v < NUMTUMORVOXELS; v++) {
      IloExpr doseConstTumor(env);
      for (b = 0; b < NUMBEAMLETS; b++) {
         doseConstTumor+=D[v][b]*w[b];
      }
      model.add(doseConstTumor >= vDose[v]);
      doseConstTumor.end();
   }
}

/* Create a dose-minimizing objective function */
static void minimizeDose(IloInt numVoxels, IloNumArray2 D, IloNumVarArray w, IloExpr dose, string name)
{
   cout << "minimizing " << name << endl;

   IloInt b, v;
   IloNum doseAggregate;

   // Efficient objective function read-in method:
   for (b = 0; b < NUMBEAMLETS; b++) {
      doseAggregate = 0;
      for (v = 0; v < numVoxels; v++) {
         doseAggregate += D[v][b]/numVoxels;
      }
      dose += doseAggregate*w[b];
   }

}
