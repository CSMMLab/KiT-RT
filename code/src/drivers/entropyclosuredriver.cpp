/*!
 * \file entropyclosuredriver.cpp
 * \brief The main subroutines for driving entropyclosure routines.
 * \author S. Schotthöfer
 *
 */
//#include <iostream>
//#include "drivers/entropyclosuredriver.h"
//
// using namespace std;
//
//
// EntropyClosureDriver::EntropyClosureDriver(Config *settings){
//  _driverSettings = settings;
//  _startTime = 0.0;
//  _timeIter = 0;
//}
//
// void EntropyClosureDriver::StartSolver(){
//  //_startTime = //Get StartTime
//
//  // _settings->Set_StartTime(StartTime);
//
//  /*--- Main external loop of the solver. Runs for the number of time steps required. ---*/
//
//  //if (rank == MASTER_NODE) Later for MPI
//    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;
//
//  //if (rank == MASTER_NODE) Later for MPI
//  {
//    cout << endl <<"Simulation Run using the Entropy Closure Driver" << endl;
//  }
//
//  /*--- Run the problem until the number of time iterations required is reached. ---*/
//  //Currently only steady state
//  {
//
//    /*--- Perform some preprocessing before starting the time-step simulation. ---*/
//
//    Preprocess(_timeIter);
//
//    /*--- Run a time-step iteration of the single-zone problem. ---*/
//
//    Run();
//
//    /*--- Perform some postprocessing on the solution before the update ---*/
//
//    Postprocess();
//
//    /*--- Monitor the computations after each iteration. ---*/
//
//    Monitor(_timeIter);
//
//    /*--- Output the solution in files. ---*/
//
//    Output(_timeIter);
//
//    /*--- If the physical time convergence criteria has been met, terminate the simulation. ---*/
//
//    //if (StopCalc) break; Comes with physical time implementation
//
//    _timeIter++;
//
//  }
//
//}
//
// void EntropyClosureDriver::Run(){
//
//  /* ---- Solve dual problem ---- */
//  _optimizer->Solve(_moments);
//  _lambda = _optimizer->GetSolution();
//
//  /* ---- Reconstruct kineticDensity ---- */
//  _kineticDensity = _solver->ReconstructKineticDensity(_lambda);
//
//  /* ---- Solve Moment System ---- */
//  _solver->Solve(_lambda);
//  _moments = _solver->GetSolution();
//
//}
