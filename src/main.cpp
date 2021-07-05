// /*! @file: main.cpp
//  *  @brief Main method to call the KiT-RT solver suite
//  *  @author: J. Kusch, S. Schotthöfer, P. Stammer, J. Wolters, T. Xiao
//  *  @version: 0.1
//  */

// #include <Python.h>
// #define PY_ARRAY_UNIQUE_SYMBOL KITRT_ARRAY_API
// #include <mpi.h>
// #include <string>

// #include "common/config.h"
// #include "common/io.h"
// #include "solvers/solverbase.h"

// #include "toolboxes/datageneratorbase.h"

// #ifdef BUILD_GUI
// #include <QApplication>

// #include "mainwindow.h"
// #endif

// int main( int argc, char** argv ) {
// #ifdef BUILD_GUI
//     MPI_Init( &argc, &argv );
//     QApplication app( argc, argv );
//     MainWindow mw;
//     mw.show();
//     return app.exec();
// #else
//     MPI_Init( &argc, &argv );
//     wchar_t* program = Py_DecodeLocale( argv[0], NULL );
//     Py_SetProgramName( program );

//     std::string filename = ParseArguments( argc, argv );

//     // CD  Load Settings from File
//     Config* config = new Config( filename );

//     // Print input file and run info to file
//     PrintLogHeader( filename );

//     if( config->GetDataGeneratorMode() ) {
//         // Build Data generator
//         DataGeneratorBase* datagen = DataGeneratorBase::Create( config );
//         // Generate Data and export
//         datagen->ComputeTrainingData();
//     }
//     else {
//         // Build solver
//         SolverBase* solver = SolverBase::Create( config );

//         // Run solver and export
//         solver->Solve();
//         solver->PrintVolumeOutput();
//     }

//     MPI_Finalize();
//     return EXIT_SUCCESS;
// #endif
// }







/*! @file: main.cpp
 *  @brief Main method to call the KiT-RT solver suite
 *  @author: J. Kusch, S. Schotthöfer, P. Stammer, J. Wolters, T. Xiao
 *  @version: 0.1
 */

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL KITRT_ARRAY_API
#include <mpi.h>
#include <random>
#include <string>

#include "common/config.h"
#include "common/io.h"
#include "solvers/solverbase.h"

#include "toolboxes/datageneratorbase.h"

#ifdef BUILD_GUI
#include <QApplication>

#include "mainwindow.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <chrono>

#define PRINTF2(fp, ...) {printf(__VA_ARGS__);fprintf(fp,__VA_ARGS__);} // for mlmc_test


void regression(int, float *, float *, float &a, float &b);

float mlmc(int Lmin, int Lmax, int N0, float eps,
           void (*mlmc_l)(int, int, double *,Config*),Config* config,
           float alpha_0, float beta_0, float gamma_0,
           int *Nl, float *Cl) {

  double sums[7], suml[3][21];
  float  ml[21], Vl[21], NlCl[21], x[21], y[21],
         alpha, beta, gamma, sum, theta;
  int    dNl[21], L, converged;

  int    diag = 0;  // diagnostics, set to 0 for none 

  //
  // check input parameters
  //

  if (Lmin<2) {
    fprintf(stderr,"error: needs Lmin >= 2 \n");
    exit(1);
  }
  if (Lmax<Lmin) {
    fprintf(stderr,"error: needs Lmax >= Lmin \n");
    exit(1);
  }

  if (N0<=0 || eps<=0.0f) {
    fprintf(stderr,"error: needs N>0, eps>0 \n");
    exit(1);
  }

  //
  // initialisation
  //

  alpha = fmax(0.0f,alpha_0);
  beta  = fmax(0.0f,beta_0);
  gamma = fmax(0.0f,gamma_0);
  theta = 0.25f;             // MSE split between bias^2 and variance

  L = Lmin;
  converged = 0;

  for(int l=0; l<=Lmax; l++) {
    Nl[l]   = 0;
    Cl[l]   = powf(2.0f,(float)l*gamma);
    NlCl[l] = 0.0f;

    for(int n=0; n<3; n++) suml[n][l] = 0.0;
  }

  for(int l=0; l<=Lmin; l++) dNl[l] = N0;

  //
  // main loop
  //

  while (!converged) {

    //
    // update sample sums
    //
    
    for (int l=0; l<=L; l++) {
      if (diag) printf(" %d ",dNl[l]);

      if (dNl[l]>0) {
        mlmc_l(l,dNl[l],sums,config); // config added
        suml[0][l] += (float) dNl[l];
        suml[1][l] += sums[1];
        suml[2][l] += sums[2];
        NlCl[l]    += sums[0];  // sum total cost
      }
    }
    if (diag) printf(" \n");

    //
    // compute absolute average, variance and cost,
    // correct for possible under-sampling,
    // and set optimal number of new samples
    //



    sum = 0.0f;

    for (int l=0; l<=L; l++) {
      //std::cout<<"l ="<<"\t"<<l<<"\n";
      ml[l] = fabs(suml[1][l]/suml[0][l]);
      Vl[l] = fmaxf(suml[2][l]/suml[0][l] - ml[l]*ml[l], 0.0f);
      //std::cout<<"suml[2][l]/suml[0][l]"<<"\t"<<"ml[l]*ml[l]"<<"\n";
      //std::cout<<suml[2][l]/suml[0][l]<<"\t"<<ml[l]*ml[l]<<"\n";
      //std::cout<<"ml[l]"<<"\t"<<"Vl[l]"<<"\n";
      //std::cout<<ml[l]<<"\t"<<Vl[l]<<"\n";

      if (gamma_0 <= 0.0f) Cl[l] = NlCl[l] / suml[0][l];

      if (l>1) {
        ml[l] = fmaxf(ml[l],  0.5f*ml[l-1]/powf(2.0f,alpha));
        Vl[l] = fmaxf(Vl[l],  0.5f*Vl[l-1]/powf(2.0f,beta));
      }


      sum += sqrtf(Vl[l]*Cl[l]);
    }

    for (int l=0; l<=L; l++) {
      dNl[l] = ceilf( fmaxf( 0.0f, 
                       sqrtf(Vl[l]/Cl[l])*sum/((1.0f-theta)*eps*eps)
                     - suml[0][l] ) );
      
      //std::cout<<"l= "<<"\t"<<l<<"\n";
      //std::cout<<"Cl[l]"<<"\t"<<"Vl[l]"<<"\n";
      //std::cout<<Cl[l]<<"\t"<<Vl[l]<<"\n";
      //std::cout<<"S1"<<"\t"<<"suml[0][l]"<<"\t"<<"dNl[l]"<<"\n";
      //std::cout<<(Vl[l]/Cl[l])*sum/((1.0f-theta)*eps*eps)<<"\t"<<suml[0][l]<<"\t"<<dNl[l]<<"\n";
    }
 
    //
    // use linear regression to estimate alpha, beta, gamma if not given
    //

    if (alpha_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(ml[l]);
      }
      regression(L,x,y,alpha,sum);
      alpha = fmax(alpha,0.5f);
      
      if (diag) printf(" alpha = %f \n",alpha);
    }

    if (beta_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = - log2f(Vl[l]);
      }
      regression(L,x,y,beta,sum);
      beta = fmax(beta,0.5f);

      if (diag) printf(" beta = %f \n",beta);
    }

     if (gamma_0 <= 0.0f) {
      for (int l=1; l<=L; l++) {
        x[l-1] = l;
        y[l-1] = log2f(Cl[l]);
      }
      regression(L,x,y,gamma,sum);
      gamma = fmax(gamma,0.5f);

      if (diag) printf(" gamma = %f \n",gamma);
    }

    

    //
    // if (almost) converged, estimate remaining error and decide 
    // whether a new level is required
    //

    

    sum = 0.0;
      for (int l=0; l<=L; l++)
        sum += fmaxf(0.0f, (float)dNl[l]-0.01f*suml[0][l]);

        //std::cout<<sum<<"\n";

    if (sum==0) {
      if (diag) printf(" achieved variance target \n");

      converged = 1;
      float rem = ml[L] / (powf(2.0f,alpha)-1.0f);

      if (rem > sqrtf(theta)*eps) {
        if (L==Lmax)
          printf("*** failed to achieve weak convergence *** \n");
        else {
          converged = 0;
          L++;
          Vl[L] = Vl[L-1]/powf(2.0f,beta);
          Cl[L] = Cl[L-1]*powf(2.0f,gamma);

          if (diag) printf(" L = %d \n",L);

          sum = 0.0f;
          for (int l=0; l<=L; l++) sum += sqrtf(Vl[l]*Cl[l]);
          for (int l=0; l<=L; l++)
            dNl[l] = ceilf( fmaxf( 0.0f, 
                            sqrtf(Vl[l]/Cl[l])*sum/((1.0f-theta)*eps*eps)
                          - suml[0][l] ) );
        }
      }
    }
  }

  //
  // finally, evaluate multilevel estimator and set outputs
  //

  float P = 0.0f;
  for (int l=0; l<=L; l++) {   
    P    += suml[1][l]/suml[0][l];
    Nl[l] = suml[0][l];
    Cl[l] = NlCl[l] / Nl[l];
  }

  return P;
}

// linear regression routine

void regression(int N, float *x, float *y, float &a, float &b){

  float sum0=0.0f, sum1=0.0f, sum2=0.0f, sumy0=0.0f, sumy1=0.0f;

  for (int i=0; i<N; i++) {
    sum0  += 1.0f;
    sum1  += x[i];
    sum2  += x[i]*x[i];

    sumy0 += y[i];
    sumy1 += y[i]*x[i];
  }

  a = (sum0*sumy1 - sum1*sumy0) / (sum0*sum2 - sum1*sum1);
  b = (sum2*sumy0 - sum1*sumy1) / (sum0*sum2 - sum1*sum1);
}

void mlmc_test(void (*mlmc_l)(int, int, double *,Config*), Config* config,int M,int N,int L,
               int N0, float *Eps, int Lmin, int Lmax, FILE *fp) {

//
// first, convergence tests
//
  // current date/time based on current system
  time_t now = time(NULL);
  char *date = ctime(&now);
  int len = strlen(date);
  date[len-1] = ' ';

  PRINTF2(fp,"\n");
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"*** MLMC file version 0.9     produced by              ***\n");
  PRINTF2(fp,"*** C++ mlmc_test on %s         ***\n",date);
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"\n");
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"*** Convergence tests, kurtosis, telescoping sum check ***\n");
  PRINTF2(fp,"*** using N =%7d samples                           ***\n",N);
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"\n l   ave(Pf-Pc)    ave(Pf)   var(Pf-Pc)    var(Pf)");
  PRINTF2(fp,"    kurtosis     check        cost \n--------------------------");
  PRINTF2(fp,"-------------------------------------------------------------\n");

  double sums[7];
  float *cost = (float *)malloc((L+1)*sizeof(float));
  float *del1 = (float *)malloc((L+1)*sizeof(float));
  float *del2 = (float *)malloc((L+1)*sizeof(float));
  float *var1 = (float *)malloc((L+1)*sizeof(float));
  float *var2 = (float *)malloc((L+1)*sizeof(float));
  float *chk1 = (float *)malloc((L+1)*sizeof(float));
  float *kur1 = (float *)malloc((L+1)*sizeof(float));



  for (int l=0; l<=L; l++) {
    mlmc_l(l,N,sums,config); //config added

    for (int m=0; m<7; m++) sums[m] = sums[m]/N;

    if (M>0)
      cost[l] = powf((float)M,(float)l);
    else
      cost[l] = sums[0];
    del1[l] = sums[1];
    del2[l] = sums[5];
    var1[l] = fmax(sums[2]-sums[1]*sums[1], 1e-10);
    var2[l] = fmax(sums[6]-sums[5]*sums[5], 1e-10);

    kur1[l]  = (      sums[4]
                - 4.0*sums[3]*sums[1]
                + 6.0*sums[2]*sums[1]*sums[1]
                - 3.0*sums[1]*sums[1]*sums[1]*sums[1] )
             / (var1[l]*var1[l]);

    if (l==0)
      chk1[l] = 0.0f;
    else
      chk1[l] = sqrtf((float) N) * 
                fabsf(  del1[l]  +       del2[l-1]  -       del2[l] )
         / (3.0f*(sqrtf(var1[l]) + sqrtf(var2[l-1]) + sqrtf(var2[l])));

    PRINTF2(fp,"%2d  %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e \n",
    l,del1[l],del2[l],var1[l],var2[l],kur1[l],chk1[l],cost[l]);
  }

//
// print out a warning if kurtosis or consistency check looks bad
//

  if (kur1[L] > 100.0f) {
    PRINTF2(fp,"\n WARNING: kurtosis on finest level = %f \n",kur1[L]);
    PRINTF2(fp," indicates MLMC correction dominated by a few rare paths; \n");
    PRINTF2(fp," for information on the connection to variance of sample variances,\n");
    PRINTF2(fp," see http://mathworld.wolfram.com/SampleVarianceDistribution.html \n");
  }

  float max_chk = 0.0f;
  for (int l=0; l<=L; l++) max_chk = fmaxf(max_chk,chk1[l]);
  if (max_chk > 1.0f) {
    PRINTF2(fp,"\n WARNING: maximum consistency error = %f \n",max_chk);
    PRINTF2(fp," indicates identity E[Pf-Pc] = E[Pf] - E[Pc] not satisfied \n");
  }

//
// use linear regression to estimate alpha, beta, gamma
//

  float alpha, beta, gamma, foo;
  float *x = (float *)malloc(L*sizeof(float));
  float *y = (float *)malloc(L*sizeof(float));

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = - log2f(fabsf(del1[l]));
  } 
  regression(L,x,y,alpha,foo);

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = - log2f(var1[l]);
  } 
  regression(L,x,y,beta,foo);

  for (int l=1; l<=L; l++) {
    x[l-1] = l;
    y[l-1] = log2f(cost[l]);
  } 
  regression(L,x,y,gamma,foo);

  PRINTF2(fp,"\n******************************************************\n");
  PRINTF2(fp,"*** Linear regression estimates of MLMC parameters ***\n");
  PRINTF2(fp,"******************************************************\n");
  PRINTF2(fp,"\n alpha = %f  (exponent for MLMC weak convergence)\n",alpha);
  PRINTF2(fp," beta  = %f  (exponent for MLMC variance) \n",beta);
  PRINTF2(fp," gamma = %f  (exponent for MLMC cost) \n",gamma);

//
// second, mlmc complexity tests
//

  PRINTF2(fp,"\n");
  PRINTF2(fp,"***************************** \n");
  PRINTF2(fp,"*** MLMC complexity tests *** \n");
  PRINTF2(fp,"***************************** \n\n");
  PRINTF2(fp,"  eps       value   mlmc_cost   std_cost  savings     N_l \n");
  PRINTF2(fp,"--------------------------------------------------------- \n");
 
  int i=0;
  int   *Nl = (int *)malloc((Lmax+1)*sizeof(int));
  float *Cl = (float *)malloc((Lmax+1)*sizeof(float));

  while (Eps[i]>0) {
    float eps = Eps[i++];

    float P = mlmc(Lmin,Lmax,N0,eps,mlmc_l, config, alpha,beta,gamma, Nl,Cl);

    float std_cost = 0.0f, mlmc_cost = 0.0f, theta=0.25f;

    for (int l=0; l<=Lmax; l++) {
      if (Nl[l]>0) {
        // printf(" l, Cl, cost = %d  %f  %f \n",l,Cl[l],cost[l]);
        mlmc_cost += Nl[l]*Cl[l];
        if (l<=L) {
          std_cost = var2[l]*cost[l] / ((1.0f-theta)*eps*eps);
	}
        else
          std_cost = var2[L]*Cl[l] / ((1.0f-theta)*eps*eps);
      }
    }

    PRINTF2(fp,"%.4f  %.4e  %.3e  %.3e  %7.2f ",
	    eps, P, mlmc_cost, std_cost, std_cost/mlmc_cost);
    for (int l=0; Nl[l]>0; l++) PRINTF2(fp,"%9d",Nl[l]);
    PRINTF2(fp,"\n");
  }
  PRINTF2(fp,"\n");
}

void mlmc_test_100(void (*mlmc_l)(int, int, double *,Config*), Config* config,float val,
                   int N0, float *Eps, int Lmin, int Lmax, FILE *fp){

  // current date/time based on current system
  time_t now = time(NULL);
  char *date = ctime(&now);
  int len = strlen(date);
  date[len-1] = ' ';

  PRINTF2(fp,"\n");
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"*** MLMC file version 0.9     produced by              ***\n");
  PRINTF2(fp,"*** C++ mlmc_test on %s         ***\n",date);
  PRINTF2(fp,"**********************************************************\n");
  PRINTF2(fp,"\n");
  PRINTF2(fp,"***************************************** \n");
  PRINTF2(fp,"*** MLMC errors from 100 calculations *** \n");
  PRINTF2(fp,"***************************************** \n");

  if (isnan(val)) {
    PRINTF2(fp,"\n Exact value unknown \n");
  }
  else {
    PRINTF2(fp,"\n Exact value: %f \n",val);
  }

  int   i = 0;
  int   *Nl = (int *)malloc((Lmax+1)*sizeof(int));
  float *Cl = (float *)malloc((Lmax+1)*sizeof(float));

  while (Eps[i]>0) {
    float eps = Eps[i++];
    PRINTF2(fp,"\n eps = %.3e \n-----------------\n",eps); 

    for(int j=0; j<100; j++) {
      float P = mlmc(Lmin,Lmax,N0,eps,mlmc_l, config, 0.0f,0.0f,0.0f, Nl,Cl);
      PRINTF2(fp," %.5e ",P);
      if (j%5==4) PRINTF2(fp,"\n");
    }
  }
}

//level makes reference to the level of refinement ( quadrature), NnMCSamples number of iterations in the montecarlo, 
void MLquad_l(int level, int nMCSamples, double *sums, Config* config ) {    

  std::random_device rd;       // Will be used to obtain a seed for the random number engine
  std::mt19937 gen( rd() );    // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis( 1.0, 2.0 );
  double Pf; //Value of the functionu at level l
  double Pc; //Value of the function at level l-1
  double dP; //dPf - dPc
  std::vector<double> dosisf;
  std::vector<double> dosisc; 

  for (int k=0; k<7; k++) sums[k] = 0.0;    
  

  if(level==0){
      
    for( unsigned k = 0; k < nMCSamples; ++k ) {

      double density = dis( gen );
      //double density = 1.0 + k*(1.0/nMCSamples);//double density = dis( gen ); 

      config->SetQuadOrder(pow(2, 3+level));       

      SolverBase* solver = SolverBase::Create( config ); //function to be evaluated at level l 

      // solver->SetNewQN(pow(2, 3+level));  //replace by config->SetQuadOrder           
      solver->SetDensity( density );            
      solver->SolveQuietly(); //solver->Solve();

      dosisf = solver->GetDosis();       
      Pf = *std::max_element( dosisf.begin(), dosisf.end() );
      // Pf= std::accumulate( dosisf.begin(), dosisf.end(),0.0);

      delete solver;

      //std::cout<< Pf << "\t" << density << "\n";

      dP= Pf;

      sums[0] += level+1;     // add number of timesteps as cost
      sums[1] += dP;
      sums[2] += dP*dP;
      sums[3] += dP*dP*dP;
      sums[4] += dP*dP*dP*dP;
      sums[5] += Pf;
      sums[6] += Pf*Pf;

    }

  }else{       

    for( unsigned k = 0; k < nMCSamples; ++k ) {

      double density = dis( gen );
      //double density = 1.0 + k*(1.0/nMCSamples);//double density = dis( gen ); 

      config->SetQuadOrder(pow(2, 3+level));           

      SolverBase* solver = SolverBase::Create( config ); //function to be evaluated at level l 

      // solver->SetNewQN(pow(2, 3+level));             
      solver->SetDensity( density );            
      solver->SolveQuietly(); //solver->Solve();

      dosisf = solver->GetDosis();       
      Pf = *std::max_element( dosisf.begin(), dosisf.end() );
      // Pf= std::accumulate( dosisf.begin(), dosisf.end(),0.0);

      delete solver;

      config->SetQuadOrder(pow(2, 2+level));

      SolverBase* solverc = SolverBase::Create( config ); //function to be evaluated at level l 

      // solverc->SetNewQN(pow(2, 2+level));             
      solverc->SetDensity( density );            
      solverc->SolveQuietly(); //solverc->Solve();

      dosisc = solverc->GetDosis();       
      Pc = *std::max_element( dosisc.begin(), dosisc.end() );
      // Pc= std::accumulate( dosisc.begin(), dosisc.end(),0.0); //nicht einfach dosis sonder dosis.dx (/X)
      

      delete solverc;

      dP= Pf - Pc;

      //std::cout<< Pf << "\t"<< Pc << "\t"<< dP << "\t" << density<< "\n";

      

      sums[0] += level+1;     // add number of timesteps as cost
      sums[1] += dP;
      sums[2] += dP*dP;
      sums[3] += dP*dP*dP;
      sums[4] += dP*dP*dP*dP;
      sums[5] += Pf;
      sums[6] += Pf*Pf;            
    } 
      
  } 
    
}

void MLmesh_l(int level, int nMCSamples, double *sums, Config* config ) {  

  std::string listofmeshes [] = {"meshes/1DMesh200.su2","meshes/1DMesh400.su2","meshes/1DMesh800.su2","meshes/1DMesh1600.su2","meshes/1DMesh3200.su2",
  "meshes/1DMesh6400.su2","meshes/1DMesh12800.su2","meshes/1DMesh25600.su2"};

  std::random_device rd;       // Will be used to obtain a seed for the random number engine
  std::mt19937 gen( rd() );    // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis( 1.0, 2.0 );
  double Pf; //Value of the functionu at level l
  double Pc; //Value of the function at level l-1
  double dP; //dPf - dPc
  std::vector<double> dosisf;
  std::vector<double> dosisc; 

  for (int k=0; k<7; k++) sums[k] = 0.0;    
  

  if(level==0){
      
    for( unsigned k = 0; k < nMCSamples; ++k ) {

      double density = dis( gen );
      //double density = 1.0 + k*(1.0/nMCSamples);//double density = dis( gen );

      config->SetMeshFile(listofmeshes[level]);            

      SolverBase* solver = SolverBase::Create( config ); //function to be evaluated at level l 

      // solver->SetNewMesh(listofmeshes[level]);             
      solver->SetDensity( density );            
      solver->SolveQuietly(); //solver->Solve();

      dosisf = solver->GetDosis();       
      Pf = *std::max_element( dosisf.begin(), dosisf.end() );
      // Pf= std::accumulate( dosisf.begin(), dosisf.end(),0.0);

      delete solver;

      //std::cout<< Pf << "\t" << density << "\n";

      dP= Pf;

      sums[0] += level+1;     // add number of timesteps as cost
      sums[1] += dP;
      sums[2] += dP*dP;
      sums[3] += dP*dP*dP;
      sums[4] += dP*dP*dP*dP;
      sums[5] += Pf;
      sums[6] += Pf*Pf;

    }

  }else{       

    for( unsigned k = 0; k < nMCSamples; ++k ) {

      double density = dis( gen );
      //double density = 1.0 + k*(1.0/nMCSamples);//double density = dis( gen );  

      config->SetMeshFile(listofmeshes[level]);          

      SolverBase* solver = SolverBase::Create( config ); //function to be evaluated at level l 

      // solver->SetNewMesh(listofmeshes[level]);             
      solver->SetDensity( density );            
      solver->SolveQuietly(); //solver->Solve();

      dosisf = solver->GetDosis();       
      Pf = *std::max_element( dosisf.begin(), dosisf.end() );
      // Pf= std::accumulate( dosisf.begin(), dosisf.end(),0.0);

      delete solver;

      config->SetMeshFile(listofmeshes[level-1]);

      SolverBase* solverc = SolverBase::Create( config ); //function to be evaluated at level l 

      // solverc->SetNewMesh(listofmeshes[level-1]);             
      solverc->SetDensity( density );            
      solverc->SolveQuietly(); //solverc->Solve();

      dosisc = solverc->GetDosis();       
      Pc = *std::max_element( dosisc.begin(), dosisc.end() );
      // Pc= std::accumulate( dosisc.begin(), dosisc.end(),0.0); //nicht einfach dosis sonder dosis.dx (/X)
      

      delete solverc;

      dP= Pf - Pc;

      //std::cout<< Pf << "\t"<< Pc << "\t"<< dP << "\t" << density<< "\n";

      

      sums[0] += level+1;     // add number of timesteps as cost
      sums[1] += dP;
      sums[2] += dP*dP;
      sums[3] += dP*dP*dP;
      sums[4] += dP*dP*dP*dP;
      sums[5] += Pf;
      sums[6] += Pf*Pf;            
    } 
      
  } 
    
}

int main( int argc, char** argv ) {
#ifdef BUILD_GUI
    MPI_Init( &argc, &argv );
    QApplication app( argc, argv );
    MainWindow mw;
    mw.show();
    return app.exec();
#else

    MPI_Init( &argc, &argv );
    wchar_t* program = Py_DecodeLocale( argv[0], NULL );
    Py_SetProgramName( program );
    std::string filename = ParseArguments( argc, argv );
    // CD  Load Settings from File
    Config* config = new Config( filename );
    // Print input file and run info to file
    //PrintLogHeader( filename );
    if( config->GetDataGeneratorMode() ) {
        // Build Data generator
        DataGeneratorBase* datagen = DataGeneratorBase::Create( config );
        // Generate Data and export
        datagen->ComputeTrainingData();
    }
    else {

    // std::cout<< config->GetNQuadPoints()<<"\n";
    // std::string listofmeshes [6] = {"meshes/1DMesh2000.su2","meshes/1DMesh2500.su2","meshes/1DMesh3000.su2","meshes/1DMesh3200.su2","meshes/1DMesh1000.su2","meshes/1DMesh1250.su2"};

    // config->SetNQuadPoints(10);
    // config->SetQuadOrder(20);
    // config->SetMeshFile(listofmeshes[1]); 

    // std::cout<< config->GetNQuadPoints()<<"\n";

    config->SetQuadOrder(20);    
    SolverBase* solver = SolverBase::Create( config);
    // solver->SetDensity( 1.5); 
    solver->SetDensity2( 1, 0.5 ); 
    solver->SetDensity3( 1, 0.25, 0.75 ); 
    solver->SolveQuietly(); //solver->Solve();
    std::vector<double> data = solver->GetDensityVector();


    int wallstart=0;
    int wallend=0;

    for(int i=0;i< data.size();i++){
      std::cout << data[i] << "\t";
      if(i>0){
        if(data[i]>0 and data[i-1]==0) {wallstart=i;}
        if(data[i-1]>0 and data[i]==0) {wallend=i;}
      }
    }

    std::cout << "\n";

    std::cout << wallstart << "\n";
    std::cout << wallend << "\n";
    std::cout << data.size() << "\n";

    delete solver;
   
    // std::string Meshname = "meshes/1DMesh1250.su2";
    // solver->SetNewMesh(Meshname);
    // std::cout<<solver->GetNQ()<<"\n";
    // solver->SetNewQN(20);
    // std::cout<<solver->GetNQ()<<"\n";
    // solver->SolveQuietly();
    // std::vector<double> dosis = solver->GetDosis();
    // std::cout<< *std::max_element( dosis.begin(), dosis.end() )<<std::endl;

    // std::cout<< config->GetNQuadPoints()<<"\n";

    // delete solver;


      
    // double sums[7];        
    // for(int n=0; n<7; n++) sums[n] = 0.0; 

    // auto start = std::chrono::high_resolution_clock::now();       
    // MLmesh_l(5,5, sums, config );
    // auto stop = std::chrono::high_resolution_clock::now();

    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    // std::cout<< "Time taken by function: " << duration.count() << " microseconds" << endl;
      

        
        

      /*
      std::vector<double> dosisf;
      double density;
      double Pf;
      SolverBase* solver = SolverBase::Create( config );

      for( unsigned k = 0; k < 5; ++k ) {
          double density = 1.2 + k*0.1;       

          solver->SetNewQN(18);             
          solver->SetDensity( density );            
          solver->SolveQuietly(); //solver->Solve();
          dosisf = solver->GetDosis()       
          Pf = *std::max_element( dosisf.begin(), dosisf.end() );

          std::cout<< Pf << "\t" << solver->GetDensity()<< "\n"; 
      }

      delete solver;

      */    
      //for(int n=0; n<7; n++) std::cout<<sums[n]<<"\n";      

      // FILE *fp;
      // fp = fopen("KIT-RT mlmc1","w");
      
      
      // float Eps[11];
      // float Eps2[] = { 50, 10, 1, 0.0 };        

      // memcpy(Eps,Eps2,sizeof(Eps2));        

      // int M  = 3;     // refinement cost factor
      // int N0 = 15;    // initial samples on each level
      // int Lmin = 2;   // minimum refinement level
      // int Lmax = 6;   // maximum refinement level
      // int N  = 5;    // samples for convergence tests
      // int L  = 6;     // levels for convergence tests
      // int NL = 5;    //Number of sampels for level L (for mlcl())      
      // float Cl = 3.0; //Cost of level L (for mlmc())
      // float eps = 50.0;  //epsilon (for mlmc())              

      // // // //double P = mlmc(Lmin,Lmax,N0,eps,MLmesh_l,config,0.0,0.0,0.0,&NL,&Cl);
      // // // //std::cout<<P;

      // // // mlmc_test(MLquad_l,config,M,N,L,N0,Eps,Lmin,Lmax,fp);
      // mlmc_test(MLmesh_l,config,M,N,L,N0,Eps,Lmin,Lmax,fp);  
      

    }
    MPI_Finalize();
    return EXIT_SUCCESS;
#endif
}
