//
// Created by chinsp on 29/10/21.
//

#include "../../include/solvers/refinesnsolver.h"
#include "common/config.h"
#include "common/io.h"
#include "common/mesh.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/problembase.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/textprocessingtoolbox.h"
#include "quadratures/quadraturebase.h"
#include "common/globalconstants.h"

// externals
#include "spdlog/spdlog.h"
#include <iostream>
#include <mpi.h>
#include <time.h>
#include <stdio.h>

using namespace std;

/*
    Idee:   Define 2 quadrature objects
            + Vector that gives refinement points for each coarse point
            + external source for fine quad and coarse quad

    Refinement: choose which direction has to be refined
                interchange direction with refined directions ( bool vector -> 0 to 1 )
                (also refine neighbors of refined triangles)

    Choose directions to refine:   -> depending on external source
    @todo   if ext source is same for all energies do calculation before loop
            if ext source is isotropic not necessary to do all the calculations -> simplification
*/

RefineSNSolver::RefineSNSolver( Config* settings ) : SNSolver ( settings ) {

    // information for quadratures
    _nqC = _quadrature->GetNq();
    _quadPointsC = _quadrature->GetPoints();
    _weightsC = _quadrature->GetWeights();
    _neighC = _quadrature->GetNeighbours();

    _nqF = _quadrature->GetNqRefined();
    _quadPointsF = _quadrature->GetPointsRefined();
    _weightsF = _quadrature->GetWeightsRefined();
    // _neighF = _quadrature->GetNeighboursRefined();

    _quadC2F = _quadrature->GetRefineVector();

    // setup external source ( in fine and coarse grid )
    _QextC = _Q;
    _settings->SetNQuadPoints(_nqF);
    _QextF = _problem->GetExternalSource( _energies );
    _solF = _problem->SetupIC();
    _settings->SetNQuadPoints(_nqC);

    // setup helper vectors
    _quadIDs = VectorU(_nqF, 0);
    _refine = VectorU(_nqC, 0);
    _refineIter = VectorU(_nqC, 0);

    _beta = _settings->GetBetaAS();
    _sigmaAS = _settings->GetSigmaAS();
}

void RefineSNSolver::Solve(){

    // --- Pre Solver Processing ---

    start = clock(); // time tic for complete run

    PrepareVolumeOutput();
    DrawPreSolverOutput();
    SolverPreprocessing();

    Times = fopen((_settings->GetLogDir() + "/TimesCompareRefine.txt").c_str(), "a+");
    radFlux = fopen((_settings->GetLogDir() + "/radFluxRefine.txt").c_str(), "a+");
    Refine = fopen((_settings->GetLogDir() + "/RefineVector.txt").c_str(), "a+");
    Info = fopen((_settings->GetLogDir() + "/Infos.txt").c_str(), "a+");
    ASFile = fopen((_settings->GetLogDir() + "/ASKernel.txt").c_str(), "a+");

    // --- Solver start ---
    for ( unsigned idx_energy = 0; idx_energy < _nEnergies; idx_energy++ ){

        startiter = clock(); // time tic for one iteration

        // --- Preprocessing --> Setup new Quadrature; Adapt _sol; Adapt _Qext; ---
        IterPreprocessing( idx_energy );

        // --- Update ---
        t = clock();
        FluxUpdate();
        FVMUpdate( idx_energy );

        if ( _sigmaAS != 0.0)
            AddArtificialScattering( idx_energy ); // Add Artificial Scattering Term

        t4 = clock() - t;

        // --- Postprocessing --> Calculate Radient Flux ---
        IterPostprocessing( idx_energy );

        // --- File Outputs ---
        clock_t enditer = clock() - startiter;
        fprintf(Times, "Step 1 -- Check if Refinement needed:\t %f \n", t1/(float)CLOCKS_PER_SEC );
        fprintf(Times, "Step 2 -- Refinement:\t \t \t %f \n", t2/(float)CLOCKS_PER_SEC );
        fprintf(Times, "Step 2.1 -- Setup Quadrature:\t \t %f \n", t21/(float)CLOCKS_PER_SEC );
        fprintf(Times, "Step 2.2 -- Setup Interpolationmatrix:\t %f \n", t22/(float)CLOCKS_PER_SEC );
        fprintf(Times, "Step 2.3 -- Interpolate next solution:\t %f \n", t23/(float)CLOCKS_PER_SEC );
        fprintf(Times, "Step 2.4 -- Compute ScatteringKernel:\t %f \n", t24/(float)CLOCKS_PER_SEC );
        fprintf(Times, "Step 3 -- Determine next ext Source:\t %f \n", t3/(float)CLOCKS_PER_SEC );
        fprintf(Times, "Step 4 -- SN-Update: \t \t \t %f \n", t4/(float)CLOCKS_PER_SEC );
        fprintf(Times, "\n Complete Iteration %u: \t %f sec. \n \n", idx_energy, enditer/(float)CLOCKS_PER_SEC);

        fprintf(radFlux, "Iteration %u \n", idx_energy);
        for (unsigned i = 0; i< _nCells; i++)
            if ( _fluxNew[i] > 10e-5)
                fprintf(radFlux, "\t Flux[%u] = %f \n", i, _fluxNew[i] );

        // --- Solver Output ---
        WriteVolumeOutput( idx_energy );
        WriteScalarOutput( idx_energy );
        PrintScreenOutput( idx_energy );
        PrintHistoryOutput( idx_energy );
        PrintVolumeOutput( idx_energy );

    }   // end of iteration

    // --- Post Solver Processing  ---
    DrawPostSolverOutput();

    clock_t end = clock() - start;
    fprintf(Times, "Complete Run: %f \n \n", end/(float)CLOCKS_PER_SEC );
    fclose(Times);
}

void RefineSNSolver::IterPreprocessing( unsigned idx_energy ){

    // save Information from last iteration
    if ( idx_energy == 0){ // in first Iteration always use values of refinement bc of IC is fine
        _nqOld = _nqF; _quadPointsOld =_quadPointsF; _quadIDsOld = _quadIDs; _refineIterOld = _refineIter;
    }
    else{ _nqOld = _nq; _quadPointsOld = _quadPoints; _quadIDsOld = _quadIDs; _refineIterOld = _refineIter; }


    // --- Step 1 --- Determine if and where Refinement is necessary (depending on external source)

    t = clock(); // tic for Step 1
    _refine = VectorU( _nqC, 0); // reset Vector

    // loop over COARSE external Source ( if _Qext[idx_energy][..][idx_quad] != 0 have to refine in direction idx_quad )
    for (unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++){
        if( _QextC.size() == 1u ) {                     // constant source for all energies
            if( _QextC[0][idx_cell].size() == 1u )          // isotropic source
                if( _QextC[0][idx_cell][0] == 0.0) continue;    // Isotropic source is 0 -> nothing to do
                else { _refine = VectorU(_nqC,1); break; }      // Isotropic Source is not 0 -> Refinement
            else{                                          // anisotropic source
                for (unsigned idx_quad = 0; idx_quad < _nqC; idx_quad++){
                    if( _QextC[0][idx_cell][idx_quad] == 0.0)
                        continue;                               // Source[idx_quad] = 0 -> not refine this point
                    else {
                        _refine[idx_quad] = 1;                  // Source[idx_quad] != 0 -> refine this point
                        for ( unsigned id = 0; id < _neighC[idx_quad].size(); id ++)
                            _refine[ _neighC[idx_quad][id]] = 1;     // also refine neighbours of this point
                    }
                }
            }
        }
        else {                                          // Same as above but for different sources for different energies
            if( _QextC[idx_energy][idx_cell].size() == 1u )
                if( _QextC[idx_energy][idx_cell][0] == 0.0) continue;
                else { _refine = VectorU(_nqC,1); break; }
            else{
                for (unsigned idx_quad = 0; idx_quad < _nqC; idx_quad++){
                    if( _QextC[idx_energy][idx_cell][idx_quad] == 0.0) continue;
                    else {
                        _refine[idx_quad] = 1;
                        for ( unsigned id = 0; id < _neighC[idx_quad].size(); id ++)
                            _refine[ _neighC[idx_quad][id]] = 1;     // also refine neighbours of this point
                    }
                }
            }
        }
    }
    // --> _refine in {0,1} -> need refinement in this iteration or not

    // Determine if we have to change new quadrature (depending on last Iterations)
    fprintf(Refine, "\n Iteration %u \n", idx_energy);
    fprintf(Refine, "\t Point \t Refine \t Change \t RefineIter \n ");
    for (unsigned id = 0; id < _nqC; id++){
        fprintf(Refine, "\t %u \t %u", id, _refine[id]);
        if ( _refine[id] == 0 ){                            // NO refinement needed for quadpoint ID
            if ( _refineIter[id] >= _settings->GetRefineIter() )
            { _refineIter[id] = 0; _refine[id] = 1; }      // reset if enough iterations -> have to change quadratre
            else if ( _refineIter[id] == 0 ) {}                  // have nothing to change (was not refined, have not to be refined)
            else { _refineIter[id] += 1; }                      // was refined, dont need to be refined again -> nothing to change
        }
        else if (_refine[id] == 1 ){                    // NEED refinement for quadpoint ID
            if ( _refineIter[id] != 0) { _refine[id] = 0; }     // was also refined in last Iter -> nothing to change
            _refineIter[id] = 1;                                // RESET counter to 1
        }
        else {ErrorMessages::Error( "Refinement went wrong!", CURRENT_FUNCTION ); break;}

        fprintf(Refine, "\t\t %u \t\t %u \n", _refine[id], _refineIter[id]);
    }

    if ( idx_energy == 1 ) // in step 1 compute complete new quadrature depending on where to refine, bc in step 0 we used complete fine quad
        _refine = 1;

    t1 = clock() - t; // toc for step 1
    // --> _refine in  {0,1} -> do we have to change something


    // --- Step 2 --- Refinement
    tt = clock();
    if ( idx_energy == 0 ){ // in first iteration always fine quadrature
        t21 = 0; t22 = 0; t23 = 0; t24 = 0;
        _nq = _nqF;
        _quadPoints = _quadPointsF;
        _weights = _weightsF;
        _quadIDs = VectorU( _nq );
        for (unsigned id = 0; id < _nq; id++){ _quadIDs[ id ] = id + _nqC;  }
        _sol = _solF;
        _solNew = _sol; _psiDx = _sol; _psiDy = _sol;
        _scatteringKernel = GetScatteringKernel(_nq, _weights);
        fprintf(Info, "Iteration: %u \n\t Use Complete Fine Quadrature \n\t Nr of Quad Points: \t %u \n", idx_energy, _nq);
    }
    else if ( idx_energy != 0 && _refine == 0 ){ // NO need to change quadrature in this Iter -> go on with standard SN
        t21 = 0; t22 = 0; t23 = 0; t24 = 0; // set toc values
        _quadIDs = _quadIDsOld;         // quadIDs are the same as in last Iteration

        fprintf(Info, "Iteration: %u \n\t Nothing to Change from last Iteration\n\t Nr of Quad Points: \t %u \n", idx_energy, _nqOld);
    }
    else{
        // --- Step 2.1 --- Refining the quadrature
        t = clock();

        if ( isZero( _refineIter ) ){   // no direction is refined -> use coarse quadrature
            _nq = _nqC;
            _quadPoints = _quadPointsC;
            _weights = _weightsC;
            _quadIDs = VectorU( _nq );
            for (unsigned id = 0; id < _nq; id++){ _quadIDs[ id ] = id;  }
            fprintf(Info, "Iteration: %u \n\t Use Complete Coarse Quadrature\n\t Nr of Quad points: \t %u \n", idx_energy, _nq);
        }
        else if ( nonZeros( _refineIter ) == _nqC ){ // all directions are refined -> use fine quadrature
            _nq = _nqF;
            _quadPoints = _quadPointsF;
            _weights = _weightsF;
            _quadIDs = VectorU( _nq );
            for (unsigned id = 0; id < _nq; id++){ _quadIDs[ id ] = id + _nqC;  }
            fprintf(Info, "Iteration: %u \n\t Use Complete Fine Quadrature\n\t Nr of Quad Points: \t %u \n", idx_energy, _nq);
        }
        else{   // just some directions have to be coarsened or refined
            _nq = 0;
            // reset(_quadPoints);
            // reset(_weights);
            // reset(_quadIDs);
            _quadPoints = VectorVector( _nqF, Vector( 3 )); // reset calculation vectors : _nqF is maximum of possible quadPoints
            _weights = Vector( _nqF );
            _quadIDs = VectorU( _nqF );

            for (unsigned idx_quad = 0; idx_quad < _nqC; idx_quad++){
                if ( _refineIter[ idx_quad ] == 0 ){         // use coarse quadrature point
                    _quadPoints[ _nq ] = _quadPointsC[ idx_quad ];
                    _weights[ _nq ] = _weightsC[ idx_quad ];
                    _quadIDs[ _nq ] = idx_quad;
                    _nq += 1;
                }
                else{                                       // use the corresponding points from fine quadrature
                    for ( unsigned id = 0; id < _quadC2F[idx_quad].size(); id++){
                        _quadPoints[ _nq + id ] = _quadPointsF[ _quadC2F[idx_quad][id] ];
                        _weights[ _nq + id ] = _weightsF[ _quadC2F[idx_quad][id] ];
                        _quadIDs[ _nq + id ] = _nqC + idx_quad + id;
                    }
                    _nq += _quadC2F[idx_quad].size();
                }
            }
            // truncate vectors to _nq
            _quadPoints.resize( _nq );
            _weights.resize( _nq );
            _quadIDs.resize( _nq );

            fprintf( Info, "Iteration %u \n\t Nr of Quad Points: \t %u \n\t Nr of Refined Points: \t %u \n\t Refined Points: ", idx_energy, _nq, (unsigned) nonZeros(_refineIter));
            for( unsigned id = 0; id < _nqC; id++)
                if ( _refineIter[id] != 0)
                    fprintf(Info, "\t [%u]", id);
            fprintf(Info , "\n");
        }
        t21 = clock() - t;

        // --- Step 2.2 and 2.3 --- Interpolation according to Jarrell 2010
        InterpolateFlux();

        // --- Step 2.2 --- Computation of Interpolation Matrix
        t = clock();

        _solInter = Matrix(_nq, _nqOld, 0.0);
        for (unsigned iterNew = 0; iterNew < _nq; iterNew ++){
            unsigned check = 0;
            for (unsigned iterOld = 0; iterOld < _nqOld; iterOld ++){ // these are the points that havent changed in new quadrature
                if ( _quadIDs[ iterNew ] == _quadIDsOld[ iterOld ] ){
                    _solInter( iterNew, iterOld ) = 1;
                    check = 1;
                    break;
                }
            }
            if ( check == 1) continue;
            else{
                CalculateInterpolation( iterNew ); // Interpolation between three nearest values in old quadrature
                //CalculateInterpolation2( iterNew );
            }
        }
        t22 = clock() - t;

        // --- Step 2.3 --- Interpolation of next _sol variable
        t = clock();
        _solNew = VectorVector( _nCells, Vector( _nq ));
        for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
            _solNew[idx_cell] = _solInter * _sol[idx_cell];
        }
        _sol = _solNew;
        t23 = clock() - t;

        // --- Step 2.4 --- Computation of new scattering Kernel
        t = clock();
        _scatteringKernel = GetScatteringKernel( _nq, _weights );
        t24 = clock() - t;

        // set slope variables ->  only sizes are relevant ( are computetd in FluxUpdate )
        _psiDx = _sol;
        _psiDy = _sol;
    }
    t2 = clock() - tt;
    // --- end of if else loop ---

    // --- Step 3 --- Set external Source dependent on quadrature points

    t = clock();
    for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
        if( _Q.size() == 1u ) {                   // constant source for all energies
            if( _Q[0][idx_cell].size() == 1u ) { }  // isotropic source
            else{
                _Q[0][idx_cell] = Vector( _nq );
                for ( unsigned  id = 0; id < _nq; id++){
                    if ( _quadIDs[ id ] < _nqC ){    // i.e. point of coarse quadrature
                        _Q[0][idx_cell][ id ] = _QextC[0][idx_cell][ _quadIDs[id] ];
                    }
                    else {
                        _Q[0][idx_cell][ id ] = _QextF[0][idx_cell][ _quadIDs[id] - _nqC ];
                    }
                }
            }
        }
        else {
            if( _Q[0][idx_cell].size() == 1u ) {}   // isotropic source
            else{
                _Q[idx_energy][idx_cell] = Vector( _nq );
                for ( unsigned  id = 0; id < _nq; id++){
                    if ( _quadIDs[ id ] < _nqC ) {// i.e. point of coarse quadrature
                        _Q[idx_energy][idx_cell][ id ] = _QextC[idx_energy][idx_cell][ _quadIDs[id] ];
                    }
                    else {
                        _Q[idx_energy][idx_cell][ id ] = _QextF[idx_energy][idx_cell][ _quadIDs[id] - _nqC ];
                    }
                }
            }
        }
    }
    t3 = clock() - t;
}


void RefineSNSolver::CalculateInterpolation( unsigned quadIter ){
    std::vector<double> distance1( _nqOld );
    std::vector<unsigned> ids( _nqOld );

    for( unsigned id = 0; id < _nqOld; id++){
        distance1[ id ] = norm( _quadPoints[ quadIter ] - _quadPointsOld[ id ]); // calculate distance to every point in last iteration
    }
    std::iota( ids.begin(), ids.end(), 0);
    std::sort( ids.begin(), ids.end(), [&](int i,int j){return distance1[i] < distance1[j];} ); // sort distance vector
    double ssum = 0;
    if ( distance1[ ids [0] ] < 10e-10 ){   // if new point and old point are nearly the same
        _solInter( quadIter, ids [ 0 ]) = 1;
    }
    else{
        for ( unsigned it = 0 ; it < 3 ; it++){
            ssum +=  1 /  pow ( distance1[ ids[ it ] ] , 2 );    // sum over three nearest values
        }
        for ( unsigned it = 0; it < 3; it++){   // setup Interpolation Matrix
            _solInter( quadIter, ids[ it ] ) = ( 1 /  pow ( distance1[ ids[ it ] ] , 2 ) );
            _solInter( quadIter, ids[ it ] ) /= ssum;
        }
    }
}

void RefineSNSolver::InterpolateFlux(){

    for ( unsigned id = 0; id < _nqC; id ++){
        if ( _refineIter[ id ] == 0 && _refineIterOld[ id ] == 0) // was coarse, is coarse -> nothing to interpolate
        {}
        else if ( _refineIter[ id ] != 0 && _refineIterOld[ id ] != 0 ) // was fine, is fine -> nothing to interpolate
        {}
        else if ( _refineIter[ id ] == 0 && _refineIterOld[ id ] != 0 ) // was fine, is coarse -> fine to coarse mapping
        {} // --> need information from master triangle of coarse triangles
        else if ( _refineIter[ id ] != 0 && _refineIterOld[ id ] == 0 ) // was coarse, is fine -> coarse to fine mapping
        {}
    }
}

Matrix RefineSNSolver::GetScatteringKernel( unsigned nq, Vector w) {
    // unsigned nq = _quad->GetNq();
    // auto w      = _quad->GetWeights();
    Matrix kernel( nq, nq );
    for( unsigned i = 0; i < nq; ++i )
        for( unsigned j = 0; j < nq; ++j ) kernel( i, j ) = w[j] / ( 4 * M_PI );

    // scale kernel to ensure mass conservation
    double tmp;
    for( unsigned i = 0; i < nq; ++i ) {
        tmp = 0.0;
        for( unsigned j = 0; j < nq; ++j ) {
            tmp += kernel( i, j );
        }
        for( unsigned j = 0; j < nq; ++j ) {
            kernel( i, j ) /= tmp;
        }
    }
    return kernel;
}

void RefineSNSolver::AddArtificialScattering( unsigned iter ) {

    _ASKernel = Matrix (_nq, _nq, 0.0);

    fprintf ( ASFile, "\n Iteration %u \n \n ", iter);

    // Calculate Artificial Scattering Kernel
    for ( unsigned id = 0; id < _nq; id ++){
        double sum = 0;
        for ( unsigned jd = 0; jd < _nq; jd ++){
            double mu = 0;
            for ( unsigned i = 0; i < _quadPoints[0].size(); i++){
                mu += _quadPoints[id][i] * _quadPoints[jd][i];
            }
            double expTerm = ( - ( ( mu - 1 ) * ( mu - 1 ) ) * _nq * _nq ) / ( _beta * _beta );
            _ASKernel( id, jd ) = _nq / _beta * exp( expTerm );
            sum += _ASKernel( id, jd );
        }
        for ( unsigned jd = 0; jd < _nq; jd++){
            _ASKernel( id , jd ) /= sum;        // Normalize
            if ( _ASKernel ( id , jd ) > 10e-5)
                fprintf(ASFile, "\t ASKernel( %u, %u) = %f \n", id, jd, _ASKernel(id, jd));
        }
        fprintf(ASFile, "sum of row %u: %f \n", id, sum);
    }

    // Add Artificiel Scattering Terms to actual solution
    for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
        _solNew[idx_cell] += _dE * _sigmaAS * ( _ASKernel * _sol[idx_cell] - _sol[idx_cell] );
    }
}