//
// Created by chinsp on 29/10/21.
//

#include "../../include/solvers/firstcollisioncsdsnsolver.h"

#include "common/config.h"
#include "common/io.h"
#include "common/mesh.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/textprocessingtoolbox.h"
#include "entropies/entropybase.h"
#include "fluxes/numericalflux.h"
#include "optimizers/optimizerbase.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/icru.h"
#include "problems/problembase.h"
#include "quadratures/quadraturebase.h"


// externals
#include "spdlog/spdlog.h"
#include <iostream>
#include <mpi.h>
#include <time.h>
#include <stdio.h>


FirstCollisionCSDSNSolver::FirstCollisionCSDSNSolver( Config * _settings ) : Solver ( _settings ){

    _FPMethod = 2;
    _RT = false;

    /////// Initialize dose and density for CSD
    _dose = std::vector<double>( _settings->GetNCells(), 0.0);

    // recompute energy (pseudo-time) steps
    _energies = Vector( _nEnergies, 0.0 );
    double _energyMin = 1e-4 * 0.511;
    _energyMax = _settings->GetMaxEnergyCSD();

    _dE = ComputeTimeStep( _settings->GetCFL() );
    _nEnergies = unsigned ( ( _energyMax - _energyMin ) / _dE );
    _energies.resize( _nEnergies );
    for ( unsigned id = 0; id < _nEnergies; id ++){
        _energies[ id ] = _energyMin + ( _energyMax - _energyMin ) / ( _nEnergies - 1 ) * id;
    }

    /////// Variables for coarse SN Setting

    switch (_settings->GetQuadName()){
        case QUAD_GaussLegendre1D: break;
        default: ErrorMessages::Error("Please choose Gauss Legendre 1D quadrature!", "FirstCollision CSD Solver");
    }
    _nqC = _quadrature->GetNq();
    _quadPointsC = _quadrature->GetPoints();
    _weightsC    = _quadrature->GetWeights();

    ////// Variables for fine SN Setting

    _quadratureFine = QuadratureBase::Create( _settings->GetQuadName(), _settings->GetQuadOrderFine() );
    _nqF = _quadratureFine->GetNq();
    _quadPointsF = _quadratureFine->GetPoints();
    _weightsF = _quadratureFine->GetWeights();


    // Setup Variables for Uncollided Calculation
    if( _settings->GetIsLocalRefine() ){
        Refinement();
        std::cout << "Refinemnt Check" << std::endl;
    }
    else if( _settings->GetIsArtificialScattering() ){
        _nq = _nqC;
        _quadPoints = _quadPointsC;
        _weights = _weightsC;

        _sigmaAS = _settings->GetSigmaAS();
        _betaAS = _settings->GetBetaAS();
        GetASTransportCoefficients(); // ->Get _xiAS;
    }
    else {
        _nq = _nqF;
        _quadPoints = _quadPointsF;
        _weights = _weightsF;
    }

    // Setup Interpolation Matrix (Checked)
    /*
    FILE * PTS = fopen((_settings->GetLogDir() + "/CellPts").c_str(), "a+");
    for( unsigned i = 0; i <_nCells; i++)
        fprintf(PTS, "%f\t%f\n", _mesh->GetCellMidPoints()[i][0], _mesh->GetCellMidPoints()[i][1]);
    */
    //if (  _settings->GetIsLocalRefine() == false ){
    if( _settings->GetIsArtificialScattering() ){ // here both quadratures are the same
        _MatInterp = Matrix(_nqC, _nqC, 0.0);
        for ( unsigned id = 0; id < _nqC; id++){ _MatInterp(id, id) = 1.0; }
    }
    else { // Interpolate between the two nearest mu  Values
        _MatInterp = Matrix(_nqC, _nq );
        for (unsigned id = 0; id < _nqC; id ++){
            Vector dist1 = Vector( _nq, 0.0 );
            VectorU ids = VectorU( _nq );
            for ( unsigned jd = 0; jd < _nq; jd++ ){
                dist1[ jd ] = fabs( _quadPointsC[ id ][0] - _quadPoints[ jd ][0] );
                if (dist1[ jd ] < 10e-12) dist1[jd] = 0;
            }
            std::iota( ids.begin(), ids.end(), 0);
            std::sort( ids.begin(), ids.end(), [&](int i,int j){return dist1[i] < dist1[j];} );
            double min1 = dist1[ ids [0] ];
            double min2 = dist1[ ids [1] ];
            _MatInterp( id, ids[0] ) = min2 / ( min1 + min2 );
            _MatInterp( id, ids[1] ) = min1 / ( min1 + min2 );
        }
    }
    //}

    // Setup Laplace Beltrami Matrices (Checked )

    _LC = SetupLaplaceBeltrami( _quadPointsC, _weightsC );
    _LFC = SetupLaplaceBeltrami( _quadPoints, _weights );

    // determine transport coefficients of Heney-Greenstein (not used)
    double g = 0.98; // Heney-Greenstein parameter
    _xi  = Matrix( 3, _nEnergies );
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        _xi( 1, n ) = 1.0 - g;
        _xi( 2, n ) = 4.0 / 3.0 - 2.0 * g + 2.0 / 3.0 * g * g;
    }
    // Stopping power
    _s = Vector( _nEnergies, 1.0 );
    _sigmaS = VectorVector( _nEnergies, Vector( _nCells, 50 ));
    _sigmaT = VectorVector( _nEnergies, Vector( _nCells, 0.01 ));

    if ( _RT ){
        //// Get Data from ICRU database
        Vector muC(_nqC);
        Vector mu(_nq);
        for( unsigned i = 0; i< _nqC; i++){ muC[i] = _quadPointsC[i][0]; }
        for( unsigned i = 0; i< _nq; i++){ mu[i] = _quadPoints[i][0]; }

        // FP Variante

        // Get TransportCoefficients -> _xi ;
        ICRU database( muC, _energies, _settings );
        _xi = Matrix(4, _nEnergies);
        database.GetTransportCoefficients( _xi );
        ////// Get Stopping Power from database
        database.GetStoppingPower( _s ); // dim _s = _nEnergies

        /*
        /// CSD variante

        Vector MUVec( _nqC * _nqC);
        Vector MUVec2( _nqC * _nq);
        for( unsigned i = 0; i < _nqC; i++){
            for( unsigned j = 0; j < _nqC; j++){
                MUVec[ i * _nqC + j ] = muC[i] * muC[j];
            }
            for( unsigned j = 0; j < _nq; j++){
                MUVec2[ i * _nq + j ] = muC[i] * mu[j];
            }
        }
        ICRU database2( MUVec, _energies, _settings);
        Matrix total;
        database2.GetAngularScatteringXS( total, _sigmaTECoarse);
        database2.GetStoppingPower(_s);
        ICRU database3( MUVec2, _energies, _settings);
        Matrix total2;
        database3.GetAngularScatteringXS( total2, _sigmaTE);

        _sigmaSECoarse = std::vector<Matrix>( _nEnergies, Matrix(_nqC, _nqC ));
        _sigmaSE = std::vector<Matrix>( _nEnergies, Matrix(_nqC, _nq ));
        for( unsigned n = 0; n < _nEnergies; ++n ) {
            for( unsigned i = 0; i < _nqC; ++i ) {
                for( unsigned j = 0; j < _nqC; ++j ) {
                    _sigmaSECoarse[n]( i, j ) = total( i * _nqC + j, n );
                }
                for( unsigned j = 0; j < _nq; ++j ) {
                    _sigmaSE[n]( i, j ) = total2( i * _nq + j, n );
                }
            }
        }
        std::cout << "Check" << std::endl;

        _scatteringKernelFC = blaze::DiagonalMatrix<Matrix>(_nq, _nq);
        for ( unsigned i = 0; i < _nq; i ++){
            _scatteringKernelFC(i,i) = _weights[i];
        }
        _scatteringKernel = blaze::DiagonalMatrix<Matrix>(_nqC, _nqC);
        for ( unsigned i = 0; i < _nqC; i ++){
            _scatteringKernel(i,i) = _weightsC[i];
        }
        std::cout << "Check2" << std::endl;

        FILE * Sigma = fopen((_settings->GetLogDir() + "/SigmaTE").c_str(), "a+");
        for( unsigned i = 0; i < _nEnergies; i++){
            fprintf(Sigma, "%f \t %f \t %u \t %u \n", _sigmaTE[i], _sigmaTECoarse[i], _sigmaSE[i].nonZeros(), _sigmaSECoarse[i].nonZeros());
        }
        */
        /*
         FILE * StopPwr = fopen((_settings->GetLogDir() + "/StopPwr").c_str(), "a+");
         for ( unsigned i = 0; i< _nEnergies; i++)
             fprintf( StopPwr, "%f \t %f \n", _energies[i], _s[i]);
         fclose(StopPwr);
         */
        std::cout << "Check" << std::endl;

    }

    /*
    ///// External Source
    _QextC = _Q; // save external source in coarse quadrature (needed for local refinement)
    _settings->SetNQuadPoints( _nq );
    _Q = _problem->GetExternalSource( _energies );  // get external source in fine quadrature
    _QextF = _Q; // save external source in fine quadrature (needed for local refinement)
    _settings->SetNQuadPoints( _nqC );
    */

    ///// Initial condition ( is zero ) (Checked)
    _solFC = VectorVector( _nCells, Vector( _nq, 0.0 ));
    _sol = VectorVector( _nCells, Vector( _nqC, 0.0));
    //_sol = _problem->SetupIC(); // This is in coarse quadrature
    _solNew = _sol;
    _solNewFC = _solFC; // I am just here to get the right size!

    ///// First Collision Source initialize
    _QFirstCollision = VectorVector( _nCells, Vector( _nqC ));

    // Density
    _density = std::vector<double>( _nCells, 1.0 );
    std::vector<double> densities( 399, 0.03 );    // water layer
    std::vector<double> densities2( 200, 1.04 );    // muscle layer
    std::vector<double> densities3( 400, 1.85 );    // water layer

    // std::vector<double> bone( 167, 1.85 );
    // std::vector<double> lung( 665, 0.03);    // maybe this is just the lung tissue and after that there should be air?

    // Concatenate layers
    densities.insert( densities.end(), densities2.begin(), densities2.end() );
    densities.insert( densities.end(), densities3.begin(), densities3.end() );
    //_density = densities;
    std::cout << " Check Build " << std::endl;
}

void FirstCollisionCSDSNSolver::Solve(){

    PrepareVolumeOutput();
    DrawPreSolverOutput();
    SolverPreprocessing(); // --> CSD transformations
    std:: cout << "Check Solver Preprocessing" << std::endl;
    _maxIter = _nEnergies - 1;
    for ( unsigned idx_energy = 0; idx_energy < _maxIter; idx_energy++){
        /*
        if ( _settings->GetIsLocalRefine() == true ){
            IterRefinement( idx_energy);
            std::cout << "Check Local Refinement; NQ = " << _nq << std::endl;
        }
        */
        // Preprocessing
        IterPreprocessing( idx_energy );
        std::cout << "Check Preproc" << std::endl;
        // Uncollided Flux
        FluxUpdateUncollided();
        std::cout << "Check flux uncoll" << std::endl;
        FVMUpdateUncollided( idx_energy );
        std::cout << "Check FVM uncoll" << std::endl;

        /*
        // Artificial Scattering
        if ( _settings->GetIsArtificialScattering() == true ){
            AddArtificialScattering( idx_energy );
            std::cout << "Check Artificial Scattering" << std::endl;
        }
        */
        ComputeFirstCollisionSource( idx_energy );
        std::cout << "Check FC Source" << std::endl;
        // Collided Flux
        FluxUpdateCollided( );
        std::cout << "Check flux coll" << std::endl;
        FVMUpdateCollided( idx_energy );
        std::cout << "Check fvm coll" << std::endl;
        // Post Processing
        IterPostprocessing( idx_energy );
        std::cout << "Check postproc" << std::endl;

        // Output Generators
        WriteVolumeOutput( idx_energy );
        std::cout << "Check Write Volume" << std::endl;
        WriteScalarOutput( idx_energy );
        std::cout << "Check Write Scalar" << std::endl;
        PrintScreenOutput( idx_energy );
        PrintHistoryOutput( idx_energy );
        PrintVolumeOutput( idx_energy );
        std::cout << "Check Write" << std::endl;
    }

    // Post Solver Processing
    DrawPostSolverOutput();

}

////// --------- Solver Functions -------- /////

void FirstCollisionCSDSNSolver::SolverPreprocessing(){

    _energiesOrig = _energies;

    // setup incoming BC on left hard coded IC
    _sol    = VectorVector( _nCells, Vector( _nqC, 0.0 ) );
    _solNew = _sol;
    for( unsigned k = 0; k < _nq; ++k ) { // -> not used for RT case
        if( _quadPoints[k][0] > 0 && !_RT ) _sol[0][k] = 1e5 * exp( -10.0 * pow( 1.0 - _quadPoints[k][0], 2 ) );
    }

    // hard coded Boundaries:
    _boundaryCells[0] = BOUNDARY_TYPE::DIRICHLET;
    _boundaryCells[_nCells - 1 ]  = BOUNDARY_TYPE::DIRICHLET;

    // Setup Identity Matrices
    _identityFC = blaze::IdentityMatrix<double>( _nq );
    _identity = blaze::IdentityMatrix<double>( _nqC );

// --- Transformation for CSD solution ---
    // do substitution from psi to psiTildeHat (cf. Dissertation Kerstion Kuepper, Eq. 1.23)
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nqC; ++k ) {
            _sol[j][k] = _sol[j][k] * _density[j] * _s[_nEnergies - 1];    // note that _s[_nEnergies - 1] is stopping power at highest energy
        }
        for( unsigned k = 0; k < _nq; ++k ) {
            _solFC[j][k] = _solFC[j][k] * _density[j] * _s[_nEnergies - 1];    // note that _s[_nEnergies - 1] is stopping power at highest energy
        }
    }

    // store transformed energies ETilde instead of E in _energies vector (cf. Dissertation Kerstion Kuepper, Eq. 1.22)
    double tmp   = 0.0;
    _energies[0] = 0.0;
    for( unsigned n = 1; n < _nEnergies; ++n ) {
        tmp          = tmp + _dE * 0.5 * ( 1.0 / _s[n] + 1.0 / _s[n - 1] );
        _energies[n] = tmp;
    }
    // store transformed energies ETildeTilde instead of ETilde in _energies vector (cf. Dissertation Kerstion Kuepper, Eq. 1.25)
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        _energies[n] = _energies[_nEnergies - 1] - _energies[n];
    }

    // determine minimal density for CFL computation
    double densityMin = _density[0];
    for( unsigned j = 1; j < _nCells; ++j ) {
        if( densityMin > _density[j] ) densityMin = _density[j];
    }

    // cross sections do not need to be transformed to ETilde energy grid since e.g. TildeSigmaT(ETilde) = SigmaT(E(ETilde))
}

void FirstCollisionCSDSNSolver::IterPreprocessing( unsigned idx_energy ){

    _dE = fabs( _energies[ idx_energy + 1]  - _energies[ idx_energy ]);

    double xi1 = _xi( 1, _nEnergies - idx_energy - 1 );
    double xi2 = _xi( 2, _nEnergies - idx_energy - 1 );
    double xi3 = _xi( 3, _nEnergies - idx_energy - 1 );

    // setup coefficients in FP step
    if( _FPMethod == 1 ) {
        alpha  = 0.0; // was 0
        alpha2 = xi1 / 2.0;
        beta   = 0.0;
    }
    else if( _FPMethod == 2 ) {
        alpha  = xi1 / 2.0 + xi2 / 8.0;
        alpha2 = 0.0;
        beta   = xi2 / 8.0 / xi1;
    }
    else if( _FPMethod == 3 ) {
        alpha  = xi2 * ( 27.0 * xi2 * xi2 + 5.0 * xi3 * xi3 - 24.0 * xi2 * xi3 ) / ( 8.0 * xi3 * ( 3.0 * xi2 - 2.0 * xi3 ) );
        beta   = xi3 / ( 6.0 * ( 3.0 * xi2 - 2.0 * xi3 ) );
        alpha2 = xi1 / 2.0 - 9.0 / 8.0 * xi2 * xi2 / xi3 + 3.0 / 8.0 * xi2;
    }

    // write BC for water phantom ( from Olbrant ) for uncollided Calculation
    if( _RT ) {
        for( unsigned k = 0; k < _nq; ++k ) {
            if( _quadPoints[k][0] > 0 ) {
                _solFC[0][k] = 1e5 * exp( -200.0 * pow( 1.0 - _quadPoints[k][0], 2 ) ) *
                               exp( -50.0 * pow( _energyMax - _energiesOrig[_nEnergies - idx_energy - 1], 2 ) ) * _density[0] *
                               _s[_nEnergies - idx_energy - 1];
            }
        }
    }

    // add FP scattering term implicitly to collided part
    IL = _identity - beta * _LC;
    ILFC = _identityFC - beta * _LFC;

    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // _sol[j] = blaze::solve( _identity - _dE * alpha2 * _LC, _sol[j] );
        _sol[j] = IL * blaze::solve( IL - _dE * _sigmaS[idx_energy][j] * alpha * _LC, _sol[j] );
    }


}

void FirstCollisionCSDSNSolver::FluxUpdateUncollided(){

    _psiDx = VectorVector( _nCells, Vector( _nq, 0.0 ) );
    _psiDy = VectorVector( _nCells, Vector( _nq, 0.0 ) );

    if( _reconsOrder > 1 ) {
        _mesh->ReconstructSlopesU( _nq, _psiDx, _psiDy, _solFC );    // unstructured reconstruction
    }

#pragma omp parallel for
    for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
        // Dirichlet Boundaries stay at IC
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        //_solNewFC[idx_cell] = ConstructFluxSN( idx_cell, true );
        for( unsigned idx_quad = 0; idx_quad < _nq ; ++idx_quad ) {
            _solNewFC[idx_cell][idx_quad] = 0.0;
            // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
            for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); ++idx_neighbor ) {
                // store flux contribution on psiNew_sigmaS to save memory
                if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                    continue; // adiabatic
                    // flux[idx_quad] +=
                    //    _g->Flux( quadPoints[idx_quad], sol[idx_cell][idx_quad], sol[idx_cell][idx_quad], _normals[idx_cell][idx_neighbor] );
                else {
                    switch( _reconsOrder ) {
                        // first order solver
                        case 1:
                            _solNewFC[idx_cell][idx_quad] += _g->Flux( _quadPoints[idx_quad],
                                                                       _solFC[idx_cell][idx_quad] / _density[idx_cell],
                                                                       _solFC[_neighbors[idx_cell][idx_neighbor]][idx_quad] / _density[_neighbors[idx_cell][idx_neighbor]],
                                                                       _normals[idx_cell][idx_neighbor] );
                            break;
                            // second order solver
                        default:
                            ErrorMessages::Error("Choose ReconsOrder = 1!", "FluxUpdate FCS");
                    }
                }
            }
        }
    }
}


void FirstCollisionCSDSNSolver::FVMUpdateUncollided( unsigned idx_energy ){

#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;

        // loop over all ordinates
        for( unsigned idx_quad = 0; idx_quad < _nq; ++idx_quad ) {
            // time update angular flux with numerical flux and total scattering cross section
            _solNewFC[idx_cell][idx_quad] = _solFC[idx_cell][idx_quad] - ( _dE / _areas[idx_cell] ) * _solNewFC[idx_cell][idx_quad] ;
            _solNewFC[idx_cell][idx_quad] -= _sigmaT[idx_energy][idx_cell] * _solFC[idx_cell][idx_quad];
            _solNewFC[idx_cell][idx_quad] -= _sigmaS[idx_energy][idx_cell] * _solFC[idx_cell][idx_quad];
            // only for CSD without FP
            // _solNewFC[idx_cell][idx_quad] -=  _dE * _sigmaTE[_nEnergies - idx_energy - 1] * _solFC[idx_cell][idx_quad];
            //_solNewFC[idx_cell][idx_quad] -= _dE * _xi( 0, _nEnergies - idx_energy - 1 ) * _solFC[idx_cell][idx_quad];
            // + external source ( not included )
        }
        if ( _settings->GetIsArtificialScattering()){
            // Add Artificiel Scattering Terms to actual solution
            _solFC[ idx_cell ] += M_PI * _xiAS[1] * _sigmaAS * _LFC * _solFC[ idx_cell ];
        }

        /*
        Vector sol1 ( _nq );
        sol1 = blaze::solve( ILFC, _solNewFC[idx_cell] ); // @todo ODER _solFC[idx_cell] -> noch zu testen
        _solNewFC[idx_cell] += _dE * alpha * _LFC * sol1;
        */
    }
}

void FirstCollisionCSDSNSolver::ComputeFirstCollisionSource( unsigned idx_energy ){

    /*
   // Version mit Cross Sections CSD
    for ( unsigned id = 0; id < _nCells; id ++){
        if( _boundaryCells[id] == BOUNDARY_TYPE::DIRICHLET ) continue;
        _QFirstCollision[id] = _sigmaSE[ _nEnergies - idx_energy -1 ] *  _scatteringKernelFC  *  _solNewFC[id];
    }
    //_QFirstCollision = VectorVector(_nCells, Vector (_nqC, 0.0));
    */


    // FP CSD version
    for (unsigned id = 0; id < _nCells; id ++){
        Vector sol =  _solNewFC[id];
        Vector sol1 = blaze::solve( ILFC, sol );
        _QFirstCollision[id] =_MatInterp *  (  _sigmaS[idx_energy][id] * alpha * _LFC * sol1  + _sigmaS[idx_energy][id] * _solFC[id] );
    }

    /*
     for( unsigned i = 0; i <_nCells; i++){
         Vector sol1 = _solNewFC[i];
         _QFirstCollision[i] = _MatInterp *  ( ILFC * blaze::solve( ILFC - _dE * alpha * _LFC, sol1 ) - sol1 );
     }
     */
    /*
    Matrix ILFC = _identityFC - beta * _LFC;
     // real version
     for (unsigned id = 0; id < _nCells; id ++){
         Vector sol1 = blaze::solve( ILFC,  _solNewFC[id] );
         _QFirstCollision[id] =  _MatInterp * alpha * _LFC * sol1;
     }
     */
    /*
    // Version mit \xi 0
    Matrix IL = _identityFC - beta * _LFC;
    for ( unsigned id = 0; id < _nCells; id ++){
        if( _boundaryCells[id] == BOUNDARY_TYPE::DIRICHLET ) continue;
        Vector solTest = _solNewFC[id] + _dE * _xi(0, _nEnergies - idx_energy - 1 ) / 4 / M_PI * _solFC[id]; // add again total XS
        //Vector sol = blaze::solve(_identityFC - alpha2 * _LFC, solTest );
        Vector sol = blaze::solve( IL, solTest );
        _QFirstCollision[id] = _MatInterp * alpha * _LFC * sol;
    }
    */


    /*
    FILE * FCSource = fopen((_settings->GetLogDir() + "/FCSource.txt").c_str(), "a+");
    FILE * SOL = fopen((_settings->GetLogDir() + "/PsiNew.txt").c_str(), "a+");
    // std::cout << "FC Kernel size: (" << _scatteringKernelFC.rows() << "," << _scatteringKernelFC.columns() << ")" << std::endl;
    // std::cout << "solNewFC size: (" << _solNewFC[0].size() << ")" << std::endl;
    // std::cout << "QFirstCollision size: (" << _QFirstCollision.size() << " , " << _QFirstCollision[0].size() << ")" << std::endl;
    if ( idx_energy == 100){

    for (unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
    //    _QFirstCollision[ idx_cell ] = _sigmaSE[ _nEnergies - idx_energy -1 ] * _scatteringKernelFC * _solNewFC[idx_cell];
        for( unsigned idx_quad = 0; idx_quad < _nqC; idx_quad ++){
            //if ( _QFirstCollision[idx_cell][idx_quad] 0)
                fprintf(FCSource, "FCSource(%u,%u) = \t %f \n", idx_cell, idx_quad, _QFirstCollision[idx_cell][idx_quad] );
            //if ( _solNewFC[idx_cell][idx_quad] > 10e-12)
                fprintf(SOL, "Psi(%u,%u) = \t %f \n", idx_cell, idx_quad, _solNewFC[idx_cell][idx_quad] );
        }
    }

    }
    */
}


void FirstCollisionCSDSNSolver::FluxUpdateCollided( ){

    _psiDx = VectorVector( _nCells, Vector( _nqC, 0.0 ) );
    _psiDy = VectorVector( _nCells, Vector( _nqC, 0.0 ) );
    if( _reconsOrder > 1 ) {
        _mesh->ReconstructSlopesU( _nqC, _psiDx, _psiDy, _sol );    // unstructured reconstruction
    }
#pragma omp parallel for
    for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
        // Dirichlet Boundaries stayd
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // _solNew[idx_cell] = ConstructFluxSN( idx_cell, false );
        for( unsigned idx_quad = 0; idx_quad < _nqC ; ++idx_quad ) {
            _solNew[idx_cell][idx_quad] = 0.0;
            // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
            for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); ++idx_neighbor ) {
                // store flux contribution on psiNew_sigmaS to save memory
                if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                    continue; // adiabatic
                    // flux[idx_quad] +=
                    //    _g->Flux( quadPoints[idx_quad], sol[idx_cell][idx_quad], sol[idx_cell][idx_quad], _normals[idx_cell][idx_neighbor] );
                else {
                    switch( _reconsOrder ) {
                        // first order solver
                        case 1:
                            _solNew[idx_cell][idx_quad] += _g->Flux( _quadPointsC[idx_quad],
                                                                     _sol[idx_cell][idx_quad] / _density[idx_cell],
                                                                     _sol[_neighbors[idx_cell][idx_neighbor]][idx_quad] / _density[_neighbors[idx_cell][idx_neighbor]],
                                                                     _normals[idx_cell][idx_neighbor] );
                            break;
                            // second order solver
                        default:
                            ErrorMessages::Error("Choose ReconsOrder = 1!", "FluxUpdate FCS");
                    }
                }
            }
        }
    }
}

void FirstCollisionCSDSNSolver::FVMUpdateCollided( unsigned idx_energy ){

#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // loop over all ordinates
        for( unsigned idx_quad = 0; idx_quad < _nqC; ++idx_quad ) {
            // time update angular flux with numerical flux
            _solNew[ idx_cell ][ idx_quad ] = _sol[idx_cell][idx_quad] - ( _dE / _areas[idx_cell] ) * _solNew[idx_cell][idx_quad] ;
            _solNew[ idx_cell ][ idx_quad ] -= _sigmaT[idx_energy][idx_cell] * _sol[idx_cell][idx_quad];
            // _solNew[ idx_cell ][ idx_quad ] -= _dE * _sigmaTECoarse[_nEnergies -idx_energy - 1] * _sol[ idx_cell ][ idx_quad ];
        }
        // scattering -> was added before
        // _solNew[ idx_cell ] += _dE * ( _sigmaSECoarse[_nEnergies - idx_energy - 1] * _scatteringKernel * _sol[ idx_cell ] );
        // First Collision Source
        _solNew[ idx_cell ] += _dE * _QFirstCollision[ idx_cell ];
    }
}
void FirstCollisionCSDSNSolver::IterPostprocessing( ){}
void FirstCollisionCSDSNSolver::IterPostprocessing( unsigned iter ){
    _sol = _solNew;
    _solFC = _solNewFC;
    ComputeRadFlux( iter ); // -> Get new _fluxNew
}




///// -------- Flux Constructor ------ /////

Vector FirstCollisionCSDSNSolver::ConstructFluxSN( unsigned idx_cell, bool FC ){

    unsigned nq;
    VectorVector sol;
    VectorVector quadPoints;
    if ( FC ){ // i.e. FC Calculation
        nq = _nq;
        quadPoints = _quadPoints;
        sol = _solFC;
    }
    else { // i.e. Coarse Calculation
        nq = _nqC;
        quadPoints = _quadPointsC;
        sol = _sol;
    }
    Vector flux( nq, 0.0 );
    double psiL;
    double psiR;

    for( unsigned idx_quad = 0; idx_quad < nq ; ++idx_quad ) {
        // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
        for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); ++idx_neighbor ) {
            // store flux contribution on psiNew_sigmaS to save memory
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                continue; // adiabatic
                // flux[idx_quad] +=
                //    _g->Flux( quadPoints[idx_quad], sol[idx_cell][idx_quad], sol[idx_cell][idx_quad], _normals[idx_cell][idx_neighbor] );
            else {
                switch( _reconsOrder ) {
                    // first order solver
                    case 1:
                        flux[idx_quad] += _g->Flux( quadPoints[idx_quad],
                                                    sol[idx_cell][idx_quad] / _density[idx_cell],
                                                    sol[_neighbors[idx_cell][idx_neighbor]][idx_quad] / _density[_neighbors[idx_cell][idx_neighbor]],
                                                    _normals[idx_cell][idx_neighbor] );
                        break;
                        // second order solver
                    case 2:
                        // left status of interface
                        psiL =  sol[idx_cell][idx_quad] +
                                _psiDx[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[idx_cell][0] ) +
                                _psiDy[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[idx_cell][1] );
                        // right status of interface
                        psiR =  sol[_neighbors[idx_cell][idx_neighbor]][idx_quad] +
                                _psiDx[_neighbors[idx_cell][idx_neighbor]][idx_quad] *
                                ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[_neighbors[idx_cell][idx_neighbor]][0] ) +
                                _psiDy[_neighbors[idx_cell][idx_neighbor]][idx_quad] *
                                ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[_neighbors[idx_cell][idx_neighbor]][1] );
                        // positivity check (if not satisfied, deduce to first order)
                        if( psiL < 0.0 || psiR < 0.0 ) {
                            psiL = sol[idx_cell][idx_quad];
                            psiR = sol[_neighbors[idx_cell][idx_neighbor]][idx_quad];
                        }
                        // flux evaluation
                        flux[idx_quad] += _g->Flux( quadPoints[idx_quad],
                                                    psiL / _density[idx_cell],
                                                    psiR / _density[_neighbors[idx_cell][idx_neighbor]],
                                                    _normals[idx_cell][idx_neighbor] );
                        break;
                        // higher order solver
                    case 3:
                        std::cout << "higher order is WIP" << std::endl;
                        break;
                        // default: first order solver
                    default:
                        flux[idx_quad] += _g->Flux( quadPoints[idx_quad],
                                                    sol[idx_cell][idx_quad] / _density[ idx_cell ],
                                                    sol[_neighbors[idx_cell][idx_neighbor]][idx_quad] / _density[_neighbors[idx_cell][idx_neighbor]],
                                                    _normals[idx_cell][idx_neighbor] );
                }
            }
        }
    }
    return flux;
}

//// ------------ Helper Functions for SN ---------- ////
void FirstCollisionCSDSNSolver::ComputeRadFlux( ){}
void FirstCollisionCSDSNSolver::ComputeRadFlux( unsigned iter ){

    _flux = _fluxNew;
    //FILE * radFluxfile = fopen((_settings->GetLogDir() + "/radFluxFC").c_str(), "a+");
    //fprintf(radFluxfile, "\n Iteration: \n ");
    //fprintf(radFluxfile, "\t uncoll \t coll \t \t cumulative \n \n ");
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        // flux from uncollided
        double fluxUC = blaze::dot( _solFC[idx_cell], _weights );
        // flux from collided
        double fluxC = blaze::dot( _sol[idx_cell], _weightsC );
        // cumulative flux
        _fluxNew[idx_cell] = fluxUC + fluxC;

        //_fluxNew[idx_cell] = blaze::dot( _sol[idx_cell], _weights );
        //fprintf(radFluxfile, "\t %f \n", (float) _fluxNew[idx_cell] * 1000);
        /*
        if( iter == _maxIter- 1 || iter == 700 || iter == 600 ){
            std::cout << idx_cell << ", " << _fluxNew[idx_cell] << std::endl;
            fprintf(radFluxfile, "%f\t%f\t%f\n", (float) fluxUC , (float) fluxC , (float) _fluxNew[idx_cell] );
        }
        */
    }
    if(( _settings->GetVolumeOutputFrequency() != 0 && iter % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
       ( iter == _maxIter - 1 ) ){
        FILE * DoseFile = fopen( (_settings->GetLogDir() + "/DoseFile" + std::to_string(iter)).c_str(), "a+");
        FILE * FluxFile = fopen( (_settings->GetLogDir() + "/FluxFile" + std::to_string(iter)).c_str(), "a+");
        double dose;

        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++){
            if ( idx_cell > 0){
                dose = 0.5 * _dE *( _fluxNew[idx_cell] * _s[_nEnergies - iter - 1] + _flux[idx_cell] * _s[_nEnergies - iter ] ) /
                       _density[idx_cell];
            }
            else{
                dose = _dE * _fluxNew[idx_cell] * _s[_nEnergies - iter - 1] / _density[idx_cell];
            }
            fprintf(DoseFile, "%f\n", (float) dose);
            fprintf(FluxFile, "%f\n", (float) _fluxNew[idx_cell]);
        }
        fclose(DoseFile);
        fclose(FluxFile);
    }
}

void FirstCollisionCSDSNSolver::GetASTransportCoefficients() {

    // Calculate Artificial Scattering Kernel
    Matrix muM( _nq, _nq);
    Vector muV( _nq * _nq );
    double sum;
    Matrix ASKernel(_nq, _nq, 0.0);
    Vector ASV( _nq *_nq );
    for ( unsigned id = 0; id < _nq; id ++){
        for ( unsigned jd = 0; jd < _nq; jd ++){
            double mu =_quadPoints[id][0] * _quadPoints[jd][0] ;
            muM( id, jd ) = mu;
            double expTerm = ( - ( ( mu - 1 ) * ( mu - 1 ) ) * _nq * _nq ) / ( _betaAS * _betaAS );
            ASKernel( id, jd ) = _nq / _betaAS * exp( expTerm );
            sum += ASKernel( id, jd );
        }
        for ( unsigned jd = 0; jd < _nq; jd++){
            ASKernel( id , jd ) /= sum;        // Normalize per Row
        }
    }
    // Get Matrices in Vector Format
    for( unsigned i = 0; i < muM.rows(); ++i ) {
        for( unsigned j = 0; j < muM.columns(); ++j ) {
            muV[i * muM.columns() + j] = std::fabs( muM( i, j ) );
            ASV[i * ASKernel.columns() + j] = ASKernel(i,j);
        }
    }

    _xiAS = Vector(3);
    for( unsigned n = 0; n < _xiAS.size(); n++){
        _xiAS[n] = 0.0;
        for ( unsigned k = 0; k <_nq; k++)
            _xiAS[n] -= 0.5 * ( pow( 1.0 - muV[k + 1], n ) * ASV[k + 1] + pow( 1.0 - muV[k], n ) * ASV[k] ) * ( muV[k + 1] - muV[k] );
    }
}

void FirstCollisionCSDSNSolver::Refinement(){
    std::cout << "Refinemnt Entered" << std::endl;

    _nq = _nqC * _nqF;
    _quadPoints = VectorVector( _nq, Vector( 3 ) );
    _weights = Vector( _nq );

    unsigned minID = _nqC;
    unsigned maxID = 0;
    double weightsum = 0.0;

    unsigned ind = 0;
    for( unsigned id = 0; id < _nqC; id ++ ){
        if (  _quadPointsC[id][0] < 0.7 ){ // these are the points that are not refined
            _quadPoints[ ind ] = _quadPointsC[ id ];
            _weights[ ind ] = _weightsC[ id ];
            ind += 1;
        }
        else{
            if( id < minID ){ minID = id; }
            if( id > maxID ){ maxID = id; }
            weightsum += _weightsC[id];
        }
    }
    std::cout << "Check Determine" << std::endl;
    std::cout << "MinID: " << minID << "MaxId: " << maxID << std::endl;

    double ws = 0.0;
    for( unsigned id = 0; id < _nqF; id ++ ){
        ws += _weightsF[id]; // sum of actual weights
    }

    double qMax, qMin;
    if (maxID == _nqC - 1 ){ qMax = 1.0; }
    else{ qMax = _quadPointsC[maxID + 1][0]; }
    if (minID == 0 ){ qMin = -1.0; }
    else{ qMin = _quadPointsC[minID - 1][0]; }
    std::cout << "QMin: " << qMin << "QMax: " << qMax << std::endl;

    for ( unsigned id = 0; id < _nqF; id++){
        _quadPoints[ind] = _quadPointsF[id] * (qMax - qMin) / 2 + ( qMax + qMin ) / 2;
        _weights[ind] = _weightsF[id] / ws * weightsum; // adapt weights
        ind += 1;
    }
    _nq = ind;
    _quadPoints.resize( _nq );
    _weights.resize( _nq );
    /*
    for ( unsigned id = 0; id < _nqC; id++){
        if (  _quadPointsC[id][0] > 0 ){
             double mind = 20.0;
             double mind2 = 20.0;
             unsigned ID, ID2;

            for ( unsigned j = 0; j < _nq - saveind; j++){
                if( fabs( _quadPointsC[ id ][0] - _quadPoints[ j ][0]) < mind){
                    mind = fabs( _quadPointsC[ id ][0] - _quadPoints[ j ][0]);
                    ID = j;
                }
                else if( fabs( _quadPointsC[ id ][0] - _quadPoints[ j ][0]) < mind2){
                    mind2 = fabs( _quadPointsC[ id ][0] - _quadPoints[ j ][0]);
                    ID2 = j;
                }
            }
            _MatInterp( id, ID ) = mind2 / ( mind + mind2 );
            _MatInterp( id, ID2 ) = mind / ( mind + mind2 );
        }
    }
    _MatInterp.resize( _nqC, _nq);
    */
}

Matrix FirstCollisionCSDSNSolver::SetupLaplaceBeltrami( VectorVector p, Vector w ){
    unsigned nq = p.size();
    Vector mu( nq );
    for( unsigned i = 0; i < nq ; i++){ mu[i] = p[i][0]; }
    Matrix LC = Matrix( nq, nq, 0.0);
    double DMinus = 0.0;
    double DPlus = 0.0;
    for ( unsigned id = 0; id < nq; id++){
        DMinus = DPlus;
        DPlus = DMinus - 2 * mu[ id ] * w[ id ];
        if ( id > 0){
            LC( id , id - 1 ) = DMinus / ( mu[ id ] - mu[ id -1 ] ) / w[ id ];
            LC( id, id )      = - DMinus / ( mu[ id ] - mu[ id - 1] ) / w[ id ];
        }
        if( id < _nqC - 1 ){
            LC( id , id + 1 ) = DPlus / ( mu[ id + 1 ] - mu[ id ] ) / w[ id ];
            LC( id, id ) += - DPlus / ( mu[ id + 1 ] - mu[ id ] ) / w[ id ];
        }
    }
    return LC;
}

void FirstCollisionCSDSNSolver::FluxUpdate(){}             // not used
void FirstCollisionCSDSNSolver::FVMUpdate( unsigned ){}    // not used


// -------- Output generators ( copied from SN_Solver ) -------------

void FirstCollisionCSDSNSolver::PrepareVolumeOutput() {
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();

    _outputFieldNames.resize( nGroups );
    _outputFields.resize( nGroups );

    // Prepare all OutputGroups ==> Specified in option VOLUME_OUTPUT
    for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
        // Prepare all Output Fields per group

        // Different procedure, depending on the Group...
        switch( _settings->GetVolumeOutput()[idx_group] ) {
            case MINIMAL:
                // Currently only one entry ==> rad flux
                _outputFields[idx_group].resize( 1 );
                _outputFieldNames[idx_group].resize( 1 );

                _outputFields[idx_group][0].resize( _nCells );
                _outputFieldNames[idx_group][0] = "radiation flux density";
                break;

            case MEDICAL:
                _outputFields[idx_group].resize( 2 );
                _outputFieldNames[idx_group].resize( 2 );

                // Dose
                _outputFields[idx_group][0].resize( _nCells );
                _outputFieldNames[idx_group][0] = "dose";
                // Normalized Dose
                _outputFields[idx_group][1].resize( _nCells );
                _outputFieldNames[idx_group][1] = "normalized dose";
                break;

            default: ErrorMessages::Error( "Volume Output Group not defined for First Collision SN Solver!", CURRENT_FUNCTION ); break;
        }
    }
}


void FirstCollisionCSDSNSolver::WriteVolumeOutput( unsigned idx_pseudoTime ) {
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();
    double maxDose;
    //FILE * OutputFlux = fopen((_settings->GetLogDir() + "/OutputFlux").c_str(), "a+" );
    //FILE * OutputDose = fopen((_settings->GetLogDir() + "/OutputDose").c_str(), "a+" );
    //FILE * OutputDoseN = fopen((_settings->GetLogDir() + "/OutputDoseN").c_str(), "a+" );
    if( ( _settings->GetVolumeOutputFrequency() != 0 && idx_pseudoTime % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
        ( idx_pseudoTime == _maxIter - 1 ) || ( idx_pseudoTime == 1640) ) { // need sol at last iteration
        for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
            switch( _settings->GetVolumeOutput()[idx_group] ) {
                case MINIMAL:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = _fluxNew[idx_cell];
                        // fprintf(OutputFlux, "%f\n", _fluxNew[idx_cell] );
                    }
                    break;

                case MEDICAL:
                    // Compute Dose
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        if( idx_cell > 0 ) {
                            _outputFields[idx_group][0][idx_cell] +=
                                    0.5 * _dE *
                                    ( _fluxNew[idx_cell] * _s[_nEnergies - idx_pseudoTime - 1] + _flux[idx_cell] * _s[_nEnergies - idx_pseudoTime] ) /
                                    _density[idx_cell];    // update dose with trapezoidal rule
                        }
                        else {
                            _outputFields[idx_group][0][idx_cell] +=
                                    _dE * _fluxNew[idx_cell] * _s[_nEnergies - idx_pseudoTime - 1] / _density[idx_cell];
                        }
                        // fprintf(OutputDose, "%f\n", _outputFields[idx_group][0][idx_cell] );
                    }
                    // Compute normalized dose
                    _outputFields[idx_group][1] = _outputFields[idx_group][0];

                    maxDose = *std::max_element( _outputFields[idx_group][0].begin(), _outputFields[idx_group][0].end() );

                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][1][idx_cell] /= maxDose;
                        // fprintf(OutputDoseN, "%f\n", _outputFields[idx_group][1][idx_cell] );
                    }
                    break;
                default: ErrorMessages::Error( "Volume Output Group not defined for First Collision SN Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
}



///////////// following is from constructor without FP

// if (_RT) Get scattering XS and Total XS from database for CSD
/* not used bc CSD Solver not good
    Matrix muMatrixC( _nqC, _nqC );
    Matrix muMatrixF( _nqC, _nqF );
    for ( unsigned id = 0; id < muMatrixC.rows(); id ++){
        for ( unsigned jd = 0; jd < muMatrixC.columns(); jd ++){
            double inner = 0;
            for ( unsigned k = 0; k < 3; k++){ inner += _quadPointsC[id][k] * _quadPointsC[jd][k]; }
            muMatrixC( id, jd ) = abs(inner); // is abs right here??
        }
        for ( unsigned jd = 0; jd < muMatrixF.columns(); jd ++){
            double inner = 0;
            for ( unsigned k = 0; k < 3; k++){ inner += _quadPointsC[id][k] * _quadPointsF[jd][k]; }
            muMatrixF( id, jd ) = abs(inner);
        }
    }
    Vector angleVecC( muMatrixC.rows() * muMatrixC.columns() );
    Vector angleVecF( muMatrixF.rows() * muMatrixF.columns() );
    for ( unsigned id = 0; id < muMatrixC.rows(); id ++){
        for ( unsigned jd = 0; jd < muMatrixC.columns(); jd ++){ angleVecC[ id * muMatrixC.columns() + jd ] = muMatrixC( id, jd ); }
        for ( unsigned jd = 0; jd < muMatrixF.columns(); jd ++){ angleVecF[ id * muMatrixF.columns() + jd ] = muMatrixF( id, jd ); }
    }
    ICRU database( angleVecC, _energies, _settings );
    ICRU database2( angleVecF, _energies, _settings );
    Matrix total;
    Matrix total2;
    database.GetAngularScatteringXS( total, _sigmaTECoarse );
    database2.GetAngularScatteringXS( total2, _sigmaTEFine );
    _sigmaSECoarse = std::vector<Matrix>( _energies.size(), Matrix( muMatrixC.rows(), muMatrixC.columns(), 0.0 ) );
    _sigmaSEFine   = std::vector<Matrix>( _energies.size(), Matrix( muMatrixF.rows(), muMatrixF.columns(), 0.0 ) );
    for( unsigned n = 0; n < _energies.size(); ++n ) {
        for( unsigned id = 0; id < muMatrixC.rows(); ++id ) {
            for( unsigned jd = 0; jd < muMatrixC.columns(); ++jd ) { _sigmaSECoarse[n]( id, jd ) = total( id * muMatrixC.columns() + jd, n ); }
            for( unsigned jd = 0; jd < muMatrixF.columns(); ++jd ) { _sigmaSEFine[n]( id, jd ) = total( id * muMatrixF.columns() + jd, n ); }
        }
    }
*/

/* not used (CSD unstable)
    for( unsigned n = 0; n <_nEnergies; n ++){
        double sigmaS = 0.0;
        for ( unsigned id = 0; id < _nCells; id ++){
            sigmaS += _sigmaS[n][id] / ( 4 * M_PI );
        }
        _sigmaSE[n] = sigmaS; _sigmaSECoarse[n] = sigmaS; _sigmaSEFine[n] = sigmaS;
    }
*/

/*
    ///// Recompute Scattering Kernel (Multiplication with weights)
    _scatteringKernel.resize( _nqC ); // Coarse
    for ( unsigned id = 0; id < _nqC; id ++){ _scatteringKernel( id, id ) = _weightsC[ id ]; }

    _scatteringKernelFC.resize( _nq ); // fine -> recompute for local refinement
    for ( unsigned id = 0; id < _nq; id ++){ _scatteringKernel( id, id ) = _weights[ id ]; }

*/


/*
void FirstCollisionCSDSNSolver::IterRefinement( unsigned idx_energy ){

    // --- Step 0 --- Save Values from last Iteration
    _nqOld = _nq;
    _quadIDsOld = _quadIDs;
    VectorU _refineOld = _refine;

    // --- Step 1 --- Determine if and where Refinement is necessary (depending on external source)

    DetermineRefinement( idx_energy ); // ---> Get Vector _refine
    std::cout << "Check: Refinement determined" << std::endl;

    // Determine if we have to change new quadrature (depending on last Iteration)

    unsigned change = 0;
    for ( unsigned id = 0 ; id < _nqC; id++){
        if ( _refineOld[ id ] != 0 && _refine[ id ] == 0 ) { change = 1; }       // was refined in last iter, has not to be refined now
        else if ( _refineOld[ id ] == 0 && _refine[ id ] != 0 ) { change = 1; }  // was not refined in last iter, has to be refined now
    }

    // --- Step 2 --- Refinement

    // --- NO need to change quadrature in this Iteration-> only adapt external Source, if not same for all energies
    if ( change == 0 ){
        // Adapt external Source
        if (_Q.size() != 1u ){
            for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
                // Isotropic
                if( _Q[idx_energy][idx_cell].size() == 1u ) {
                    if ( _refine == 0) { _Q[idx_energy] = _QextC[idx_energy]; }
                    else { _Q[idx_energy] = _QextF[idx_energy]; } }
                else{
                    _Q[idx_energy][idx_cell] = Vector( _nq );
                    for ( unsigned  id = 0; id < _nq; id++){
                        // if _quadPoint[id] is a coarse point
                        if ( _quadIDs[ id ] < _nqC ) { _Q[idx_energy][idx_cell][ id ] = _QextC[idx_energy][idx_cell][ _quadIDs[id] ]; }
                        else    { _Q[idx_energy][idx_cell][ id ] = _QextF[idx_energy][idx_cell][ _quadIDs[id] - _nqC ]; }
                    }
                }
            }
        }
        std::cout << "Check: Quadarture not changed" << std::endl;
        return;
    }

    // --- Step 2.1 --- Refining the quadrature
    RefineQuadrature(); // ----> Get _nq, _quadPoints, _weights, _quadIDs
    std::cout << "Check: Quadrature adapted" << std::endl;

    // --- Step 2.2 --- Interpolation of next _solFC Variable
    InterpolateSolution(); // ----> Get _solFC
    std::cout << "Check: Solution interpolated" << std::endl;

    // --- Step 2.3 --- Computation of new scattering Kernel ( scattering XS )
    _scatteringKernelFC = 0.0;
    for( unsigned id = 0; id < _nq; id ++){
        _scatteringKernelFC( id, id ) = _weights(id);
        if ( _quadIDs[ id ] < _nqC)
            for ( unsigned jd = 0 ; jd < _nqC; jd ++){ _sigmaSE[idx_energy]( jd, id ) = _sigmaSECoarse[idx_energy]( jd, _quadIDs[ id ] ); }
        else
            for ( unsigned jd = 0 ; jd < _nqC; jd ++){ _sigmaSE[idx_energy]( jd, id ) = _sigmaSEFine[idx_energy]( jd, _quadIDs[ id ] - _nqC ); }
    }

    // _scatteringKernelFC = KERNEL->GetScatteringKernelFirstCollision( _nq, _weights );
    // std::cout << "Check: Scattering Kernel calculated" << std::endl;

    // --- Step 3 --- Set external Source dependent on quadrature points

    for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
        if( _Q.size() == 1u ) {                   // constant source for all energies
            if( _Q[0][idx_cell].size() == 1u ) { }  // isotropic source
            else{
                _Q[0][idx_cell] = Vector( _nq );
                for ( unsigned  id = 0; id < _nq; id++){
                    if ( _quadIDs[ id ] < _nqC )  // i.e. point of coarse quadrature
                            { _Q[0][idx_cell][ id ] = _QextC[0][idx_cell][ _quadIDs[id] ]; }
                    else    { _Q[0][idx_cell][ id ] = _QextF[0][idx_cell][ _quadIDs[id] - _nqC ]; }
                }
            }
        }
        else {
            if( _Q[idx_energy][idx_cell].size() == 1u ) {}   // isotropic source
            else{
                _Q[idx_energy][idx_cell] = Vector( _nq );
                for ( unsigned  id = 0; id < _nq; id++){
                    if ( _quadIDs[ id ] < _nqC ) // i.e. point of coarse quadrature
                            { _Q[idx_energy][idx_cell][ id ] = _QextC[idx_energy][idx_cell][ _quadIDs[id] ]; }
                    else    { _Q[idx_energy][idx_cell][ id ] = _QextF[idx_energy][idx_cell][ _quadIDs[id] - _nqC ];
                    }
                }
            }
        }
    }
    std::cout << "Check: External source adapted" << std::endl;

}
*/

//// ------------ Helper Functions for Refinement ---------- ////

/*
void FirstCollisionCSDSNSolver::DetermineRefinement( unsigned idx_energy ){

    VectorU refine( _nqC, 0); // --> refine in {0,1} -> need refinement in this iteration or not

    // loop over COARSE external Source ( if _Qext[idx_energy][..][idx_quad] != 0 have to refine in direction idx_quad )
    for (unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++){
        if( _QextC.size() == 1u ) {                     // constant source for all energies
            if( _QextC[0][idx_cell].size() == 1u )          // isotropic source
                if( _QextC[0][idx_cell][0] == 0.0) continue;    // Isotropic source is 0 -> nothing to do
                else { refine = VectorU( _nqC, 1 ); break; }      // Isotropic Source is not 0 -> Refinement
            else{                                          // anisotropic source
                for (unsigned idx_quad = 0; idx_quad < _nqC; idx_quad++){
                    if( _QextC[0][idx_cell][idx_quad] == 0.0)
                        continue;                               // Source[idx_quad] = 0 -> not refine this point
                    else {
                        refine[idx_quad] = 1;                  // Source[idx_quad] != 0 -> refine this point
                        for ( unsigned id = 0; id < _neigh[idx_quad].size(); id ++)
                            refine[ _neigh[idx_quad][id]] = 1;     // also refine neighbours of this point
                    }
                }
            }
        }
        else {                                          // Same as above but for different sources for different energies
            if( _QextC[idx_energy][idx_cell].size() == 1u )
                if( _QextC[idx_energy][idx_cell][0] == 0.0) continue;
                else { refine = VectorU( _nqC, 1 ); break; }
            else{
                for (unsigned idx_quad = 0; idx_quad < _nqC; idx_quad++){
                    if( _QextC[idx_energy][idx_cell][idx_quad] == 0.0) continue;
                    else {
                        refine[idx_quad] = 1;
                        for ( unsigned id = 0; id < _neigh[idx_quad].size(); id ++)
                            refine[ _neigh[idx_quad][id]] = 1;     // also refine neighbours of this point
                    }
                }
            }
        }
    }

    for (unsigned id = 0; id < _nqC; id++){
        if ( refine[id] == 0 ){                            // NO refinement needed for quadpoint ID
            if ( _refine[id] >= _settings->GetRefineIter() ) { _refine[id] = 0;}            // reset if enough iterations
            else if ( _refine[id] != 0 ) { _refine[id] += 1;  }         // was also refined in last Iter; update _refine
        }
        else if ( refine[id] == 1 ){  _refine[id] = 1; }               // NEED refinement for quadpoint ID -- RESET counter to 1
        else {ErrorMessages::Error( "Refinement went wrong!", CURRENT_FUNCTION ); break;}
    }
}
*/

/*
void FirstCollisionCSDSNSolver::RefineQuadrature( ){

     if ( _refine == 0 ){   // no direction is refined -> use coarse quadrature
        _nq = _nqC;
        _quadPoints = _quadPointsC;
        _weights = _weightsC;
        _quadIDs = VectorU( _nq );
        for (unsigned id = 0; id < _nq; id++){ _quadIDs[ id ] = id;  }
    }
    else if ( _refine.nonZeros() == _refine.size() ){ // all directions are refined -> use fine quadrature
        _nq = _nqF;
        _quadPoints = _quadPointsF;
        _weights = _weightsF;
        _quadIDs = VectorU( _nq );
        for (unsigned id = 0; id < _nq; id++){ _quadIDs[ id ] = id + _nqC;  }
    }
    else{   // just some directions have to be coarsened or refined
        _nq = 0;
        _quadPoints = VectorVector( _nqF, Vector( 3 )); // reset calculation vectors : _nqF is maximum of possible quadPoints
        _weights = Vector( _nqF );
        _quadIDs = VectorU( _nqF );

        for (unsigned idx_quad = 0; idx_quad < _nqC; idx_quad++){
            if ( _refine[ idx_quad ] == 0 ){         // use coarse quadrature point
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
    }
}
*/

/*
void FirstCollisionCSDSNSolver::InterpolateSolution(){

    Matrix MatInterp( _nq, _nqOld, 0.0);
    std::vector<double> distance1( _nqOld );
    std::vector<unsigned> ids( _nqOld );
    std::cout << "NQOld: " << _nqOld << ", NQ: " << _nq << std::endl;

    // Setup Interpolation Matrix
    for ( unsigned iterNew = 0; iterNew < _nq; iterNew ++){
        unsigned check = 0;
        double ssum = 0;
        // Interpolation between three nearest values in old quadrature
        for( unsigned iterOld = 0; iterOld < _nqOld; iterOld++){
            // determine points that havent changed from last Iteration
            if ( _quadIDs[ iterNew ] == _quadIDsOld[ iterOld ] ){ MatInterp( iterNew, iterOld ) = 1.0;  goto next; }
            // calculate distance to every point in last iteration
            if ( _quadIDsOld[ iterOld ] < _nqC ){
                    distance1[ iterOld ] = norm( _quadPoints[ iterNew ] - _quadPointsC[ _quadIDsOld[ iterOld ] ]); }
            else    {
                distance1[ iterOld ] = norm( _quadPoints[ iterNew ] - _quadPointsF[ _quadIDsOld[ iterOld ] - _nqC ]); }
        }
        // sort distance vector
        std::iota( ids.begin(), ids.end(), 0);
        std::sort( ids.begin(), ids.end(), [&](int i,int j){return distance1[i] < distance1[j];} );
        //std::cout << "Iter: " << iterNew << "-- ids[0] = " << ids[0] << "-- ids[1] = " << ids[1]  << "-- ids[2] = " << ids[2]<< std::endl;
        if ( distance1[ ids [0] ] < 10e-10 ) { MatInterp( iterNew , ids [ 0 ]) = 1.0; } // if new point and old point are nearly the same
        else{
            for ( unsigned it = 0 ; it < 3 ; it++){ ssum +=  1.0 /  pow ( distance1[ ids[ it ] ] , 2 ); }   // sum over three nearest values
            for ( unsigned it = 0; it < 3; it++){   // setup Interpolation Matrix
                MatInterp( iterNew, ids[ it ] ) = 1.0 /  pow ( distance1[ ids[ it ] ] , 2 );
                MatInterp( iterNew, ids[ it ] ) /= ssum;
            }
        }
        next:{}
    }
    // Interpolate Solution
    for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
        _solFC[idx_cell] = MatInterp * _solFC[idx_cell];
    }
}

*/