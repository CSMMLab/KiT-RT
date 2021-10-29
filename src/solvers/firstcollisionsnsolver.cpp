//
// Created by chinsp on 29/10/21.
//

#include "../../include/solvers/firstcollisionsnsolver.h"

#include "common/config.h"
#include "common/io.h"
#include "common/mesh.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/textprocessingtoolbox.h"
#include "entropies/entropybase.h"
#include "fluxes/numericalflux.h"
#include "optimizers/optimizerbase.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/problembase.h"
#include "quadratures/quadraturebase.h"


// externals
#include "spdlog/spdlog.h"
#include <iostream>
#include <mpi.h>
#include <time.h>
#include <stdio.h>


FirstCollisionSNSolver::FirstCollisionSNSolver( Config * _settings ) : Solver ( _settings ){

    localizedIC = false; // Source only in positive x-direction

    /////// Variables for coarse SN Setting

    _nqC = _quadrature->GetNq();
    _quadPointsC = _quadrature->GetPoints();
    _weightsC    = _quadrature->GetWeights();
    _neigh = _quadrature->GetNeighbours();

    FILE * PTS = fopen((_settings->GetLogDir() + "/QuadPoints").c_str(), "a+");
    for( unsigned i = 0; i < _quadPointsC.size(); i++){
        fprintf(PTS, "%f\t%f\t%f\n", _quadPointsC[i][0], _quadPointsC[i][1], _quadPointsC[i][2]);
    }
    fclose(PTS);
    ////// Variables for fine SN Setting

    /// if Local Refinement is wanted
    if ( _settings->GetIsLocalRefine() == true){
        /// Refinement actual only works for Icosahedron Quadrature !!!
        switch (_settings->GetQuadName()){
            case QUAD_Icosahedron_Mid:
                _nqF = _quadrature->GetNqRefined();
                _quadPointsF = _quadrature->GetPointsRefined();
                _weightsF = _quadrature->GetWeightsRefined();
                _quadC2F = _quadrature->GetRefineVector();
                break;
            case QUAD_GaussLegendreTensorized:
                _settings->SetQuadOrder( 2 * _settings->GetQuadOrder() );
                _quadratureFine = QuadratureBase::Create( _settings ); // Gauss legendre needs settings class in construction
                _settings->SetQuadOrder( 1/2 * _settings->GetQuadOrder() );
                _nqF = _quadratureFine->GetNq();
                _quadPointsF = _quadratureFine->GetPoints();
                _quadPointsSphere = _quadrature->GetPointsSphere();
                _quadPointsSphereFine = _quadratureFine->GetPointsSphere();
                _weightsF = _quadratureFine->GetWeights();
                break;
            default: ErrorMessages::Error( "Local Refinement does not exist for chosen quadrature!", CURRENT_FUNCTION );
        }
    }
    else{ /// if NO Local Refinement wanted
        unsigned o;
        switch( _settings->GetQuadName() ){
            case QUAD_GaussLegendreTensorized:
                o = _settings->GetQuadOrder();
                _settings->SetQuadOrder( _settings->GetQuadOrderFine() );
                _quadratureFine = QuadratureBase::Create( _settings ); // Gauss legendre needs settings class in construction
                _settings->SetQuadOrder( o );
                break;
            default:
                _quadratureFine = QuadratureBase::Create( _settings->GetQuadName(), _settings->GetQuadOrderFine() );
        }
        _nqF = _quadratureFine->GetNq();
        _quadPointsF = _quadratureFine->GetPoints();
        _weightsF = _quadratureFine->GetWeights();
    }



    ///// External Source
    _QextC = _Q; // save external source in coarse quadrature (needed for local refinement)
    _settings->SetNQuadPoints( _nqF );
    _QextF = _problem->GetExternalSource( _energies );  // get external source in fine quadrature
    _settings->SetNQuadPoints( _nqC );

    ///// Set Initials

    if ( _settings->GetIsLocalRefine() == true){     /// if Local Refinement is wanted
        _nq = _nqC;
        _quadPoints = _quadPointsC;
        _weights = _weightsC;
        _refine = VectorU( _nqC, 0 );
        _quadIDs = VectorU( _nqC, 0);
        for ( unsigned id = 0; id < _nqC; id ++) { _quadIDs[ id ] = id; }
        _Q = _QextC;
    }
    else if( _settings->GetIsArtificialScattering() == true ){
        _nq = _nqF;
        _quadPoints = _quadPointsF;
        _weights = _weightsF;
        _Q = _QextF;
    }
    else{
        _nq = _nqF;
        _quadPoints = _quadPointsF;
        _weights = _weightsF;
        _Q = _QextF;
    }
    if ( _settings->GetIsLocalRefine() == true ){
        RefineTest();
    }

    std::cout << "NQ: " << _nqC << " , " << _nq << std::endl;
    std::cout << "NCells: " << _nCells << std::endl;

    ///// Scattering Kernel
    KERNEL = ScatteringKernel::CreateScatteringKernel( _settings->GetKernelName(), _quadrature );
    _scatteringKernel   = KERNEL->GetScatteringKernel();                                   // Scattering Kernel in coarse mesh
    _scatteringKernelFC = KERNEL->GetScatteringKernelFirstCollision( _nq, _weights);       // Scattering Kernel from Fine to Coarse mesh
    // ScatteringKernel * K = ScatteringKernel::CreateScatteringKernel( _settings->GetKernelName(), _quadratureFine );
    // _scK = K->GetScatteringKernel();

    ///// Initial condition ( completely in uncollided part )
    _sol = VectorVector( _nCells, Vector( _nqC, 0.0 ));
    //_sol = _problem->SetupIC(); // This is in coarse quadrature
    _settings->SetNQuadPoints(_nq );
    _solFC = _problem->SetupIC();
    _settings->SetNQuadPoints(_nqC);

    if ( localizedIC ){
        for( unsigned i = 0; i < _nq; i++ ){
            if ( _quadPoints[i][0] < 0.0 || fabs(_quadPoints[i][1]) > 0.5  || fabs(_quadPoints[i][2] > 0.5) ){ // um x-Achse
                for( unsigned j = 0; j < _nCells; j++ )
                    _solFC[j][i] = 0.0;
            }
        }
    }
    ///// First Collision Source initialize
    _QFirstCollision = VectorVector( _nCells, Vector( _nqC ));

    std::cout << "Check Build" << std::endl;
}

void FirstCollisionSNSolver::Solve(){

    FILE * Times = fopen((_settings->GetLogDir()+"/TimeFile").c_str(), "a+");
    clock_t start, end, startiter, enditer,t, t1, t2, t3, t4, t5, t6, t7, t8;

    start = clock();
    PrepareVolumeOutput();
    DrawPreSolverOutput();

    for ( unsigned idx_energy = 0; idx_energy < _nEnergies; idx_energy++){
        startiter = clock();
        /*
        if ( _settings->GetIsLocalRefine() == true ){
            IterRefinement( idx_energy);
            std::cout << "Check Local Refinement; NQ = " << _nq << std::endl;
        }
        */
        t = clock();
        IterPreprocessing( idx_energy );
        t1 = clock() -t;

        t = clock();
        FluxUpdateUncollided();
        t2 = clock() - t;

        t = clock();
        FVMUpdateUncollided( idx_energy );
        t3 = clock() - t;

        if ( _settings->GetIsArtificialScattering() == true ){
            t = clock();
            AddArtificialScattering( idx_energy );
            t4 = clock() - t;
        }

        t = clock();
        ComputeFirstCollisionSource( idx_energy );
        t8 = clock() - t;

        t = clock();
        FluxUpdateCollided( );
        t5 = clock() -t;

        t = clock();
        FVMUpdateCollided( idx_energy );
        t6 = clock() - t;

        /*
        if ( _settings->GetIsArtificialScattering() == true ){
            t = clock();
            AddArtificialScatteringCoarse( idx_energy );
            t4 = clock() - t;
        }
        */
        t = clock();
        IterPostprocessing( );

        WriteVolumeOutput( idx_energy );
        WriteScalarOutput( idx_energy );
        PrintScreenOutput( idx_energy );
        PrintHistoryOutput( idx_energy );
        PrintVolumeOutput( idx_energy );

        t7 = clock() - t;

        enditer = clock() - startiter;

        fprintf(Times, "PreProcessing: \t %f \n", t1 / (float)CLOCKS_PER_SEC);
        fprintf(Times, "Flux Update U: \t %f \n", t2 / (float)CLOCKS_PER_SEC);
        fprintf(Times, "FVM  Update U: \t %f \n", t3 / (float)CLOCKS_PER_SEC);
        fprintf(Times, "Artfic Scattr: \t %f \n", t4 / (float)CLOCKS_PER_SEC);
        fprintf(Times, " FC Source   : \t %f \n", t8 / (float)CLOCKS_PER_SEC);
        fprintf(Times, "Flux Update C: \t %f \n", t5 / (float)CLOCKS_PER_SEC);
        fprintf(Times, "FVM  Update C: \t %f \n", t6 / (float)CLOCKS_PER_SEC);
        fprintf(Times, "PostProcessin: \t %f \n", t7 / (float)CLOCKS_PER_SEC);

        fprintf(Times, "\n Iteration %u: \t %f \n \n", idx_energy, enditer / (float)CLOCKS_PER_SEC);

    }

    // Post Solver Processing
    DrawPostSolverOutput();

    end = clock() - start;

    fprintf(Times, "\n Complete Run: \t %f \n \n", end / (float)CLOCKS_PER_SEC);
    fprintf(Times, "\n per Iter: \t %f \n \n", end / (float)CLOCKS_PER_SEC / (float)_nEnergies );
}


////// --------- Solver Functions -------- /////

void FirstCollisionSNSolver::IterRefinement( unsigned idx_energy ){
    std::cout << "Entered Refinement" << std::endl;
    // --- Step 0 --- Save Values from last Iteration
    _nqOld = _nq;
    _quadIDsOld = _quadIDs;
    VectorU _refineOld = _refine;
    //_refine.reset();
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

    // --- Step 2.3 --- Computation of new scattering Kernel
    _scatteringKernelFC = KERNEL->GetScatteringKernelFirstCollision( _nq, _weights );
    std::cout << "Check: Scattering Kernel calculated" << std::endl;

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

void FirstCollisionSNSolver::IterPreprocessing( unsigned idx_energy ){

    if ( _settings -> GetIsLocalRefine() == true){
        /// Adapt external source to correct size
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
                        else    { _Q[idx_energy][idx_cell][ id ] = _QextF[idx_energy][idx_cell][ _quadIDs[id] - _nqC ]; }
                    }
                }
            }
        }
    }
    // Source Term !!!!!
    if ( idx_energy == 0){
        for (unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++){
            if( _Q.size() == 1u ) {                   // constant source for all energies
                if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                    _solFC[idx_cell] += _Q[0][idx_cell][0];
                else
                    _solFC[idx_cell] += _Q[0][idx_cell];
            }
            else {
                if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                    _solFC[idx_cell] += _Q[idx_energy][idx_cell][0];
                else
                    _solFC[idx_cell] += _Q[idx_energy][idx_cell];
            }
        }
    }
}

void FirstCollisionSNSolver::FluxUpdateUncollided(){

    _psiDx = VectorVector( _nCells, Vector( _nq, 0.0 ) );
    _psiDy = VectorVector( _nCells, Vector( _nq, 0.0 ) );
    _solNewFC = _solFC;
    if( _reconsOrder > 1 ) {
        _mesh->ReconstructSlopesU( _nq, _psiDx, _psiDy, _solFC );    // unstructured reconstruction
    }

    double psiL, psiR;
// #pragma omp parallel for
    for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){

        // Dirichlet Boundaries stay at IC
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) { continue; }

        // _solNewFC[idx_cell] = ConstructFluxSN( idx_cell, true );
        for( unsigned idx_quad = 0; idx_quad < _nq; ++idx_quad ) {
            // Reset temporary variable
            _solNewFC[idx_cell][idx_quad] = 0;
            // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
            for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); ++idx_neighbor ) {
                // store flux contribution on psiNew_sigmaS to save memory
                if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells ){
                    _solNewFC[idx_cell][idx_quad] +=
                            _g->Flux( _quadPoints[idx_quad], _solFC[idx_cell][idx_quad], _solFC[idx_cell][idx_quad], _normals[idx_cell][idx_neighbor] );
                }
                else {
                    switch( _reconsOrder ) {
                        // first order solver
                        case 1:
                            _solNewFC[idx_cell][idx_quad] += _g->Flux( _quadPoints[idx_quad],
                                                                       _solFC[idx_cell][idx_quad],
                                                                       _solFC[_neighbors[idx_cell][idx_neighbor]][idx_quad],
                                                                       _normals[idx_cell][idx_neighbor] );
                            break;
                            // second order solver
                        case 2:
                            // left status of interface
                            psiL = _solFC[idx_cell][idx_quad] +
                                   _psiDx[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[idx_cell][0] ) +
                                   _psiDy[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[idx_cell][1] );
                            // right status of interface
                            psiR = _solFC[_neighbors[idx_cell][idx_neighbor]][idx_quad] +
                                   _psiDx[_neighbors[idx_cell][idx_neighbor]][idx_quad] *
                                   ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[_neighbors[idx_cell][idx_neighbor]][0] ) +
                                   _psiDy[_neighbors[idx_cell][idx_neighbor]][idx_quad] *
                                   ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[_neighbors[idx_cell][idx_neighbor]][1] );
                            // positivity check (if not satisfied, deduce to first order)
                            if( psiL < 0.0 || psiR < 0.0 ) {
                                psiL = _solFC[idx_cell][idx_quad];
                                psiR = _solFC[_neighbors[idx_cell][idx_neighbor]][idx_quad];
                            }
                            // flux evaluation
                            _solNewFC[idx_cell][idx_quad] += _g->Flux( _quadPoints[idx_quad], psiL, psiR, _normals[idx_cell][idx_neighbor] );
                            break;
                            // higher order solver
                        case 3:
                            std::cout << "higher order is WIP" << std::endl;
                            break;
                            // default: first order solver
                        default:
                            _solNewFC[idx_cell][idx_quad] += _g->Flux( _quadPoints[idx_quad],
                                                                       _solFC[idx_cell][idx_quad],
                                                                       _solFC[_neighbors[idx_cell][idx_neighbor]][idx_quad],
                                                                       _normals[idx_cell][idx_neighbor] );
                    }
                }
            }
        }
    }
}


void FirstCollisionSNSolver::FVMUpdateUncollided( unsigned idx_energy ){

// #pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;

        // Flux Update and Absorption Term
        // loop over all ordinates
        for( unsigned idx_quad = 0; idx_quad < _nq; ++idx_quad ) {
            // time update angular flux with numerical flux and total scattering cross section
            _solNewFC[idx_cell][idx_quad] = _solFC[idx_cell][idx_quad] - ( _dE / _areas[idx_cell] ) * _solNewFC[idx_cell][idx_quad]
                                            - _dE * _sigmaT[idx_energy][idx_cell] * _solFC[idx_cell][idx_quad];
        }
        //_solFC[ idx_cell ] += _dE * _sigmaS[ idx_energy ][ idx_cell ] * _scK * _solFC[idx_cell];
    }
}

void FirstCollisionSNSolver::ComputeFirstCollisionSource( unsigned idx_energy ){

    // Q = _sigmaS * int_S^2 ( Sigma * psi_u ) d Omega'
    // isotropic scattering -> Sigma independent of Omega'
    FILE * FCSource = fopen((_settings->GetLogDir() + "/FCSource.txt").c_str(), "a+");
    // std::cout << "FC Kernel size: (" << _scatteringKernelFC.rows() << "," << _scatteringKernelFC.columns() << ")" << std::endl;
    // std::cout << "solNewFC size: (" << _solNewFC[0].size() << ")" << std::endl;
    // std::cout << "QFirstCollision size: (" << _QFirstCollision.size() << " , " << _QFirstCollision[0].size() << ")" << std::endl;

    for (unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
        _QFirstCollision[ idx_cell ] = _sigmaS[ idx_energy ][ idx_cell ] * _scatteringKernelFC * _solNewFC[idx_cell];

        for( unsigned idx_quad = 0; idx_quad < _nqC; idx_quad ++)
            if ( _QFirstCollision[idx_cell][idx_quad] > 10e10){
                ErrorMessages::Error("infty in _FC Source","FCSource");
                fprintf(FCSource, "FCSource(%u,%u) = \t %f \n", idx_cell, idx_quad, _QFirstCollision[idx_cell][idx_quad] );
            }

    }


    // following only possible for isotropic scattering (if anisotropic we need different scattering calculation !!!! )
    // aim: subtract scattered particles from uncollided flux
    /*
    Matrix SCKernelFine(_nqF, _nqF);
    for ( unsigned id = 0; id < _nqF; id ++)
        for ( unsigned jd = 0; jd < _nqF; jd ++)
            SCKernelFine(id, jd) = _scatteringKernelFC(0, jd);
    for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++)
        _solNewF[ idx_cell ] -= _sigmaS[ idx_energy ][ idx_cell ] * SCKernelFine * _solNewF[ idx_cell ];
    */
}


void FirstCollisionSNSolver::FluxUpdateCollided( ){

    _psiDx = VectorVector( _nCells, Vector( _nqC, 0.0 ) );
    _psiDy = VectorVector( _nCells, Vector( _nqC, 0.0 ) );
    _solNew = _sol;
    if( _reconsOrder > 1 ) {
        _mesh->ReconstructSlopesU( _nqC, _psiDx, _psiDy, _sol );    // unstructured reconstruction
    }
    double psiL, psiR;
#pragma omp parallel for
    for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
        // Dirichlet Boundaries stayd
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // Loop over all ordinates
        for( unsigned idx_quad = 0; idx_quad < _nqC; ++idx_quad ) {
            // Reset temporary variable
            _solNew[idx_cell][idx_quad] = 0;
            // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
            for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); ++idx_neighbor ) {
                // store flux contribution on psiNew_sigmaS to save memory
                if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                    _solNew[idx_cell][idx_quad] +=
                            _g->Flux( _quadPointsC[idx_quad], _sol[idx_cell][idx_quad], _sol[idx_cell][idx_quad], _normals[idx_cell][idx_neighbor] );
                else {
                    switch( _reconsOrder ) {
                        // first order solver
                        case 1:
                            _solNew[idx_cell][idx_quad] += _g->Flux( _quadPointsC[idx_quad],
                                                                     _sol[idx_cell][idx_quad],
                                                                     _sol[_neighbors[idx_cell][idx_neighbor]][idx_quad],
                                                                     _normals[idx_cell][idx_neighbor] );
                            break;
                            // second order solver
                        case 2:
                            // left status of interface
                            psiL = _sol[idx_cell][idx_quad] +
                                   _psiDx[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[idx_cell][0] ) +
                                   _psiDy[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[idx_cell][1] );
                            // right status of interface
                            psiR = _sol[_neighbors[idx_cell][idx_neighbor]][idx_quad] +
                                   _psiDx[_neighbors[idx_cell][idx_neighbor]][idx_quad] *
                                   ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[_neighbors[idx_cell][idx_neighbor]][0] ) +
                                   _psiDy[_neighbors[idx_cell][idx_neighbor]][idx_quad] *
                                   ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[_neighbors[idx_cell][idx_neighbor]][1] );
                            // positivity check (if not satisfied, deduce to first order)
                            if( psiL < 0.0 || psiR < 0.0 ) {
                                psiL = _sol[idx_cell][idx_quad];
                                psiR = _sol[_neighbors[idx_cell][idx_neighbor]][idx_quad];
                            }
                            // flux evaluation
                            _solNew[idx_cell][idx_quad] += _g->Flux( _quadPointsC[idx_quad], psiL, psiR, _normals[idx_cell][idx_neighbor] );
                            break;
                            // higher order solver
                        case 3:
                            std::cout << "higher order is WIP" << std::endl;
                            break;
                            // default: first order solver
                        default:
                            _solNew[idx_cell][idx_quad] += _g->Flux( _quadPointsC[idx_quad],
                                                                     _sol[idx_cell][idx_quad],
                                                                     _sol[_neighbors[idx_cell][idx_neighbor]][idx_quad],
                                                                     _normals[idx_cell][idx_neighbor] );
                    }
                }
            }
        }
    }
}

void FirstCollisionSNSolver::FVMUpdateCollided( unsigned idx_energy ){

    // std::cout << "scattering kernel size: " << _scatteringKernel.rows() << " , " << _scatteringKernel.columns() << std::endl;
    // std::cout << "sol size: " << _sol.size() << " , " << _sol[0].size() << std::endl;
    // std::cout << "solNew size: " << _solNew.size() << " , " << _solNew[0].size() << std::endl;

#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // loop over all ordinates
        for( unsigned idx_quad = 0; idx_quad < _nqC; ++idx_quad ) {
            // time update angular flux with numerical flux and total scattering cross section
            _solNew[ idx_cell ][ idx_quad ] = _sol[idx_cell][idx_quad] - ( _dE / _areas[idx_cell] ) * _solNew[idx_cell][idx_quad] -
                                              _dE * _sigmaT[idx_energy][idx_cell] * _sol[idx_cell][idx_quad];
        }
        // compute scattering effects
        _solNew[ idx_cell ] += _dE * _sigmaS[idx_energy][idx_cell] * _scatteringKernel * _sol[idx_cell];    // multiply scattering matrix with psi

        // Add First Collision Source
        _solNew[ idx_cell ] += _dE * _sigmaS[ idx_energy ][ idx_cell ] * _scatteringKernelFC * _solNewFC[idx_cell];
    }
}

void FirstCollisionSNSolver::IterPostprocessing( ){

    _sol = _solNew;
    _solFC = _solNewFC; // was set in FV Update Uncollided
    ComputeRadFlux();
}


//// ------------ Helper Functions for SN ---------- ////

void FirstCollisionSNSolver::ComputeRadFlux(){

    FILE * radFluxfile = fopen((_settings->GetLogDir() + "/radFluxFC.txt").c_str(), "a+");
    fprintf(radFluxfile, "\n Iteration: \n ");
    fprintf(radFluxfile, "\t uncoll \t coll \t \t cumulative \n \n ");
    double f1 = 0, f2 = 0, f3 = 0;
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        // flux from uncollided
        double fluxUC;
        fluxUC = blaze::dot( _solFC[idx_cell], _weights );
        // flux from collided
        double fluxC = blaze::dot( _sol[idx_cell], _weightsC );
        // cumulative flux
        _fluxNew[idx_cell] = fluxUC + fluxC;

        f1 += fluxUC;
        f2 += fluxC;
        f3 += _fluxNew[idx_cell];
        //_fluxNew[idx_cell] = blaze::dot( _sol[idx_cell], _weights );
        //fprintf(radFluxfile, "\t %f \n", (float) _fluxNew[idx_cell] * 1000);
    }

    fprintf(radFluxfile, "\t %f \t %f \t %f \n", (float) f1 , (float) f2 , (float) f3 );
    fclose(radFluxfile);
}


//// ------------ Helper Functions for Refinement ---------- ////


void FirstCollisionSNSolver::DetermineRefinement( unsigned idx_energy ){

    VectorU refine( _nqC, 0); // --> refine in {0,1} -> need refinement in this iteration or not

    // loop over COARSE external Source ( if _Qext[idx_energy][..][idx_quad] != 0 have to refine in direction idx_quad )
    for (unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++){
        if( _QextC.size() == 1u ) {                     // constant source for all energies
            if( _QextC[0][idx_cell].size() == 1u ) {         // isotropic source
                if( _QextC[0][idx_cell][0] == 0.0) continue;        // Isotropic source is 0 -> nothing to do
                else { refine = VectorU( _nqC, 1 ); break; }        // Isotropic Source is not 0 -> Refine whole Quadrature
            }
            else{   // anisotropic source
                for (unsigned idx_quad = 0; idx_quad < _nqC; idx_quad++){
                    if( _QextC[0][idx_cell][idx_quad] == 0.0) { continue; }  // Source[idx_quad] = 0 -> not refine this point
                    refine[idx_quad] = 1;                  // Source[idx_quad] != 0 -> refine this point
                    switch (_settings->GetQuadName() ){
                        case QUAD_Icosahedron_Mid:      // also refine neighbours of this point
                            for ( unsigned id = 0; id < _neigh[idx_quad].size(); id ++) { refine[ _neigh[idx_quad][id] ] = 1; }
                    }
                }
            }
        }
        else {                                          // Same as above but for different sources for different energies
            if( _QextC[idx_energy][idx_cell].size() == 1u ){
                if( _QextC[idx_energy][idx_cell][0] == 0.0) continue;
                else { refine = VectorU( _nqC, 1 ); break; }
            }
            else{
                for (unsigned idx_quad = 0; idx_quad < _nqC; idx_quad++){
                    if( _QextC[idx_energy][idx_cell][idx_quad] == 0.0) { continue; }
                    else{
                        refine[idx_quad] = 1;
                        switch (_settings->GetQuadName() ){
                            case QUAD_Icosahedron_Mid:             // also refine neighbours of this point
                                for ( unsigned id = 0; id < _neigh[idx_quad].size(); id ++) { refine[ _neigh[idx_quad][id]] = 1; }
                        }
                    }
                }
            }
        }
    }
    std::cout << "RefIter: " << _settings->GetRefineIter() << std::endl;
    for (unsigned id = 0; id < _nqC; id++){
        if ( refine[id] == 0 ){                            // NO refinement needed for quadpoint ID
            if ( _refine[id] >= _settings->GetRefineIter() ) { _refine[id] = 0;}            // reset if enough iterations
            else if ( _refine[id] != 0 ) { _refine[id] += 1;  }         // was also refined in last Iter; update _refine
        }
        else if ( refine[id] == 1 ){ _refine[id] = 1; }               // NEED refinement for quadpoint ID -- RESET counter to 1
        else {ErrorMessages::Error( "Refinement went wrong!", CURRENT_FUNCTION ); break;}
    }


    FILE * REF = fopen(( _settings->GetLogDir() + "/Refine Vectors").c_str(), "a+");
    fprintf(REF, "Iteration %u \n\n", idx_energy);
    for ( unsigned i = 0; i < _nqC; i++){
        if( idx_energy == 0 )
            fprintf( REF, "ID: %u:    \t Refine = %u \t _refine = %u \n", i, refine[i], _refine[i]);
    }
    fclose(REF);
}


void FirstCollisionSNSolver::RefineQuadrature( ){

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

        switch ( _settings->GetQuadName()){
            case QUAD_Icosahedron_Mid:
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
                            _quadIDs[ _nq + id ] = _nqC + _quadC2F[idx_quad][id];
                        }
                        _nq += _quadC2F[idx_quad].size();
                    }
                }
                break;
            case QUAD_GaussLegendreTensorized:
                unsigned nn = _quadrature->GetOrder();
                std::cout << "NN: " << nn << std::endl;
                Vector muVec( nn );
                Vector phiVec( 2 * nn );
                for ( unsigned i = 0; i < muVec.size(); i ++){ muVec[i] = _quadPointsSphere[ i * phiVec.size() ][0]; }
                for ( unsigned i = 0; i < phiVec.size(); i ++){ phiVec[i] = _quadPointsSphere[i][1]; }
                std::cout << "Check1" << std::endl;
                double muMax = -1.0;
                double muMin = 1.0;
                double phiMax = 0.0;
                double phiMin = 3.0;
                unsigned muMinID, muMaxID, phiMinID, phiMaxID;
                std::cout << _refine.nonZeros() << " , " << _nqC << std::endl;
                for (unsigned idx_quad = 0; idx_quad < _nqC; idx_quad++){
                    if ( _refine[idx_quad] != 0 ){
                        if (_quadPointsSphere[idx_quad][0] < muMin ) { muMin = _quadPointsSphere[idx_quad][0] ; muMinID = (unsigned) (idx_quad / phiVec.size()); }
                        if (_quadPointsSphere[idx_quad][0] > muMax ) { muMax = _quadPointsSphere[idx_quad][0] ; muMaxID = (unsigned) (idx_quad / phiVec.size()); }
                        if (_quadPointsSphere[idx_quad][1] < phiMin ) { phiMin = _quadPointsSphere[idx_quad][1] ; phiMinID = idx_quad % phiVec.size(); }
                        if (_quadPointsSphere[idx_quad][1] > phiMax ) { phiMax = _quadPointsSphere[idx_quad][1] ; phiMaxID = idx_quad % phiVec.size(); }
                    }
                }
                std::cout << "MU MIN = " << muMinID << ",  PHI MIN = " << phiMinID << "  MU MAX = " << muMaxID << ",  PHI MAX = " << phiMaxID << std::endl;
                std::cout << "MU MIN = " << muMin << ",  PHI MIN = " << phiMin << "  MU MAX = " << muMax << ",  PHI MAX = " << phiMax << std::endl;

                VectorU rref(_nqC, 1u);
                double weightsum( phiVec.size() );
                unsigned ind = 0;
                for (unsigned idx_quad = 0; idx_quad < _nqC; idx_quad++){
                    if ( _quadPointsSphere[idx_quad][0] < muMin || _quadPointsSphere[idx_quad][0] > muMax ){ rref[idx_quad] = 0; }
                    if ( _quadPointsSphere[idx_quad][1] < phiMin || _quadPointsSphere[idx_quad][1] > phiMax ){ rref[idx_quad] = 0; }
                    if ( rref[idx_quad] == 1){
                        for( unsigned jd = 0; jd < phiVec.size(); jd++)
                            if (_quadPointsSphere[idx_quad][1] == phiVec[jd] ) { weightsum += _weightsC[idx_quad]; ind += 1; }
                        weightsum /= ind;
                    }
                }
                std::cout << "Check 3" << std::endl;
                // alle points, die jetzt noch rref = 1 haben, werden entsorgt
                // calculate mu points for refined area
                QuadratureBase * quadRef = QuadratureBase::Create( QUAD_GaussLegendre1D, 2 * ( muMaxID - muMinID + 2 ));
                unsigned nq = quadRef->GetNq();
                VectorVector muV = quadRef->GetPoints();
                Vector weights = quadRef->GetWeights();
                Vector mu(nq);
                // adapt to area
                if (muMinID == 0){ muMin = -1;} else { muMin = muVec[muMinID - 1]; }
                if (muMaxID == muVec.size()-1){ muMax = 1;} else { muMax = muVec[muMaxID + 1]; }
                for ( unsigned id = 0; id < nq; id++){
                    mu[id] = ( muMax - muMin ) / 2 * muV[id][0] + ( muMax +muMin) / 2;
                    weights[id] /= weightsum;
                }
                std::cout << "Check 4" << std::endl;
                // calculate phi points
                Vector phi( 2 * nq );
                if ( phiMinID == 0 ){ phiMin = phiVec[phiVec.size() - 1] - 2* M_PI; } else { phiMin = phiVec[ phiMinID - 1] ; }
                if ( phiMaxID == phiVec.size()-1 ){ phiMax = phiVec[0] + 2*M_PI; } else { phiMax = phiVec[ phiMaxID + 1]; }
                for( unsigned i = 0; i < 2 * nq; ++i ) {
                    phi[i] = ( i + 0.5 ) * M_PI / nq;
                    phi[i] = ( phiMax - phiMin ) / 2 * phi[i] + phiMin; // adapt to area
                }
                std::cout << "Check 5" << std::endl;
                // construct whole quadrature
                ind = 0;
                for( unsigned id = 0; id < _nqC; id++){
                    if( rref[id] == 0 ) {
                        _quadPoints[ ind ] = _quadPointsC[ id ];
                        _weights[ ind ] = _weightsC[ id ];
                        _quadIDs[ ind ] = id;
                        ind +=1;
                    }
                }
                std::cout << "Check 6" << std::endl;
                for ( unsigned i = 0; i < mu.size(); i ++){
                    for( unsigned j = 0; j < phi.size(); j ++){
                        _quadPoints[ ind ][0] = sqrt( 1 - mu[j] * mu[j] ) * std::cos( phi[i] );
                        _quadPoints[ ind ][1] = sqrt( 1 - mu[j] * mu[j] ) * std::sin( phi[i] );
                        _quadPoints[ ind ][2] = mu[j];
                        _weights[ ind ] = M_PI / mu.size() * weights[j];
                        ind += 1;
                        _quadIDs[ ind ] = i + _nqC;
                    }
                }
                _nq = ind;
                std::cout << "Check 7" << std::endl;
        }
        // truncate vectors to _nq
        _quadPoints.resize( _nq );
        _weights.resize( _nq );
        _quadIDs.resize( _nq );
    }
}

void FirstCollisionSNSolver::InterpolateSolution(){

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

void FirstCollisionSNSolver::FluxUpdate(){}             // not used
void FirstCollisionSNSolver::FVMUpdate( unsigned ){}    // not used


// -------- Output generators ( copied from SN_Solver ) -------------

void FirstCollisionSNSolver::PrepareVolumeOutput() {
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

            case ANALYTIC:
                // one entry per cell
                _outputFields[idx_group].resize( 1 );
                _outputFieldNames[idx_group].resize( 1 );
                _outputFields[idx_group][0].resize( _nCells );
                _outputFieldNames[idx_group][0] = std::string( "analytic radiation flux density" );
                break;

            default: ErrorMessages::Error( "Volume Output Group not defined for First Collision SN Solver!", CURRENT_FUNCTION ); break;
        }
    }
}


void FirstCollisionSNSolver::WriteVolumeOutput( unsigned idx_pseudoTime ) {
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();
    FILE * Output = fopen((_settings->GetLogDir() + "Output").c_str(), "a+");
    FILE * OutputAnalytic = fopen((_settings->GetLogDir() + "OutputAnalytic").c_str(), "a+");
    if( ( _settings->GetVolumeOutputFrequency() != 0 && idx_pseudoTime % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
        ( idx_pseudoTime == _nEnergies - 1 ) ) { // need sol at last iteration

        for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
            switch( _settings->GetVolumeOutput()[idx_group] ) {
                case MINIMAL:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = _fluxNew[idx_cell];
                        fprintf(Output, "%f\t%f\t%f\t%u\n", _fluxNew[idx_cell], _mesh->GetCellMidPoints()[idx_cell][0], _mesh->GetCellMidPoints()[idx_cell][1], idx_cell);
                    }
                    break;

                case ANALYTIC:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        double time                           = idx_pseudoTime * _dE;
                        _outputFields[idx_group][0][idx_cell] = _problem->GetAnalyticalSolution(
                                _mesh->GetCellMidPoints()[idx_cell][0], _mesh->GetCellMidPoints()[idx_cell][1], time, _sigmaS[idx_pseudoTime][idx_cell] );
                        fprintf(OutputAnalytic, "%f\t%f\t%f\t%u\n", _outputFields[idx_group][0][idx_cell], _cellMidPoints[idx_cell][0], _cellMidPoints[idx_cell][1], idx_cell);
                    }
                    break;

                default: ErrorMessages::Error( "Volume Output Group not defined for First Collision SN Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
    fclose(Output);
    fclose(OutputAnalytic);
}



void FirstCollisionSNSolver::AddArtificialScattering( unsigned iter ) {

    Matrix ASKernel(_nq, _nq, 0.0);
    double betaAS = _settings->GetBetaAS();
    double sigmaAS = _settings->GetSigmaAS();

    // Calculate Artificial Scattering Kernel
    for ( unsigned id = 0; id < _nq; id ++){
        double sum = 0;
        for ( unsigned jd = 0; jd < _nq; jd ++){
            double mu = 0;
            for ( unsigned i = 0; i < _quadPoints[0].size(); i++){
                mu += _quadPoints[id][i] * _quadPoints[jd][i];
            }
            double expTerm = ( - ( ( mu - 1 ) * ( mu - 1 ) ) * _nq * _nq ) / ( betaAS * betaAS );
            ASKernel( id, jd ) = _nq / betaAS * exp( expTerm );
            sum += ASKernel( id, jd );
        }
        for ( unsigned jd = 0; jd < _nq; jd++){
            ASKernel( id , jd ) /= sum ;        // Normalize
        }
    }
    // Add Artificiel Scattering Terms to actual solution
    for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
        _solNewFC[idx_cell] += _dE * sigmaAS * ( ASKernel * _solFC[idx_cell] - _solFC[idx_cell] );
    }
}

void FirstCollisionSNSolver::AddArtificialScatteringCoarse( unsigned iter ) {

    Matrix ASKernel(_nqC, _nqC, 0.0);
    double betaAS = _settings->GetBetaAS();
    double sigmaAS = _settings->GetSigmaAS();

    // Calculate Artificial Scattering Kernel
    for ( unsigned id = 0; id < _nqC; id ++){
        double sum = 0;
        for ( unsigned jd = 0; jd < _nqC; jd ++){
            double mu = 0;
            for ( unsigned i = 0; i < _quadPointsC[0].size(); i++){
                mu += _quadPointsC[id][i] * _quadPointsC[jd][i];
            }
            double expTerm = ( - ( ( mu - 1 ) * ( mu - 1 ) ) * _nqC * _nqC ) / ( betaAS * betaAS );
            ASKernel( id, jd ) = _nqC / betaAS * exp( expTerm );
            sum += ASKernel( id, jd );
        }
        for ( unsigned jd = 0; jd < _nqC; jd++){
            ASKernel( id , jd )/= sum ;        // Normalize
        }
    }
    // Add Artificiel Scattering Terms to actual solution
    for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
        _solNew[idx_cell] += _dE * sigmaAS * ( ASKernel * _sol[idx_cell] - _sol[idx_cell] );
    }
}

void FirstCollisionSNSolver::RefineTest(){

    unsigned k = 0;
    _quadPoints = VectorVector( _nqF, Vector(3, 0.0 ));
    _weights = Vector( _nqF);
    _quadIDs = VectorU( _nqF);

    for( unsigned id = 0; id < _nqC; id ++){
        if ( localizedIC ){
            if ( _quadPointsC[id][0] < 0.0 || fabs(_quadPointsC[id][1]) > 0.7 || fabs(_quadPointsC[id][2]) > 0.7 ){ // um x-Achse
                _quadPoints[ k ] = _quadPointsC[ id ];
                _weights[ k ] = _weightsC[ id ];
                _quadIDs[ k ] = id;
                k += 1;
            }
            else{
                for( unsigned j = 0; j < _quadC2F[id].size(); j++ ){
                    _quadPoints[ k ]= _quadPointsF[ _quadC2F[id][j] ];
                    _weights[ k ] = _weightsF[ _quadC2F[ id ][j]];
                    _quadIDs[ k ] = _quadC2F[id][j] + _nqC;
                    k += 1;
                }
            }
        }
        else{
            if( _quadPointsC[id][0] < 0.0 || _quadPointsC[id][1] < 0.0 ){ // erster quadrant
                _quadPoints[ k ] = _quadPointsC[ id ];
                _weights[ k ] = _weightsC[ id ];
                _quadIDs[ k ] = id;
                k += 1;
            }
            else{
                for( unsigned j = 0; j < _quadC2F[id].size(); j++ ){
                    _quadPoints[ k ]= _quadPointsF[ _quadC2F[id][j] ];
                    _weights[ k ] = _weightsF[ _quadC2F[ id ][j]];
                    _quadIDs[ k ] = _quadC2F[id][j] + _nqC;
                    k += 1;
                }
            }
        }
    }

    _nq = k;
    _quadPoints.resize(_nq);
    _weights.resize(_nq);
    _quadIDs.resize(_nq);
}