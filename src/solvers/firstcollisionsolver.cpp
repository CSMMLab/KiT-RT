//
// Created by chinsp on 29/10/21.
//

#include "../../include/solvers/firstcollisionsolver.h"

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
#include "toolboxes/sphericalbase.h"

// externals
#include "spdlog/spdlog.h"
#include <iostream>
#include <mpi.h>
#include <time.h>
#include <stdio.h>


FirstCollisionSolver::FirstCollisionSolver( Config * settings ) : SolverBase( settings ){

    /////// Variables for coarse Setting
    _nq = _quadrature->GetNq();
    _quadPoints = _quadrature->GetPoints();
    _weights    = _quadrature->GetWeights();
    _quadPointsSphere = _quadrature->GetPointsSphere();
    ScatteringKernel* k = ScatteringKernel::CreateScatteringKernel( _settings->GetKernelName(), _quadrature );
    _scatteringKernel   = k->GetScatteringKernel();                             // Scattering Kernel in coarse mesh

    ////// Variables for fine SN Setting

    _nqF = _quadrature->GetNqRefined();
    _quadPointsF = _quadrature->GetPointsRefined();
    _weightsF = _quadrature->GetWeightsRefined();

    ////// MN Variables
    _LMaxDegree    = _settings->GetMaxMomentDegree();                   // max degree of spherical basis
    _basis         = SphericalBase::Create( _settings );                // Initialize Spherical Basis
    _nTotalEntries = _basis->GetBasisSize();                            // Basis Size

    _scatterMatDiag = Vector( _nTotalEntries, 0.0 );                    // Scattering Kernel Matrix (saved as Vector bc of diagonal)
    _scatterMatDiag[ 0 ] = -1; /// why is this -1 ????

    _entropy = EntropyBase::Create( _settings );    	                // Initialize Entropy
    _optimizer = OptimizerBase::Create( _settings );                    // Initialize Optimizer
    _alpha = VectorVector( _nCells, Vector( _nTotalEntries, 0.0 ) );    // Initialize Lagrange Multiplier
    _moments = VectorVector( _nq, Vector( _nTotalEntries, 0.0 ) );      // Initialize Moments at quadrature points
    ComputeMoments();                                                   // Pre-Compute Moments


    _settings->SetNQuadPoints( _nqF );
    _Q = _problem->GetExternalSource( _energies );
    _solF = _problem->SetupIC();
    _settings->SetNQuadPoints( _nq );
    _sol = VectorVector( _nCells, Vector( _nq, 0.0 )); // -> Initial condition is completely in uncollided part
    _solNew = _sol;
    _solNewF = _solF; // I am just here to get the right size!
    _scatteringKernelFC = k->GetScatteringKernelFirstCollision( _nqF, _weightsF );       // Scattering Kernel from Fine to Coarse mesh
    delete k;
    _QFirstCollision = VectorVector(_nCells, Vector( _nq ));

}

void FirstCollisionSolver::Solve(){

    // Pre Solver Processing

    PrepareVolumeOutput();
    DrawPreSolverOutput();

    for ( unsigned idx_energy = 0; idx_energy < _nEnergies; idx_energy++){

        IterPreprocessing( idx_energy );
        std::cout << "Check Preprocessing" << std::endl;
        std::cout << "Sol size:  " << _solF.size() << " , " << _solF[0].size() << std::endl;
        std::cout << "Moments size:  " << _moments.size() << " , " << _moments[0].size() << std::endl;

        FluxUpdateUncollided();
        std::cout << "Check Uncollided Flux" << std::endl;

        FVMUpdateUncollided( idx_energy );
        std::cout << "Check Uncollided FVM" << std::endl;

        ComputeFirstCollisionSource( idx_energy );
        std::cout << "Check FC Source" << std::endl;

        FluxUpdateCollided( );
        std::cout << "Check Collided Flux" << std::endl;

        FVMUpdateCollided( idx_energy );
        std::cout << "Check Collided FVM" << std::endl;

        IterPostprocessing( );
        std::cout << "Check Postprocessing" << std::endl;


        WriteVolumeOutput( idx_energy );
        WriteScalarOutput( idx_energy );
        PrintScreenOutput( idx_energy );
        PrintHistoryOutput( idx_energy );
        PrintVolumeOutput( idx_energy );

        std::cout << "Check Prints" << std::endl;
    }

    // Post Solver Processing
    DrawPostSolverOutput();

}


////// --------- Solver Functions -------- /////

void FirstCollisionSolver::IterPreprocessing( unsigned idx_energy ){

    switch ( _settings->GetFirstCollisionSolver() ) {

        case SN_SOLVER:
            break; // Nothing to do here
        case MN_SOLVER:
            // ------- Reconstruction Step ----------------
            _optimizer->SolveMultiCell( _alpha, _solF, _moments );

            // ------- Relizablity Preservation Step ----
            for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
                ComputeRealizableSolution( idx_cell );    // already parallelized
            }
            break;
        default: // SN_SOLVER
            if ( idx_energy == 0)
                ErrorMessages::Error( "Chosen First Collision Solver not supported yet -> go on with SN SOLVER", CURRENT_FUNCTION );
            break;
    }
}

void FirstCollisionSolver::FluxUpdateUncollided(){

    if( _settings->GetFirstCollisionSolver() == MN_SOLVER ){

        _psiDx = VectorVector( _nCells, Vector( _nTotalEntries, 0.0 ) );
        _psiDy = VectorVector( _nCells, Vector( _nTotalEntries, 0.0) );
        if( _reconsOrder > 1 ) {
            _mesh->ReconstructSlopesU( _nTotalEntries, _psiDx, _psiDy, _alpha );    // unstructured reconstruction
        }
        std::cout << " check reconstruction" << std::endl;
        std::cout << "psidx size:  " << _psiDx.size() << " , " << _psiDx[0].size() << std::endl;
        std::cout << "psidy size:  " << _psiDy.size() << " , " << _psiDy[0].size() << std::endl;
        std::cout << "alpha size:  " << _alpha.size() << " , " << _alpha[0].size() << std::endl;
        //#pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
            // Dirichlet Boundaries stayd
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
            _solNewF[idx_cell] = ConstructFluxMN( idx_cell );
            std::cout << "cell" << idx_cell;
        }
        std::cout << std::endl;
    }
    else{ // Default is SN_SOLVER

        _psiDx = VectorVector( _nCells, Vector( _nqF, 0.0 ) );
        _psiDy = VectorVector( _nCells, Vector( _nqF, 0.0 ) );
        if( _reconsOrder > 1 ) {
            _mesh->ReconstructSlopesU( _nqF, _psiDx, _psiDy, _solF );    // unstructured reconstruction
        }
#pragma omp parallel for
        for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
            // Dirichlet Boundaries stayd
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
            _solNewF[idx_cell] = ConstructFluxSN( idx_cell, _nqF );
        }
    }
}


void FirstCollisionSolver::FVMUpdateUncollided( unsigned idx_energy ){

    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;

        // Flux Update and Absorption Term
        if ( _settings->GetFirstCollisionSolver() == MN_SOLVER) {
            // loop over moment entries
            for( unsigned idx_system = 0; idx_system < _nTotalEntries; idx_system++ ) {
                _solNewF[idx_cell][idx_system] = _solF[idx_cell][idx_system] - ( _dE / _areas[idx_cell] ) * _solNewF[idx_cell][idx_system]
                                                 - _dE  * _sigmaT[idx_energy][idx_cell] * _solF[idx_cell][idx_system] ;
            }
        }
        else { // default SN
            // loop over all ordinates
            for( unsigned idx_quad = 0; idx_quad < _nqF; ++idx_quad ) {
                // time update angular flux with numerical flux and total scattering cross section
                _solNewF[idx_cell][idx_quad] = _solF[idx_cell][idx_quad] - ( _dE / _areas[idx_cell] ) * _solNewF[idx_cell][idx_quad]
                                               - _dE * _sigmaT[idx_energy][idx_cell] * _solF[idx_cell][idx_quad];
            }
        }

        // Source Term
        if( _Q.size() == 1u ) {                   // constant source for all energies
            if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                _solNewF[idx_cell] += _dE * _Q[0][idx_cell][0];
            else
                _solNewF[idx_cell] += _dE * _Q[0][idx_cell];
        }
        else {
            if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                _solNewF[idx_cell] += _dE * _Q[idx_energy][idx_cell][0];
            else
                _solNewF[idx_cell] += _dE * _Q[idx_energy][idx_cell];
        }
    }
}


void FirstCollisionSolver::ComputeFirstCollisionSource( unsigned idx_energy ){

    // Q = _sigmaS * int_S^2 ( Sigma * psi_u ) d Omega'
    // isotropic scattering -> Sigma independent of Omega'
    std::cout << "FC Kernel size: (" << _scatteringKernelFC.rows() << "," << _scatteringKernelFC.columns() << ")" << std::endl;
    std::cout << "solNewF size: (" << _solNewF[0].size() << ")" << std::endl;
    std::cout << "QFirstCollision size: (" << _QFirstCollision.size() << " , " << _QFirstCollision[0].size() << ")" << std::endl;
    switch ( _settings->GetFirstCollisionSolver() ) {
        case SN_SOLVER: // -> size( Q ) = _nCells * _nqC;
            for (unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
                _QFirstCollision[ idx_cell ] = _sigmaS[ idx_energy ][ idx_cell ] * _scatteringKernelFC * _solNewF[idx_cell]; // sol or solNew?????
            }
            break;
        case MN_SOLVER: // size( Q ) = _nCells * _nqC;
            for ( unsigned idx_cell = 0; idx_cell <_nCells; idx_cell ++){
                for ( unsigned idx_quad = 0; idx_quad < _nq; idx_quad ++)
                    _QFirstCollision[ idx_cell ][ idx_quad ] = _sigmaS[idx_energy][ idx_cell ] * _scatterMatDiag[ 0 ] * _solNewF[ idx_cell ][ 0 ];
            }
            break;
        default: // SN_SOLVER
            for (unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
                _QFirstCollision[ idx_cell ] = _sigmaS[ idx_energy ][ idx_cell ] * _scatteringKernelFC * _solNewF[idx_cell];
            }
            break;
    }
}


void FirstCollisionSolver::FluxUpdateCollided( ){

    _psiDx = VectorVector( _nCells, Vector( _nq, 0.0 ) );
    _psiDy = VectorVector( _nCells, Vector( _nq, 0.0 ) );
    if( _reconsOrder > 1 ) {
        _mesh->ReconstructSlopesU( _nq, _psiDx, _psiDy, _sol );    // unstructured reconstruction
    }
#pragma omp parallel for
    for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++){
        // Dirichlet Boundaries stayd
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        _solNew[idx_cell] = ConstructFluxSN( idx_cell, _nq );
    }
}

void FirstCollisionSolver::FVMUpdateCollided( unsigned idx_energy ){

    std::cout << "scattering kernel size: " << _scatteringKernel.rows() << " , " << _scatteringKernel.columns() << std::endl;
    std::cout << "sol size: " << _sol.size() << " , " << _sol[0].size() << std::endl;
    std::cout << "solNew size: " << _solNew.size() << " , " << _solNew[0].size() << std::endl;

    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // loop over all ordinates
        for( unsigned idx_quad = 0; idx_quad < _nq; ++idx_quad ) {
            // time update angular flux with numerical flux and total scattering cross section
            _solNew[idx_cell][idx_quad] = _sol[idx_cell][idx_quad] - ( _dE / _areas[idx_cell] ) * _solNew[idx_cell][idx_quad] -
                                          _dE * _sigmaT[idx_energy][idx_cell] * _sol[idx_cell][idx_quad];
        }
        // compute scattering effects
        _solNew[idx_cell] += _dE * _sigmaS[idx_energy][idx_cell] * _scatteringKernel * _sol[idx_cell];    // multiply scattering matrix with psi

        // First Collision Source
        _solNew[ idx_cell ] += _dE * _QFirstCollision[ idx_cell ];

    }
}

void FirstCollisionSolver::IterPostprocessing( ){

    _sol = _solNew;
    for ( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell ++)
        _solF[ idx_cell ] = _solNewF[ idx_cell ]; // maybe subtract first collision source because these are already scattered

    ComputeRadFlux();
}


///// -------- Flux Constructors ------ /////

Vector FirstCollisionSolver::ConstructFluxSN( unsigned idx_cell, unsigned nq ){

    VectorVector sol;
    VectorVector quadPoints;
    if ( nq == _nq ){
        quadPoints = _quadPoints;
        sol = _sol;
    }
    else if ( nq == _nqF ){
        quadPoints = _quadPointsF;
        sol = _solF;
    }
    Vector flux( nq );
    double psiL;
    double psiR;

    for( unsigned idx_quad = 0; idx_quad < nq ; ++idx_quad ) {
        // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
        for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); ++idx_neighbor ) {
            // store flux contribution on psiNew_sigmaS to save memory
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                flux[idx_quad] +=
                        _g->Flux( quadPoints[idx_quad], sol[idx_cell][idx_quad], sol[idx_cell][idx_quad], _normals[idx_cell][idx_neighbor] );
            else {
                switch( _reconsOrder ) {
                    // first order solver
                    case 1:
                        flux[idx_quad] += _g->Flux( quadPoints[idx_quad],
                                                    sol[idx_cell][idx_quad],
                                                    sol[_neighbors[idx_cell][idx_neighbor]][idx_quad],
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
                        flux[idx_quad] += _g->Flux( quadPoints[idx_quad], psiL, psiR, _normals[idx_cell][idx_neighbor] );
                        break;
                        // higher order solver
                    case 3:
                        std::cout << "higher order is WIP" << std::endl;
                        break;
                        // default: first order solver
                    default:
                        flux[idx_quad] += _g->Flux( quadPoints[idx_quad],
                                                    sol[idx_cell][idx_quad],
                                                    sol[_neighbors[idx_cell][idx_neighbor]][idx_quad],
                                                    _normals[idx_cell][idx_neighbor] );
                }
            }
        }
    }
    return flux;
}

Vector FirstCollisionSolver::ConstructFluxMN( unsigned idx_cell ) {

    //--- Integration of moments of flux ---
    double entropyL, entropyR, entropyFlux;

    Vector flux( _nTotalEntries, 0.0 );

    //--- Temporary storages of reconstructed alpha ---
    Vector alphaL( _nTotalEntries, 0.0 );
    Vector alphaR( _nTotalEntries, 0.0 );

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        entropyFlux = 0.0;    // reset temporary flux

        for( unsigned idx_neigh = 0; idx_neigh < _neighbors[idx_cell].size(); idx_neigh++ ) {
            // Left side reconstruction
            if( _reconsOrder > 1 ) {
                alphaL = _alpha[idx_cell] + _psiDx[idx_cell] * ( _interfaceMidPoints[idx_cell][idx_neigh][0] - _cellMidPoints[idx_cell][0] ) +
                         _psiDy[idx_cell] * ( _interfaceMidPoints[idx_cell][idx_neigh][1] - _cellMidPoints[idx_cell][1] );
            }
            else {
                alphaL = _alpha[idx_cell];
            }
            entropyL = _entropy->EntropyPrimeDual( blaze::dot( alphaL, _moments[idx_quad] ) );

            // Right side reconstruction
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neigh] == _nCells )
                entropyR = entropyL;
            else {
                if( _reconsOrder > 1 ) {
                    alphaR = _alpha[_neighbors[idx_cell][idx_neigh]] +
                             _psiDx[_neighbors[idx_cell][idx_neigh]] *
                             ( _interfaceMidPoints[idx_cell][idx_neigh][0] - _cellMidPoints[_neighbors[idx_cell][idx_neigh]][0] ) +
                             _psiDy[_neighbors[idx_cell][idx_neigh]] *
                             ( _interfaceMidPoints[idx_cell][idx_neigh][1] - _cellMidPoints[_neighbors[idx_cell][idx_neigh]][1] );
                }
                else {
                    alphaR = _alpha[_neighbors[idx_cell][idx_neigh]];
                }
                entropyR = _entropy->EntropyPrimeDual( blaze::dot( alphaR, _moments[idx_quad] ) );
            }

            // Entropy flux
            entropyFlux += _g->Flux( _quadPoints[idx_quad], entropyL, entropyR, _normals[idx_cell][idx_neigh] );
        }
        // Solution flux
        flux += _moments[idx_quad] * ( _weights[idx_quad] * entropyFlux );
    }
    return flux;
}

//// ------------ Helper Functions ---------- ////

void FirstCollisionSolver::ComputeMoments(){
    double my, phi;

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        my  = _quadPointsSphere[idx_quad][0];
        phi = _quadPointsSphere[idx_quad][1];

        _moments[idx_quad] = _basis->ComputeSphericalBasis( my, phi );
    }
}

void FirstCollisionSolver::ComputeRealizableSolution( unsigned idx_cell ){

    double entropyReconstruction = 0.0;
    _solF[idx_cell]               = 0;
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        // Make entropyReconstruction a member vector, s.t. it does not have to be re-evaluated in ConstructFlux
        entropyReconstruction = _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_cell], _moments[idx_quad] ) );
        _solF[idx_cell] += _moments[idx_quad] * ( _weights[idx_quad] * entropyReconstruction );
    }
}

void FirstCollisionSolver::ComputeRadFlux(){

    FILE * radFluxfile = fopen((_settings->GetLogDir() + "/radFluxFC.txt").c_str(), "a+");
    fprintf(radFluxfile, "\n Iteration: \n ");
    fprintf(radFluxfile, "\t uncoll \t coll \t \t cumulative \n \n ");
    double f1 = 0, f2 = 0, f3 = 0;
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        // flux from uncollided
        double fluxUC;
        if ( _settings->GetFirstCollisionSolver() == MN_SOLVER )
            fluxUC = _solF[idx_cell][0];
        else
            fluxUC = blaze::dot( _solF[idx_cell], _weightsF);
        // flux from collided
        double fluxC = blaze::dot( _sol[idx_cell], _weights );
        // cumulative flux
        _fluxNew[idx_cell] = fluxUC + fluxC;

        f1 += fluxUC;
        f2 += fluxC;
        f3 += _fluxNew[idx_cell];
        //_fluxNew[idx_cell] = blaze::dot( _sol[idx_cell], _weights );
        //fprintf(radFluxfile, "\t %f \n", (float) _fluxNew[idx_cell] * 1000);
    }

    fprintf(radFluxfile, "\t %f \t %f \t %f \n", (float) f1 * 1000, (float) f2 * 1000, (float) f3 * 1000);
    fclose(radFluxfile);
}

void FirstCollisionSolver::FluxUpdate(){}             // not used
void FirstCollisionSolver::FVMUpdate( unsigned ){}    // not used

// -------- Output generators ( copied from SN_Solver ) -------------

void FirstCollisionSolver::PrepareVolumeOutput() {
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


void FirstCollisionSolver::WriteVolumeOutput( unsigned idx_pseudoTime ) {
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();

    if( ( _settings->GetVolumeOutputFrequency() != 0 && idx_pseudoTime % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
        ( idx_pseudoTime == _nEnergies - 1 ) ) { // need sol at last iteration

        for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
            switch( _settings->GetVolumeOutput()[idx_group] ) {
                case MINIMAL:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = _fluxNew[idx_cell];
                    }
                    break;

                case ANALYTIC:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        double time                           = idx_pseudoTime * _dE;
                        _outputFields[idx_group][0][idx_cell] = _problem->GetAnalyticalSolution(
                                _mesh->GetCellMidPoints()[idx_cell][0], _mesh->GetCellMidPoints()[idx_cell][1], time, _sigmaS[idx_pseudoTime][idx_cell] );
                    }
                    break;

                default: ErrorMessages::Error( "Volume Output Group not defined for First Collision SN Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
}