#include "solvers/csdsnsolverfp2d.h"
#include "common/config.h"
#include "common/io.h"
#include "common/mesh.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/problembase.h"
#include "quadratures/quadraturebase.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

CSDSNSolverFP2D::CSDSNSolverFP2D( Config* settings ) : SNSolver( settings ) {
    // quadrature
    unsigned order = _quadrature->GetOrder();
    unsigned nq = _settings->GetNQuadPoints();
    _quadPoints = _quadrature->GetPoints();
    _weights    = _quadrature->GetWeights();
    _quadPointsSphere = _quadrature->GetPointsSphere();

    // transform structured quadrature
    _mu = Vector( order );
    _phi = Vector( 2 * order );
    _wp = Vector( order );
    _wa = Vector( 2 * order );
    for( unsigned j = 0; j < order; ++j ) {
        for( unsigned i = 0; i < 2 * order; ++i ) {
            _mu[j] = _quadPointsSphere[j * ( 2 * order ) + i][0];
            _phi[i] = _quadPointsSphere[j * ( 2 * order ) + i][1];

            _wp[j] = _weights[j * ( 2 * order ) + i] * order / M_PI;
            _wa[i] = M_PI / order;
        }
    }

    // dose
    _dose = std::vector<double>( _settings->GetNCells(), 0.0 );

    // set angle and energies
    _angle           = Vector( nq, 0.0 );
    _energies        = Vector( _nEnergies, 0.0 );
    _energyMin = 1e-4;
    _energyMax = 5.0;

    // write equidistant energy grid
    _dE        = ComputeTimeStep( settings->GetCFL() );
    _nEnergies = unsigned( ( _energyMax - _energyMin ) / _dE );
    std::cout<<"nEnergies = "<<_nEnergies<<std::endl;
    _energies.resize( _nEnergies );
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        _energies[n] = _energyMin + ( _energyMax - _energyMin ) / ( _nEnergies - 1 ) * n;
    }

    // set polar and azimuthal matrix L & R
    _L = Matrix( order, order, 0.0 ); // polar: mu
    _R = Matrix( 2 * order, 2 * order, 0.0 ); // azimuthal: phi

    // polar matrix (left multiplication)
    double DMinus = 0.0;
    double DPlus  = 0.0;
    for( unsigned k = 0; k < order; ++k ) {
        DMinus = DPlus;
        DPlus  = DMinus - 2 * _mu[k] * _weights[k];

        if( k > 0 ) {
            _L( k, k - 1 ) = DMinus / ( _mu[k] - _mu[k - 1] ) / _wp[k];
            _L( k, k )     = -DMinus / ( _mu[k] - _mu[k - 1] ) / _wp[k];
        }
        if( k < order - 1 ) {
            _L( k, k + 1 ) = DPlus / ( _mu[k + 1] - _mu[k] ) / _wp[k];
            _L( k, k ) += -DPlus / ( _mu[k + 1] - _mu[k] ) / _wp[k];
        }
    }

    // azimuthal matrix (right multiplication)
    for( unsigned k = 0; k < 2 * order; ++k ) {
        if( k > 0 ) {
            _R( k - 1, k ) = 1.0 / ( _phi[k] - _phi[k - 1] ) / _wa[k];
            _R( k, k ) = -1.0 / ( _phi[k] - _phi[k - 1] ) / _wa[k];
        }
        if( k < 2 * order - 1 ) {
            _R( k + 1, k ) = 1.0 / ( _phi[k + 1] - _phi[k] ) / _wa[k];
            _R( k, k ) -= 1.0 / ( _phi[k + 1] - _phi[k] ) / _wa[k];
        }
    }

    // Heney-Greenstein parameter
    double g = 0.8;

    // determine momente of Heney-Greenstein
    _xi1 = Vector(_nEnergies, 1.0 - g); // paper Olbrant, Frank (11)
    _xi2 = Vector(_nEnergies, 4.0 / 3.0 - 2.0 * g + 2.0 / 3.0 * g * g);
    _xi = Matrix(3,_nEnergies);

    _s = Vector( _nEnergies, 1.0 );

    // read medical data if radiation therapy option selected
    _RT = true;
    if(_RT){
        //_nEnergies = 1000;
        //_energies.resize(_nEnergies);
        //_xi = Matrix(6,_nEnergies);
        //double minExp = -4.0;
        //double  maxExp = 1.5;
        //for( int n = 0; n<_nEnergies; ++n){
        //    double exponent = minExp + ( maxExp - minExp ) / ( _nEnergies - 1 ) * n;
        //    _energies[n] = pow(10.0,exponent);
        //}
        ICRU database( abs(_mu), _energies );
        database.GetTransportCoefficients( _xi );
        database.GetStoppingPower( _s );
        /*// print coefficients
        std::cout<<"E = [";
        for( unsigned n = 0; n<_nEnergies; ++n){
            std::cout<<_energies[n]<<"; ";
        }
        std::cout<<"];"<<std::endl;
        std::cout<<"xi = [";
        for( unsigned n = 0; n<_nEnergies; ++n){
            std::cout<<_xi(0,n)<<" "<<_xi(1,n)<<" "<<_xi(2,n)<<" "<<_xi(3,n)<<" "<<_xi(4,n)<<" "<<_xi(5,n)<<"; ";
        }
        std::cout<<"];"<<std::endl;*/
    }

    // recompute scattering kernel. TODO: add this to kernel function
    for( unsigned p = 0; p < _nq; ++p ) {
        for( unsigned q = 0; q < _nq; ++q ) {
            _scatteringKernel( p, q ) = 0.0;
        }
        _scatteringKernel( p, p ) = _weights[p];
    }

    _density = Vector( _nCells, 1.0 );
    //exit(EXIT_SUCCESS);
}

void CSDSNSolverFP2D::Solve() {
    std::cout << "Solve Fokker-Planck with Heney-Greenstein kernel using "<<_nEnergies<<" energies " << std::endl;
    auto log      = spdlog::get( "event" );
    auto cellMids = _mesh->GetCellMidPoints();

    // setup IC and incoming BC on left
    _sol = std::vector<Vector>( _nCells, Vector( _nq, 0.0 ) );
    for( unsigned k = 0; k < _nq; ++k ) {
        if( _quadPoints[k][0] > 0 && !_RT) _sol[0][k] = 1e5 * exp( -10.0 * pow( 1.0 - _quadPoints[k][0], 2 ) );
    }
    _boundaryCells[0]           = BOUNDARY_TYPE::DIRICHLET;
    _boundaryCells[_nCells - 1] = BOUNDARY_TYPE::DIRICHLET;

    unsigned order = _quadrature->GetOrder();
    //Matrix identity( _nq, _nq, 0.0 );
    //for( unsigned k = 0; k < _nq; ++k ) identity( k, k ) = 1.0;
    Matrix identity( order, order, 0.0 ); // unit matrix 
    for( unsigned k = 0; k < order; ++k ) identity( k, k ) = 1.0;

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew( _nCells, Vector( _nq, 0.0 ) );
    double dFlux = 1e10;
    Vector fluxNew( _nCells, 0.0 );
    Vector fluxOld( _nCells, 0.0 );

    auto psiStruct = std::vector<Matrix>( _nCells, Matrix( order, 2 * order, 0.0 ) );
    auto psi1 = psiStruct;
    auto psi2 = psiStruct;

    Vector coeff(order, 0.0);
    double A = 0.0;
    double A1 = 0.0;
    double A2 = 0.0;
    double A3 = 0.0;
    double A4Plus = 0.0;
    double A4Minus = 0.0;
    double DPlus = 0.0;
    double DMinus = 0.0;
    for( unsigned k = 0; k < order; ++k ) {
        DMinus = DPlus;
        DPlus  = DMinus - 2 * _mu[k] * _wp[k];

        if (k > 0 && k < order - 1) {
            A4Plus = (sqrt(1.0 - _mu[k+1] * _mu[k+1]) - sqrt(1.0 - _mu[k] * _mu[k])) / (_mu[k+1] - _mu[k]);
            A4Minus = (sqrt(1.0 - _mu[k] * _mu[k]) - sqrt(1.0 - _mu[k-1] * _mu[k-1])) / (_mu[k] - _mu[k-1]);
            A3 = (DPlus * A4Plus - DMinus * A4Minus) / _wp[k];
            A2 = 2.0 * (1.0 - _mu[k] * _mu[k]) + A3 * sqrt(1.0 - _mu[k] * _mu[k]);
            A1 = pow(M_PI, 2) * A2 / order / (1.0 - cos(2.0 * M_PI / _nq));
            A = A1 / (1.0 - _mu[k] * _mu[k]);
            coeff[k] = A;
        }
    }

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 ) log->info( "{:10}   {:10}", "E", "dFlux" );

    // revert energies and stopping power
    Vector sSave = _s;
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        _energies[n] = _energies[_nEnergies - 1] - _energies[n];
        _s[n]        = sSave[_nEnergies - 1 - n];
    }

    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nq; ++k ) {
            fluxOld[j] += _weights[k] * _sol[j][k] * _s[0];
        }
    }

    // determine minimal density for CFL computation
    double densityMin = _density[0];
    for( unsigned j = 1; j < _nCells; ++j ) {
        if( densityMin > _density[j] ) densityMin = _density[j];
    }

    //VectorVector psi1 = _sol;

    // cross sections do not need to be transformed to ETilde energy grid since e.g. TildeSigmaT(ETilde) = SigmaT(E(ETilde))
    
    // loop over energies (pseudo-time)
    for( unsigned n = 0; n < _nEnergies - 1; ++n ) {
        _dE = fabs( _energies[n + 1] - _energies[n] );    // is the sign correct here?

        // set alpha and beta
        _alpha = _xi(1,_nEnergies-n-1) / 2.0 + _xi(2,_nEnergies-n-1) / 8.0;
        _beta  = _xi(2,_nEnergies-n-1) / 8.0 / _xi(1,_nEnergies-n-1);

        _IL = identity - _beta * _L;

        // write BC
        if(_RT){
            for( unsigned k = 0; k < _nq; ++k ) {
                if( _quadPoints[k][0] > 0 ) {
                        _sol[0][k] = 1e5 * exp( -200.0 * pow( 1.0 - _quadPoints[k][0], 2 ) ) * exp(-50*pow(_energyMax-_energies[n],2));
                }
            }
        }

// loop over all spatial cells
// scattering step
//#pragma omp parallel for
        for( unsigned j = 0; j < _nCells; ++j ) {
            if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
            //psi1[j] = blaze::solve( _IL, _sol[j] );
            //psi1[j] = _alpha * _L * psi1[j];

            psi1[j] = blaze::solve( _L, psiStruct[j] );
        }
        std::cout<<"point 1"<<std::endl;
// advection step
        for( unsigned j = 0; j < _nCells; ++j ) {
            if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
            // loop over all ordinates
            for( unsigned i = 0; i < _nq; ++i ) {
                psiNew[j][i] = 0.0;
                // loop over all neighbor cells (edges) of cell j and compute numerical fluxes
                for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[j].size(); ++idx_neighbor ) {
                    // store flux contribution on psiNew_sigmaS to save memory
                    if( _boundaryCells[j] == BOUNDARY_TYPE::NEUMANN && _neighbors[j][idx_neighbor] == _nCells )
                        continue;    // adiabatic wall, add nothing
                    else
                        psiNew[j][i] += _g->Flux( _quadPoints[i], _sol[j][i], _sol[_neighbors[j][idx_neighbor]][i], _normals[j][idx_neighbor] ) /
                                        _areas[j] / ( _density[j] * _s[n + 1] );;
                }
                // time update angular flux with numerical flux and total scattering cross section
                //psiNew[j][i] = _sol[j][i] - _dE * psiNew[j][i] + _dE * psi1[j][i]/ ( _density[j] * _s[n + 1] ) +
                //        ( _s[n] / _s[n + 1] - 1 ) * _sol[j][i];
                psiNew[j][i] = _sol[j][i] - _dE * psiNew[j][i] + ( _s[n] / _s[n + 1] - 1 ) * _sol[j][i];

                // fokker-planck
                for( unsigned jj = 0; jj < order; ++jj ) {
                    for( unsigned ii = 0; ii < 2 * order; ++ii ) {
                        psiStruct[j](jj, ii) = _sol[j][jj * ( 2 * order ) + ii];
                    }
                }


            }
        }
        std::cout<<"point 2"<<std::endl;
        for( unsigned j = 0; j < _nCells; ++j ) {
            if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
            _sol[j] = psiNew[j];
        }
        std::cout<<"point 3"<<std::endl;
        for( unsigned j = 0; j < _nCells; ++j ) {
            fluxNew[j] = 0.0;
            for( unsigned k = 0; k < _nq; ++k ) {
                fluxNew[j] += _weights[k] * _sol[j][k] * _s[n + 1];
            }
            _dose[j] += 0.5 * _dE * ( fluxNew[j] + fluxOld[j] );    // update dose with trapezoidal rule
            _solverOutput[j] = fluxNew[j];
        }
        std::cout<<"point 4"<<std::endl;
        // Save( n );
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}  {:01.5e}  {:01.5e}", _energies[n], _dE / densityMin, dFlux );
        if( std::isinf( dFlux ) || std::isnan( dFlux ) ) break;
    }
    Save( 1 );
}

void CSDSNSolverFP2D::Save() const {
    std::vector<std::string> fieldNames{ "dose", "normalized dose" };
    std::vector<std::vector<double>> dose( 1, _dose );
    std::vector<std::vector<double>> normalizedDose( 1, _dose );
    double maxDose = *std::max_element( _dose.begin(), _dose.end() );
    for( unsigned i = 0; i < _dose.size(); ++i ) normalizedDose[0][i] /= maxDose;
    std::vector<std::vector<std::vector<double>>> results{ dose, normalizedDose };
    ExportVTK( _settings->GetOutputFile(), results, fieldNames, _mesh );
}

void CSDSNSolverFP2D::Save( int currEnergy ) const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<std::vector<double>> scalarField( 1, _solverOutput );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), results, fieldNames, _mesh );
}
