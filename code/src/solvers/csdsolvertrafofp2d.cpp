#include "solvers/csdsolvertrafofp2d.h"
#include "common/config.h"
#include "common/io.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/icru.h"
#include "problems/problembase.h"
#include "quadratures/qproduct.h"
#include "quadratures/quadraturebase.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

CSDSolverTrafoFP2D::CSDSolverTrafoFP2D( Config* settings ) : SNSolver( settings ) {
    _dose = std::vector<double>( _settings->GetNCells(), 0.0 );

    // Set angle and energies
    _energies  = Vector( _nEnergies, 0.0 );    // equidistant
    _energyMin = 1e-4 * 0.511;
    _energyMax = 0.01;    // 50e0;

    // write equidistant energy grid (false) or refined grid (true)
    GenerateEnergyGrid( false );

    // create quadrature
    unsigned order    = _quadrature->GetOrder();
    unsigned nq       = _settings->GetNQuadPoints();
    _quadPoints       = _quadrature->GetPoints();
    _weights          = _quadrature->GetWeights();
    _quadPointsSphere = _quadrature->GetPointsSphere();

    unsigned orderMu;
    if( _settings->GetQuadName() == QUAD_GaussLegendreTensorized ) orderMu = order;
    if( _settings->GetQuadName() == QUAD_Product ) orderMu = 2 * order;

    // transform structured quadrature
    _mu  = Vector( orderMu );
    _phi = Vector( 2 * order );
    _wp  = Vector( orderMu );
    _wa  = Vector( 2 * order );

    // create quadrature 1D to compute mu grid
    QuadratureBase* quad1D = QuadratureBase::Create( QUAD_GaussLegendre1D, orderMu );
    Vector w               = quad1D->GetWeights();
    VectorVector muVec     = quad1D->GetPoints();

    for( unsigned k = 0; k < orderMu; ++k ) {
        _mu[k] = muVec[k][0];
        _wp[k] = w[k];
    }

    for( unsigned i = 0; i < 2 * order; ++i ) {
        _phi[i] = ( i + 0.5 ) * M_PI / order;
        _wa[i]  = M_PI / order;
    }

    // setup Laplace Beltrami matrix L in slab geometry
    _L = Matrix( nq, nq, 0.0 );

    _FPMethod = 2;

    double DMinus = 0.0;    // is beta_{n-1/2}
    double DPlus  = 0.0;    // is beta_{n+1/2}

    Vector gamma( 2 * order, 0.0 );
    double dPlus;
    double c, K;

    // Setup Coefficients (see "Advances in Discrete Ordinates Methodology",  eq (1.137-1.140))
    double dMinus = 0.0;
    DPlus         = DMinus - 2 * _mu[0] * w[0];
    for( unsigned j = 0; j < orderMu - 1; ++j ) {
        DMinus   = DPlus;
        DPlus    = DMinus - 2 * _mu[j] * w[j];
        dPlus    = ( sqrt( 1 - _mu[j + 1] * _mu[j + 1] ) - sqrt( 1 - _mu[j] * _mu[j] ) ) / ( _mu[j + 1] - _mu[j] );
        c        = ( DPlus * dPlus - DMinus * dMinus ) / _wp[j];
        K        = 2 * ( 1 - _mu[j] * _mu[j] ) + c * sqrt( 1 - _mu[j] * _mu[j] );
        gamma[j] = M_PI * M_PI * K / ( 2 * order * ( 1 - std::cos( M_PI / order ) ) );
        dMinus   = dPlus;
    }

    DPlus = 0.0;

    // implementation of 2D spherical Laplacian according book "Advances in Discrete Ordinates Methodology", equation (1.136)
    for( unsigned j = 0; j < orderMu; ++j ) {
        DMinus = DPlus;
        DPlus  = DMinus - 2 * _mu[j] * _wp[j];
        for( unsigned i = 0; i < 2 * order; ++i ) {
            if( j > 0 ) {
                _L( j * 2 * order + i, ( j - 1 ) * 2 * order + i ) = DMinus / ( _mu[j] - _mu[j - 1] ) / _wp[j];
                _L( j * 2 * order + i, j * 2 * order + i )         = -DMinus / ( _mu[j] - _mu[j - 1] ) / _wp[j];
            }
            if( i > 0 ) {
                _L( j * 2 * order + i, j * 2 * order + i - 1 ) = 1.0 / ( 1 - _mu[j] * _mu[j] ) * gamma[j] / ( _phi[i] - _phi[i - 1] ) / _wa[i];
                _L( j * 2 * order + i, j * 2 * order + i ) += -1.0 / ( 1 - _mu[j] * _mu[j] ) * gamma[j] / ( _phi[i] - _phi[i - 1] ) / _wa[i];
            }
            if( j < orderMu - 1 ) {
                _L( j * 2 * order + i, ( j + 1 ) * 2 * order + i ) = DPlus / ( _mu[j + 1] - _mu[j] ) / _wp[j];
                _L( j * 2 * order + i, j * 2 * order + i ) += -DPlus / ( _mu[j + 1] - _mu[j] ) / _wp[j];
            }
            if( i < 2 * order - 1 ) {
                _L( j * 2 * order + i, j * 2 * order + i + 1 ) = 1.0 / ( 1 - _mu[j] * _mu[j] ) * gamma[j] / ( _phi[i + 1] - _phi[i] ) / _wa[i];
                _L( j * 2 * order + i, j * 2 * order + i ) += -1.0 / ( 1 - _mu[j] * _mu[j] ) * gamma[j] / ( _phi[i + 1] - _phi[i] ) / _wa[i];
            }
        }
    }

    // Heney-Greenstein parameter
    double g = 0.8;

    // determine momente of Heney-Greenstein
    _xi1 = Vector( _nEnergies, 1.0 - g );    // paper Olbrant, Frank (11)
    _xi2 = Vector( _nEnergies, 4.0 / 3.0 - 2.0 * g + 2.0 / 3.0 * g * g );
    _xi  = Matrix( 4, _nEnergies );
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        _xi( 1, n ) = 1.0 - g;
        _xi( 2, n ) = 4.0 / 3.0 - 2.0 * g + 2.0 / 3.0 * g * g;
    }

    // initialize stopping power vector
    _s = Vector( _nEnergies, 1.0 );

    _RT = true;

    // read in medical data if radiation therapy option selected
    if( _RT ) {
        /*
        _nEnergies = 100;
        _energies.resize(_nEnergies);
        _xi = Matrix(6,_nEnergies);
        double minExp = -4.0;
        double  maxExp = 1.5;
        for( int n = 0; n<_nEnergies; ++n){
            double exponent = minExp + ( maxExp - minExp ) / ( _nEnergies - 1 ) * n;
            _energies[n] = pow(10.0,exponent);
        }*/
        ICRU database( _mu, _energies, _settings );
        database.GetTransportCoefficients( _xi );
        database.GetStoppingPower( _s );
        /*
        // print coefficients
        std::cout<<"E = [";
        for( unsigned n = 0; n<_nEnergies; ++n){
            std::cout<<_energies[n]<<"; ";
        }
        std::cout<<"];"<<std::endl;
        std::cout<<"xi = [";
        for( unsigned n = 0; n<_nEnergies; ++n){
            std::cout<<_xi(0,n)<<" "<<_xi(1,n)<<" "<<_xi(2,n)<<" "<<_xi(3,n)<<" "<<_xi(4,n)<<" "<<_xi(5,n)<<"; ";
        }
        std::cout<<"];"<<std::endl;
        */
    }

    //_density = std::vector<double>( _nCells, 1.0 );
    // exit(EXIT_SUCCESS);
}

void CSDSolverTrafoFP2D::GenerateEnergyGrid( bool refinement ) {
    _dE = ComputeTimeStep( _settings->GetCFL() );
    if( !refinement ) {
        _nEnergies = unsigned( ( _energyMax - _energyMin ) / _dE );
        _energies.resize( _nEnergies );
        for( unsigned n = 0; n < _nEnergies; ++n ) {
            _energies[n] = _energyMin + ( _energyMax - _energyMin ) / ( _nEnergies - 1 ) * n;
        }
    }
    else {
        // hard-coded positions for energy-grid refinement
        double energySwitch    = 0.58;
        double energySwitchMin = 0.03;

        // number of energies per intervals [E_min,energySwitchMin], [energySwitchMin,energySwitch], [energySwitch,E_Max]
        unsigned nEnergies1 = unsigned( ( _energyMax - energySwitch ) / _dE );
        unsigned nEnergies2 = unsigned( ( energySwitch - energySwitchMin ) / ( _dE / 2 ) );
        unsigned nEnergies3 = unsigned( ( energySwitchMin - _energyMin ) / ( _dE / 3 ) );
        _nEnergies          = nEnergies1 + nEnergies2 + nEnergies3 - 2;
        _energies.resize( _nEnergies );

        // write equidistant energy grid in each interval
        for( unsigned n = 0; n < nEnergies3; ++n ) {
            _energies[n] = _energyMin + ( energySwitchMin - _energyMin ) / ( nEnergies3 - 1 ) * n;
        }
        for( unsigned n = 1; n < nEnergies2; ++n ) {
            _energies[n + nEnergies3 - 1] = energySwitchMin + ( energySwitch - energySwitchMin ) / ( nEnergies2 - 1 ) * n;
        }
        for( unsigned n = 1; n < nEnergies1; ++n ) {
            _energies[n + nEnergies3 + nEnergies2 - 2] = energySwitch + ( _energyMax - energySwitch ) / ( nEnergies1 - 1 ) * n;
        }
    }
}

// IO
void CSDSolverTrafoFP2D::PrepareVolumeOutput() {
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

            default: ErrorMessages::Error( "Volume Output Group not defined for CSD_SN_FP_TRAFO Solver!", CURRENT_FUNCTION ); break;
        }
    }
}

void CSDSolverTrafoFP2D::WriteVolumeOutput( unsigned idx_pseudoTime ) {
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();
    double maxDose;
    if( ( _settings->GetVolumeOutputFrequency() != 0 && idx_pseudoTime % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
        ( idx_pseudoTime == _maxIter - 1 ) /* need sol at last iteration */ ) {

        for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
            switch( _settings->GetVolumeOutput()[idx_group] ) {
                case MINIMAL:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = _fluxNew[idx_cell];
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
                    }
                    // Compute normalized dose
                    _outputFields[idx_group][1] = _outputFields[idx_group][0];

                    maxDose = *std::max_element( _outputFields[idx_group][0].begin(), _outputFields[idx_group][0].end() );

                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][1][idx_cell] /= maxDose;
                    }
                    break;

                default: ErrorMessages::Error( "Volume Output Group not defined for CSD_SN_FP_TRAFO Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
}

// Solver
void CSDSolverTrafoFP2D::FVMUpdate( unsigned idx_energy ) {
// loop over all spatial cells
#pragma omp parallel for
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // loop over all ordinates
        for( unsigned i = 0; i < _nq; ++i ) {
            // time update angular flux with numerical flux and total scattering cross section
            _solNew[j][i] = _sol[j][i] - _dE * _solNew[j][i];
        }
    }
}

void CSDSolverTrafoFP2D::FluxUpdate() {
// loop over all spatial cells
#pragma omp parallel for
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // loop over all ordinates
        for( unsigned i = 0; i < _nq; ++i ) {
            _solNew[j][i] = 0.0;
            // loop over all neighbor cells (edges) of cell j and compute numerical fluxes
            for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[j].size(); ++idx_neighbor ) {
                // store flux contribution on psiNew_sigmaS to save memory
                if( _boundaryCells[j] == BOUNDARY_TYPE::NEUMANN && _neighbors[j][idx_neighbor] == _nCells )
                    continue;    // adiabatic wall, add nothing
                else
                    _solNew[j][i] += _g->Flux( _quadPoints[i],
                                               _sol[j][i] / _density[j],
                                               _sol[_neighbors[j][idx_neighbor]][i] / _density[_neighbors[j][idx_neighbor]],
                                               _normals[j][idx_neighbor] ) /
                                     _areas[j];
            }
        }
    }
}

void CSDSolverTrafoFP2D::IterPreprocessing( unsigned idx_pseudotime ) {
    unsigned n = idx_pseudotime;
    _dE        = fabs( _energies[n + 1] - _energies[n] );    // is the sign correct here?

    double xi1 = _xi( 1, _nEnergies - n - 1 );
    double xi2 = _xi( 2, _nEnergies - n - 1 );
    double xi3 = _xi( 3, _nEnergies - n - 1 );

    // setup coefficients in FP step
    if( _FPMethod == 1 ) {
        _alpha  = 0.0;
        _alpha2 = xi1 / 2.0;
        _beta   = 0.0;
    }
    else if( _FPMethod == 2 ) {
        _alpha  = xi1 / 2.0 + xi2 / 8.0;
        _alpha2 = 0.0;
        _beta   = xi2 / 8.0 / xi1;
    }
    else if( _FPMethod == 3 ) {
        _alpha  = xi2 * ( 27.0 * xi2 * xi2 + 5.0 * xi3 * xi3 - 24.0 * xi2 * xi3 ) / ( 8.0 * xi3 * ( 3.0 * xi2 - 2.0 * xi3 ) );
        _beta   = xi3 / ( 6.0 * ( 3.0 * xi2 - 2.0 * xi3 ) );
        _alpha2 = xi1 / 2.0 - 9.0 / 8.0 * xi2 * xi2 / xi3 + 3.0 / 8.0 * xi2;
    }

    _IL = _identity - _beta * _L;

    // write BC for water phantom
    if( _RT && false ) {
        for( unsigned k = 0; k < _nq; ++k ) {
            if( _quadPoints[k][0] > 0 ) {
                _sol[0][k] = 1e5 * exp( -200.0 * pow( 1.0 - _quadPoints[k][0], 2 ) ) *
                             exp( -50.0 * pow( _energyMax - _energiesOrig[_nEnergies - n - 1], 2 ) ) * _density[0] * _s[_nEnergies - n - 1];
            }
        }
    }

// add FP scattering term implicitly
#pragma omp parallel for
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
        //_sol[j] = blaze::solve( _identity - _dE * _alpha2 * _L, psiNew[j] );
        _sol[j] = _IL * blaze::solve( _IL - _dE * _alpha * _L, _sol[j] );
    }
}

void CSDSolverTrafoFP2D::IterPostprocessing() {
    // --- Update Solution ---
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
        _sol[j] = _solNew[j];
    }

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();
}

void CSDSolverTrafoFP2D::SolverPreprocessing() {
    auto log = spdlog::get( "event" );

    _densityMin = 0.1;
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _density[j] < _densityMin ) _density[j] = _densityMin;
    }

    // save original energy field for boundary conditions
    _energiesOrig = _energies;

    // setup incoming BC on left
    //_sol = VectorVector( _density.size(), Vector( _settings->GetNQuadPoints(), 0.0 ) );    // hard coded IC, needs to be changed
    // for( unsigned k = 0; k < _nq; ++k ) {
    //    if( _quadPoints[k][0] > 0 && !_RT ) _sol[0][k] = 1e5 * exp( -10.0 * pow( 1.0 - _quadPoints[k][0], 2 ) );
    //}
    // hard coded boundary type for 1D testcases (otherwise cells will be NEUMANN)
    //_boundaryCells[0]           = BOUNDARY_TYPE::DIRICHLET;
    //_boundaryCells[_nCells - 1] = BOUNDARY_TYPE::DIRICHLET;

    // setup identity matrix for FP scattering
    _identity = Matrix( _nq, _nq, 0.0 );

    for( unsigned k = 0; k < _nq; ++k ) _identity( k, k ) = 1.0;

    // angular flux at next time step
    VectorVector psiNew( _nCells, Vector( _nq, 0.0 ) );
    double dFlux = 1e10;
    Vector fluxNew( _nCells, 0.0 );
    Vector fluxOld( _nCells, 0.0 );
    // for( unsigned j = 0; j < _nCells; ++j ) {
    //    fluxOld[j] = dot( _sol[j], _weights );
    //}

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    // if( rank == 0 ) log->info( "{:10}   {:10}", "E", "dFlux" );

// do substitution from psi to psiTildeHat (cf. Dissertation Kerstion Kuepper, Eq. 1.23)
#pragma omp parallel for
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nq; ++k ) {
            _sol[j][k] = _sol[j][k] * _density[j] * _s[_nEnergies - 1];    // note that _s[_nEnergies - 1] is stopping power at highest energy
        }
    }
    // store transformed energies ETilde instead of E in _energies vector (cf. Dissertation Kerstion Kuepper, Eq. 1.25)
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
    _densityMin = _density[0];
    for( unsigned j = 1; j < _nCells; ++j ) {
        if( _densityMin > _density[j] ) _densityMin = _density[j];
    }
    // cross sections do not need to be transformed to ETilde energy grid since e.g. TildeSigmaT(ETilde) = SigmaT(E(ETilde))
}
