#include "solvers/csdsolvertrafofp2d.h"
#include "common/config.h"
#include "common/io.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/icru.h"
#include "problems/problembase.h"
#include "quadratures/qproduct.h"
#include "quadratures/quadraturebase.h"
#include "solvers/csdpn_starmap_constants.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

CSDSolverTrafoFP2D::CSDSolverTrafoFP2D( Config* settings ) : SNSolver( settings ) {
    std::cout << "FP" << std::endl;
    _dose = std::vector<double>( _settings->GetNCells(), 0.0 );

    // Set angle and energies
    _energies  = Vector( _nEnergies, 0.0 );    // equidistant
    _energyMin = 1e-4 * 0.511;
    _energyMax = _settings->GetMaxEnergyCSD();    // 5e0;

    // determine transformed energy grid for tabulated grid
    Vector E_transformed( E_trans.size(), 0.0 );
    for( unsigned i = 1; i < E_trans.size(); ++i )
        E_transformed[i] = E_transformed[i - 1] + ( E_tab[i] - E_tab[i - 1] ) / 2 * ( 1.0 / S_tab[i] + 1.0 / S_tab[i - 1] );

    // determine minimal and maximal energies
    double minE = 5e-5;
    double maxE = _settings->GetMaxEnergyCSD();
    _E_cutoff   = maxE;

    // define interpolation from energies to corresponding transformed energies \tilde{E} (without substraction of eMaxTrafo)
    Interpolation interpEToTrafo( E_tab, E_transformed );
    double eMaxTrafo = interpEToTrafo( maxE );
    double eMinTrafo = interpEToTrafo( minE );

    // Define intepolation back to original energy
    Interpolation interpTrafoToE( E_transformed, E_tab );

    // define linear grid in fully transformed energy \tilde\tilde E (cf. Dissertation Kerstion Kuepper, Eq. 1.25)
    _eTrafo        = blaze::linspace( _nEnergies, eMaxTrafo - eMaxTrafo, eMaxTrafo - eMinTrafo );
    double dETrafo = _eTrafo[1] - _eTrafo[0];

    // compute Corresponding original energies
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        _energies[n] = interpTrafoToE( eMaxTrafo - _eTrafo[n] );
    }

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

    unsigned jMinus, iMinus, jPlus, iPlus;

    // implementation of 2D spherical Laplacian according book "Advances in Discrete Ordinates Methodology", equation (1.136)
    for( unsigned j = 0; j < orderMu; ++j ) {
        DMinus = DPlus;
        DPlus  = DMinus - 2 * _mu[j] * _wp[j];
        for( unsigned i = 0; i < 2 * order; ++i ) {

            if( j > 0 ) {
                jMinus = j - 1;
            }
            else {
                jMinus = orderMu - 1;
            }
            if( i > 0 ) {
                iMinus = i - 1;
            }
            else {
                iMinus = 2 * order - 1;
            }
            if( j < orderMu - 1 ) {
                jPlus = j + 1;
            }
            else {
                jPlus = 0;
            }
            if( i < 2 * order - 1 ) {
                iPlus = i + 1;
            }
            else {
                iPlus = 0;
            }
            _L( j * 2 * order + i, jMinus * 2 * order + i ) = DMinus / ( _mu[j] - _mu[jMinus] ) / _wp[j];
            _L( j * 2 * order + i, j * 2 * order + i )      = -DMinus / ( _mu[j] - _mu[jMinus] ) / _wp[j];
            _L( j * 2 * order + i, j * 2 * order + iMinus ) = 1.0 / ( 1 - _mu[j] * _mu[j] ) * gamma[j] / ( _phi[i] - _phi[iMinus] ) / _wa[i];
            _L( j * 2 * order + i, j * 2 * order + i ) += -1.0 / ( 1 - _mu[j] * _mu[j] ) * gamma[j] / ( _phi[i] - _phi[iMinus] ) / _wa[i];
            _L( j * 2 * order + i, jPlus * 2 * order + i ) = DPlus / ( _mu[jPlus] - _mu[j] ) / _wp[j];
            _L( j * 2 * order + i, j * 2 * order + i ) += -DPlus / ( _mu[jPlus] - _mu[j] ) / _wp[j];
            _L( j * 2 * order + i, j * 2 * order + iPlus ) = 1.0 / ( 1 - _mu[j] * _mu[j] ) * gamma[j] / ( _phi[iPlus] - _phi[i] ) / _wa[i];
            _L( j * 2 * order + i, j * 2 * order + i ) += -1.0 / ( 1 - _mu[j] * _mu[j] ) * gamma[j] / ( _phi[iPlus] - _phi[i] ) / _wa[i];
        }
    }

    // Heney-Greenstein parameter
    double g = 0.8;

    // determine moments of Heney-Greenstein
    _xi = Matrix( 4, _nEnergies );

    // initialize stopping power vector
    _s = Vector( _nEnergies, 1.0 );

    _RT = true;

    // read in medical data if radiation therapy option selected
    ICRU database( _mu, _energies, _settings );
    database.GetTransportCoefficients( _xi );
    database.GetStoppingPower( _s );

    _density = std::vector<double>( _nCells, 1.0 );
    _sol     = _problem->SetupIC();
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
                        _outputFields[idx_group][0][idx_cell] += _dose[idx_cell];
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
void CSDSolverTrafoFP2D::FVMUpdate( unsigned /*idx_energy*/ ) {
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
    _dE        = _eTrafo[idx_pseudotime + 1] - _eTrafo[idx_pseudotime];

    double xi1 = _xi( 1, n );
    double xi2 = _xi( 2, n );
    double xi3 = _xi( 3, n );

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

// add FP scattering term implicitly
#pragma omp parallel for
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
        _sol[j] = _IL * blaze::solve( _IL - _dE * _alpha * _L, _sol[j] );
    }
}

void CSDSolverTrafoFP2D::IterPostprocessing( unsigned idx_pseudotime ) {
    unsigned n = idx_pseudotime;
    // --- Update Solution ---
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
        _sol[j] = _solNew[j];
    }

    for( unsigned j = 0; j < _nCells; ++j ) {
        _fluxNew[j] = dot( _sol[j], _weights );
        if( n > 0 ) {
            _dose[j] += 0.5 * _dE * ( _fluxNew[j] * _s[n] + _flux[j] * _s[n] ) / _density[j];    // update dose with trapezoidal rule
        }
        else {
            _dose[j] += _dE * _fluxNew[j] * _s[n] / _density[j];
        }
        _flux[j] = _fluxNew[j];
    }

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();
}

void CSDSolverTrafoFP2D::SolverPreprocessing() {
    auto log = spdlog::get( "event" );

    _densityMin = 0.1;
    for( unsigned j = 0; j < _nCells; ++j ) {
        _density[j] = 1.0;
        if( _density[j] < _densityMin ) _density[j] = _densityMin;
    }

    // save original energy field for boundary conditions
    _energiesOrig = _energies;

    // setup identity matrix for FP scattering
    _identity = Matrix( _nq, _nq, 0.0 );

    for( unsigned k = 0; k < _nq; ++k ) _identity( k, k ) = 1.0;

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

// do substitution from psi to psiTildeHat (cf. Dissertation Kerstion Kuepper, Eq. 1.23)
#pragma omp parallel for
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nq; ++k ) {
            _sol[j][k] = _sol[j][k] * _density[j] * _s[0];    // note that _s[_nEnergies - 1] is stopping power at highest energy
        }
    }

    // determine minimal density for CFL computation
    _densityMin = _density[0];
    for( unsigned j = 1; j < _nCells; ++j ) {
        if( _densityMin > _density[j] ) _densityMin = _density[j];
    }
    // cross sections do not need to be transformed to ETilde energy grid since e.g. TildeSigmaT(ETilde) = SigmaT(E(ETilde))
}
