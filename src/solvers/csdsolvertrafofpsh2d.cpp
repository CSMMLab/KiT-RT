#include "solvers/csdsolvertrafofpsh2d.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "fluxes/numericalflux.hpp"
#include "kernels/scatteringkernelbase.hpp"
#include "problems/problembase.hpp"
#include "quadratures/qproduct.hpp"
#include "quadratures/quadraturebase.hpp"
#include "solvers/csdpn_starmap_constants.hpp"
#include "solvers/solverbase.hpp"
#include "toolboxes/icru.hpp"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

CSDSolverTrafoFPSH2D::CSDSolverTrafoFPSH2D( Config* settings ) : SNSolver( settings ) {
    _dose = std::vector<double>( _settings->GetNCells(), 0.0 );

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

    // check what happens if we intepolate back
    Interpolation interpTrafoToE( E_transformed, E_tab );

    // define linear grid in fully transformed energy \tilde\tilde E (cf. Dissertation Kerstion Kuepper, Eq. 1.25)
    _eTrafo = blaze::linspace( _nEnergies, eMaxTrafo - eMaxTrafo, eMaxTrafo - eMinTrafo );
    // double dETrafo = _eTrafo[1] - _eTrafo[0];

    // compute Corresponding original energies
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        _energies[n] = interpTrafoToE( eMaxTrafo - _eTrafo[n] );
    }

    // evaluate corresponding stopping powers and transport coefficients
    // compute sigmaT is now done during computation in IterPreprocessing()

    // compute stopping powers
    Vector etmp = E_tab;
    Vector stmp = S_tab;
    Interpolation interpS( etmp, stmp );
    _s = interpS( _energies );

    // compute stopping power between energies for dose computation
    double dE = _eTrafo[2] - _eTrafo[1];
    Vector eTrafoMid( _nEnergies - 1 );
    for( unsigned n = 0; n < _nEnergies - 1; ++n ) {
        eTrafoMid[n] = _eTrafo[n] + dE / 2;
    }
    // compute Corresponding original energies at intermediate points
    Vector eMid( _nEnergies - 1 );
    for( unsigned n = 0; n < _nEnergies - 1; ++n ) {
        eMid[n] = interpTrafoToE( eMaxTrafo - eTrafoMid[n] );
    }
    _sMid = interpS( eMid );

    _quadPointsSphere  = _quadrature->GetPointsSphere();    // (my,phi,r) gridpoints of the quadrature in spherical cordinates
    SphericalBase* sph = SphericalBase::Create( _settings );
    Matrix Y           = blaze::zero<double>( sph->GetBasisSize(), _nq );

    for( unsigned q = 0; q < _nq; q++ ) {
        blaze::column( Y, q ) = sph->ComputeSphericalBasis( _quadPointsSphere[q][0], _quadPointsSphere[q][1] );
    }

    _O = blaze::trans( Y );
    _M = _O;    // transposed to simplify weight multiplication. We will transpose _M later again

    unsigned nSph = _M.columns();
    for( unsigned j = 0; j < nSph; ++j ) {
        blaze::column( _M, j ) *= _weights;
    }

    _M.transpose();

    _S               = Matrix( nSph, nSph, 0.0 );
    unsigned counter = 0;
    for( int l = 0; l <= settings->GetMaxMomentDegree(); ++l ) {
        for( int m = -l; m <= l; ++m ) {
            _S( counter, counter ) = double( -l * ( l + 1 ) );
            counter++;
        }
    }
    _polyDegreeBasis = settings->GetMaxMomentDegree();

    _L = _O * _S * _M;

    // determine moments of Heney-Greenstein
    _xi = Matrix( 4, _nEnergies );
    ICRU database( _mu, _energies, _settings );
    database.GetTransportCoefficients( _xi );

    std::cout << "check I: " << _M * _O << std::endl;
    // exit( 0 );
    _RT = true;

    delete sph;
}

// IO
void CSDSolverTrafoFPSH2D::PrepareVolumeOutput() {
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

void CSDSolverTrafoFPSH2D::WriteVolumeOutput( unsigned idx_pseudoTime ) {
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
                        _outputFields[idx_group][0][idx_cell] = _dose[idx_cell];
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
    if( idx_pseudoTime == _nEnergies - 2 ) {
        std::ofstream out( _settings->GetOutputFile().append( ".txt" ) );
        unsigned nx = _settings->GetNCells();

        for( unsigned j = 0; j < nx; ++j ) {
            out << _cellMidPoints[j][0] << " " << _cellMidPoints[j][1] << " " << _dose[j] << std::endl;
        }
        out.close();
    }
}

// Solver
void CSDSolverTrafoFPSH2D::FVMUpdate( unsigned /*idx_energy*/ ) {
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

void CSDSolverTrafoFPSH2D::FluxUpdate() {
// loop over all spatial cells
#pragma omp parallel for
    for( unsigned j = 0; j < _nCells; ++j ) {
        double solL;
        double solR;
        if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // loop over all ordinates
        for( unsigned i = 0; i < _nq; ++i ) {
            _solNew[j][i] = 0.0;
            // loop over all neighbor cells (edges) of cell j and compute numerical fluxes
            for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[j].size(); ++idx_neighbor ) {
                // store flux contribution on psiNew_sigmaS to save memory
                if( _boundaryCells[j] == BOUNDARY_TYPE::NEUMANN && _neighbors[j][idx_neighbor] == _nCells )
                    continue;    // adiabatic wall, add nothing
                else {
                    unsigned int nbr_glob = _neighbors[j][idx_neighbor];    // global idx of neighbor cell

                    switch( _reconsOrder ) {
                        case 1: _solNew[j][i] += _g->Flux( _quadPoints[i], _sol[j][i], _sol[nbr_glob][i], _normals[j][idx_neighbor] ); break;
                        // second order solver
                        case 2:
                            // left status of interface
                            solL = _sol[j][i] / _density[j] +
                                   _limiter[j][i] * ( _solDx[j][i] * ( _interfaceMidPoints[j][idx_neighbor][0] - _cellMidPoints[j][0] ) +
                                                      _solDy[j][i] * ( _interfaceMidPoints[j][idx_neighbor][1] - _cellMidPoints[j][1] ) );
                            solR = _sol[nbr_glob][i] / _density[nbr_glob] +
                                   _limiter[nbr_glob][i] *
                                       ( _solDx[nbr_glob][i] * ( _interfaceMidPoints[j][idx_neighbor][0] - _cellMidPoints[nbr_glob][0] ) +
                                         _solDy[nbr_glob][i] * ( _interfaceMidPoints[j][idx_neighbor][1] - _cellMidPoints[nbr_glob][1] ) );

                            // flux evaluation
                            _solNew[j][i] += _g->Flux( _quadPoints[i], solL, solR, _normals[j][idx_neighbor] ) / _areas[j];
                            break;
                            // higher order solver
                        default: ErrorMessages::Error( "Reconstruction order not supported.", CURRENT_FUNCTION ); break;
                    }
                }
            }
        }
    }
}

void CSDSolverTrafoFPSH2D::IterPreprocessing( unsigned idx_pseudotime ) {
    if( _reconsOrder > 1 ) {
        VectorVector solDivRho = _sol;
        for( unsigned j = 0; j < _nCells; ++j ) {
            solDivRho[j] = _sol[j] / _density[j];
        }
        _mesh->ComputeSlopes( _nq, _solDx, _solDy, solDivRho );
        _mesh->ComputeLimiter( _nq, _solDx, _solDy, solDivRho, _limiter );
    }
    bool old = true;
    if( old ) {
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
    else {
        _dE = _eTrafo[idx_pseudotime + 1] - _eTrafo[idx_pseudotime];
        Vector sigmaSAtEnergy( _polyDegreeBasis + 1, 0.0 );
        Vector sigmaTAtEnergy( _polyDegreeBasis + 1, 0.0 );
        // compute scattering cross section at current energy
        for( unsigned idx_degree = 0; idx_degree <= _polyDegreeBasis; ++idx_degree ) {
            // setup interpolation from E to sigma at degree idx_degree
            Interpolation interp( E_ref, blaze::column( sigma_ref, idx_degree ) );
            sigmaSAtEnergy[idx_degree] = interp( _energies[idx_pseudotime + 1] );
        }
        for( unsigned idx_degree = 0; idx_degree <= _polyDegreeBasis; ++idx_degree ) {
            sigmaTAtEnergy[idx_degree] = ( sigmaSAtEnergy[0] - sigmaSAtEnergy[idx_degree] );
        }

        Matrix Sigma     = 0.0 * _S;
        unsigned counter = 0;
        for( int l = 0; l <= _polyDegreeBasis; ++l ) {
            for( int m = -l; m <= l; ++m ) {
                Sigma( counter, counter ) = sigmaTAtEnergy[l];
                counter++;
            }
        }
// add scattering term implicitly
#pragma omp parallel for
        for( unsigned j = 0; j < _nCells; ++j ) {
            if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
            Vector u = _M * _sol[j];
            counter  = 0;
            for( int l = 0; l <= _polyDegreeBasis; ++l ) {
                for( int m = -l; m <= l; ++m ) {
                    u[counter] = u[counter] / ( 1.0 + _dE * sigmaTAtEnergy[l] );
                    counter++;
                }
            }
            _sol[j] = _O * u;
            //_sol[j] = blaze::solve( _identity + _dE * _O * Sigma * _M, _sol[j] );
        }
    }
}

void CSDSolverTrafoFPSH2D::IterPostprocessing( unsigned idx_pseudotime ) {
    unsigned n = idx_pseudotime;
    // --- Update Solution ---
    // for( unsigned j = 0; j < _nCells; ++j ) {
    //    if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
    //    _sol[j] = _solNew[j];
    //}

    // -- Compute Dose
    for( unsigned j = 0; j < _nCells; ++j ) {
        _fluxNew[j] = dot( _sol[j], _weights );
        if( n > 0 && n < _nEnergies - 1 ) {
            _dose[j] += _dE * ( _fluxNew[j] * _sMid[n] ) / _density[j];    // update dose with trapezoidal rule // diss Kerstin
        }
        else {
            _dose[j] += 0.5 * _dE * ( _fluxNew[j] * _sMid[n] ) / _density[j];
        }
        _flux[j] = _fluxNew[j];
    }

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();
}

void CSDSolverTrafoFPSH2D::SolverPreprocessing() {
    auto log = spdlog::get( "event" );

    _densityMin = 0.1;
    for( unsigned j = 0; j < _nCells; ++j ) {
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
