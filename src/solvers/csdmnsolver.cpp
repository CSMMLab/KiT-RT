#include "solvers/csdmnsolver.hpp"
#include "common/config.hpp"
#include "common/mesh.hpp"
#include "entropies/entropybase.hpp"
#include "fluxes/numericalflux.hpp"
#include "optimizers/optimizerbase.hpp"
#include "solvers/csdpn_starmap_constants.hpp"
#include "toolboxes/interpolation.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"
#include "velocitybasis/sphericalbase.hpp"

// externals
#include "spdlog/spdlog.h"

CSDMNSolver::CSDMNSolver( Config* settings ) : MNSolver( settings ) {
    // --- Initialize Dose ---
    _dose = std::vector<double>( _nCells, 0.0 );

    // --- Compute transformed energy grid ---

    // determine transformed energy grid for tabulated grid
    Vector E_transformed( E_trans.size(), 0.0 );
    for( unsigned i = 1; i < E_trans.size(); ++i )
        E_transformed[i] = E_transformed[i - 1] + ( E_tab[i] - E_tab[i - 1] ) / 2 * ( 1.0 / S_tab[i] + 1.0 / S_tab[i - 1] );

    // determine minimal and maximal energies
    double minE = 5e-5;
    double maxE = _settings->GetMaxEnergyCSD();

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

    // --- evaluate corresponding stopping powers and transport coefficients

    // compute stopping powers
    Vector etmp = E_tab;
    Vector stmp = S_tab;
    Interpolation interpS( etmp, stmp );
    _sigmaTAtEnergy = Vector( _polyDegreeBasis + 1, 0.0 );

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
}

CSDMNSolver::~CSDMNSolver() {}

void CSDMNSolver::IterPreprocessing( unsigned idx_iter ) {

    // ------- Entropy closure Step ----------------

    _optimizer->SolveMultiCell( _alpha, _sol, _momentBasis );

    // ------- Solution reconstruction step ----
    //#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            // compute the kinetic density at all grid cells
            _kineticDensity[idx_cell][idx_quad] = _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_cell], _momentBasis[idx_quad] ) );
        }
        if( _settings->GetRealizabilityReconstruction() ) ComputeRealizableSolution( idx_cell );
    }

    // ------ Compute density normalized slope limiters and cell gradients ---
    if( _reconsOrder > 1 ) {
        VectorVector solDivRho = _sol;
        for( unsigned j = 0; j < _nCells; ++j ) {
            solDivRho[j] = _kineticDensity[j] / _density[j];
        }
        _mesh->ComputeSlopes( _nq, _solDx, _solDy, solDivRho );
        _mesh->ComputeLimiter( _nq, _solDx, _solDy, solDivRho, _limiter );
    }

    // ------ evaluate scatter coefficient at current energy level
    Vector sigmaSAtEnergy( _polyDegreeBasis + 1, 0.0 );
    // compute scattering cross section at current energy
    for( unsigned idx_degree = 0; idx_degree <= _polyDegreeBasis; ++idx_degree ) {
        // setup interpolation from E to sigma at degree idx_degree
        Interpolation interp( E_ref, blaze::column( sigma_ref, idx_degree ) );
        sigmaSAtEnergy[idx_degree] = interp( _energies[idx_iter + 1] );
    }
    for( unsigned idx_degree = 0; idx_degree <= _polyDegreeBasis; ++idx_degree ) {
        _sigmaTAtEnergy[idx_degree] = ( sigmaSAtEnergy[0] - sigmaSAtEnergy[idx_degree] );
    }
}

void CSDMNSolver::SolverPreprocessing() {
    // cross sections do not need to be transformed to ETilde energy grid since e.g. TildeSigmaT(ETilde) = SigmaT(E(ETilde))
}

void CSDMNSolver::IterPostprocessing( unsigned idx_iter ) {
    // --- Update Solution ---
    //_sol = _solNew;

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();

    // --- Compute Dose ---
#pragma omp parallel for
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( idx_iter > 0 && idx_iter < _nEnergies - 1 ) {
            _dose[j] += _dE * ( _sol[j][0] * _sMid[idx_iter] ) / _density[j];    // update dose with trapezoidal rule // diss Kerstin
        }
        else {
            _dose[j] += 0.5 * _dE * ( _sol[j][0] * _sMid[idx_iter] ) / _density[j];
        }
    }
}

Vector CSDMNSolver::ConstructFlux( unsigned idx_cell ) {
    //--- Integration of moments of flux ---
    double solL, solR, kineticFlux;
    Vector flux( _nSystem, 0.0 );

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        kineticFlux = 0.0;    // reset temorary flux

        for( unsigned idx_nbr = 0; idx_nbr < _neighbors[idx_cell].size(); idx_nbr++ ) {
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_nbr] == _nCells ) {
                // Boundary cells are first order and mirror ghost cells
                solL = _kineticDensity[idx_cell][idx_quad] / _density[idx_cell];
                solR = solL;
            }
            else {
                // interior cell
                unsigned int nbr_glob = _neighbors[idx_cell][idx_nbr];    // global idx of neighbor cell
                if( _reconsOrder == 1 ) {
                    solL = _kineticDensity[idx_cell][idx_quad] / _density[idx_cell];
                    solR = _kineticDensity[nbr_glob][idx_quad] / _density[nbr_glob];
                }
                else if( _reconsOrder == 2 ) {
                    solL = _kineticDensity[idx_cell][idx_quad] / _density[idx_cell] +
                           _limiter[idx_cell][idx_quad] *
                               ( _solDx[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][0] - _cellMidPoints[idx_cell][0] ) +
                                 _solDy[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][1] - _cellMidPoints[idx_cell][1] ) );
                    solR = _kineticDensity[nbr_glob][idx_quad] / _density[nbr_glob] +
                           _limiter[nbr_glob][idx_quad] *
                               ( _solDx[nbr_glob][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][0] - _cellMidPoints[nbr_glob][0] ) +
                                 _solDy[nbr_glob][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][1] - _cellMidPoints[nbr_glob][1] ) );
                }
                else {
                    ErrorMessages::Error( "Reconstruction order not supported.", CURRENT_FUNCTION );
                }
            }
            // Kinetic flux
            kineticFlux += _g->FluxXZ( _quadPoints[idx_quad], solL, solR, _normals[idx_cell][idx_nbr] );
        }
        // Solution flux
        flux += _momentBasis[idx_quad] * ( _weights[idx_quad] * kineticFlux );
    }
    return flux;
}

void CSDMNSolver::FluxUpdate() {
// Loop over the grid cells
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        // Dirichlet Boundaries stayd
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        _solNew[idx_cell] = ConstructFlux( idx_cell );
    }
}

void CSDMNSolver::FVMUpdate( unsigned idx_energy ) {
    bool implicitScattering = true;
    // transform energy difference
    _dE = fabs( _eTrafo[idx_energy + 1] - _eTrafo[idx_energy] );
    // loop over all spatial cells
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        //  Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
            for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                int idx_sys = _basis->GetGlobalIndexBasis( idx_l, idx_k );
                if( implicitScattering ) {
                    _solNew[idx_cell][idx_sys] =
                        _sol[idx_cell][idx_sys] - ( _dE / _areas[idx_cell] ) * _solNew[idx_cell][idx_sys];            /* cell averaged flux */
                    _solNew[idx_cell][idx_sys] = _solNew[idx_cell][idx_sys] / ( 1.0 + _dE * _sigmaTAtEnergy[idx_l] ); /* implicit scattering */
                }
                else {
                    _solNew[idx_cell][idx_sys] = _sol[idx_cell][idx_sys] -
                                                 ( _dE / _areas[idx_cell] ) * _solNew[idx_cell][idx_sys]   /* cell averaged flux */
                                                 - _dE * _sol[idx_cell][idx_sys] * _sigmaTAtEnergy[idx_l]; /* scattering */
                }
            }
            // Source Term TODO
            // _solNew[idx_cell][0] += _dE * _Q[0][idx_cell][0];
        }
    }
}

void CSDMNSolver::PrepareVolumeOutput() {
    // std::cout << "Prepare Volume Output...";
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();

    _outputFieldNames.resize( nGroups );
    _outputFields.resize( nGroups );

    // Prepare all OutputGroups ==> Specified in option VOLUME_OUTPUT
    for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
        // Prepare all Output Fields per group
        unsigned glob_idx = 0;

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
            case MOMENTS:
                // As many entries as there are moments in the system
                _outputFields[idx_group].resize( _nSystem );
                _outputFieldNames[idx_group].resize( _nSystem );
                if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
                    for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
                        for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                            glob_idx = _basis->GetGlobalIndexBasis( idx_l, idx_k );
                            _outputFields[idx_group][glob_idx].resize( _nCells );
                            _outputFieldNames[idx_group][glob_idx] = std::string( "u_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
                        }
                    }
                }
                else {
                    for( unsigned idx_l = 0; idx_l <= _polyDegreeBasis; idx_l++ ) {
                        unsigned maxOrder_k = _basis->GetCurrDegreeSize( idx_l );
                        for( unsigned idx_k = 0; idx_k < maxOrder_k; idx_k++ ) {
                            _outputFields[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )].resize( _nCells );
                            _outputFieldNames[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )] =
                                std::string( "u_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
                        }
                    }
                }
                break;
            case DUAL_MOMENTS:
                // As many entries as there are moments in the system
                _outputFields[idx_group].resize( _nSystem );
                _outputFieldNames[idx_group].resize( _nSystem );

                if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
                    for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
                        for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                            _outputFields[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )].resize( _nCells );
                            _outputFieldNames[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )] =
                                std::string( "alpha_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
                        }
                    }
                }
                else {
                    for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
                        unsigned maxOrder_k = _basis->GetCurrDegreeSize( idx_l );
                        for( unsigned idx_k = 0; idx_k < maxOrder_k; idx_k++ ) {
                            _outputFields[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )].resize( _nCells );
                            _outputFieldNames[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )] =
                                std::string( "alpha_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
                        }
                    }
                }
                break;
            default: ErrorMessages::Error( "Volume Output Group not defined for CSD MN Solver!", CURRENT_FUNCTION ); break;
        }
    }
}

void CSDMNSolver::WriteVolumeOutput( unsigned idx_pseudoTime ) {
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
                case MOMENTS:
                    for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {
                        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                            _outputFields[idx_group][idx_sys][idx_cell] = _sol[idx_cell][idx_sys];
                        }
                    }
                    break;
                case DUAL_MOMENTS:
                    for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {
                        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                            _outputFields[idx_group][idx_sys][idx_cell] = _alpha[idx_cell][idx_sys];
                        }
                    }
                    break;
                default: ErrorMessages::Error( "Volume Output Group not defined for CSD MN Solver!", CURRENT_FUNCTION ); break;
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

double CSDMNSolver::NormPDF( double x, double mu, double sigma ) {
    return INV_SQRT_2PI / sigma * std::exp( -( ( x - mu ) * ( x - mu ) ) / ( 2.0 * sigma * sigma ) );
}
