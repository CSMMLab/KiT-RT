#include "solvers/csdpnsolver.hpp"
#include "common/config.hpp"
#include "common/mesh.hpp"
#include "fluxes/numericalflux.hpp"
#include "problems/problembase.hpp"
#include "solvers/csdpn_starmap_constants.hpp"
#include "toolboxes/interpolation.hpp"

// externals
#include <fstream>

#include "spdlog/spdlog.h"

CSDPNSolver::CSDPNSolver( Config* settings ) : PNSolver( settings ) {
    // --- Initialize Dose
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

CSDPNSolver::~CSDPNSolver() {}

void CSDPNSolver::IterPreprocessing( unsigned idx_iter ) {
    if( _reconsOrder > 1 ) {
        VectorVector solDivRho = _sol;
        for( unsigned j = 0; j < _nCells; ++j ) {
            solDivRho[j] = _sol[j] / _density[j];
        }
        _mesh->ComputeSlopes( _nSystem, _solDx, _solDy, solDivRho );
        _mesh->ComputeLimiter( _nSystem, _solDx, _solDy, solDivRho, _limiter );
    }

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

void CSDPNSolver::SolverPreprocessing() {}

void CSDPNSolver::IterPostprocessing( unsigned idx_iter ) {
    // --- Update Solution ---
    //_sol = _solNew;

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();

    // -- Compute Dose
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( idx_iter > 0 && idx_iter < _nEnergies - 1 ) {
            _dose[j] += _dT * ( _sol[j][0] * _sMid[idx_iter] ) / _density[j];    // update dose with trapezoidal rule // diss Kerstin
        }
        else {
            _dose[j] += 0.5 * _dT * ( _sol[j][0] * _sMid[idx_iter] ) / _density[j];
        }
    }
}

void CSDPNSolver::FluxUpdate() {
    // Loop over all spatial cells
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        Vector solL( _nSystem, 0.0 );
        Vector solR( _nSystem, 0.0 );
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;

        // Reset temporary variable psiNew
        for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {
            _solNew[idx_cell][idx_sys] = 0.0;
        }

        // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
        for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); idx_neighbor++ ) {

            // Compute flux contribution and store in psiNew to save memory
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                _solNew[idx_cell] += _g->FluxXZ( _AxPlus,
                                                 _AxMinus,
                                                 _AyPlus,
                                                 _AyMinus,
                                                 _AzPlus,
                                                 _AzMinus,
                                                 _sol[idx_cell] / _density[idx_cell],
                                                 _sol[idx_cell] / _density[idx_cell],
                                                 _normals[idx_cell][idx_neighbor] );
            else {
                unsigned int nbr_glob = _neighbors[idx_cell][idx_neighbor];    // global idx of neighbor cell
                switch( _reconsOrder ) {
                    // first order solver
                    case 1:
                        _solNew[idx_cell] +=
                            _g->FluxXZ( _AxPlus,
                                        _AxMinus,
                                        _AyPlus,
                                        _AyMinus,
                                        _AzPlus,
                                        _AzMinus,
                                        _sol[idx_cell] * ( 1.0 / _density[idx_cell] ),
                                        _sol[_neighbors[idx_cell][idx_neighbor]] * ( 1.0 / _density[_neighbors[idx_cell][idx_neighbor]] ),
                                        _normals[idx_cell][idx_neighbor] );
                        break;
                    // second order solver
                    case 2:
                        // left status of interface
                        for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {
                            solL[idx_sys] =
                                _sol[idx_cell][idx_sys] / _density[idx_cell] +
                                _limiter[idx_cell][idx_sys] *
                                    ( _solDx[idx_cell][idx_sys] * ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[idx_cell][0] ) +
                                      _solDy[idx_cell][idx_sys] * ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[idx_cell][1] ) );
                            solR[idx_sys] =
                                _sol[nbr_glob][idx_sys] / _density[nbr_glob] +
                                _limiter[nbr_glob][idx_sys] *
                                    ( _solDx[nbr_glob][idx_sys] * ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[nbr_glob][0] ) +
                                      _solDy[nbr_glob][idx_sys] * ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[nbr_glob][1] ) );
                        }
                        // flux evaluation
                        _solNew[idx_cell] +=
                            _g->FluxXZ( _AxPlus, _AxMinus, _AyPlus, _AyMinus, _AzPlus, _AzMinus, solL, solR, _normals[idx_cell][idx_neighbor] );
                        break;
                    // default: first order solver
                    default:
                        _solNew[idx_cell] += _g->FluxXZ( _AxPlus,
                                                         _AxMinus,
                                                         _AyPlus,
                                                         _AyMinus,
                                                         _AzPlus,
                                                         _AzMinus,
                                                         _sol[idx_cell] / _density[idx_cell],
                                                         _sol[_neighbors[idx_cell][idx_neighbor]] / _density[_neighbors[idx_cell][idx_neighbor]],
                                                         _normals[idx_cell][idx_neighbor] );
                }
            }
        }
    }
}

void CSDPNSolver::FVMUpdate( unsigned idx_energy ) {
    bool implicitScattering = true;
    // transform energy difference
    _dT = fabs( _eTrafo[idx_energy + 1] - _eTrafo[idx_energy] );

    // loop over all spatial cells
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        //  Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
            for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                int idx_sys = GlobalIndex( idx_l, idx_k );
                if( implicitScattering ) {
                    _solNew[idx_cell][idx_sys] =
                        _sol[idx_cell][idx_sys] - ( _dT / _areas[idx_cell] ) * _solNew[idx_cell][idx_sys];            /* cell averaged flux */
                    _solNew[idx_cell][idx_sys] = _solNew[idx_cell][idx_sys] / ( 1.0 + _dT * _sigmaTAtEnergy[idx_l] ); /* implicit scattering */
                }
                else {
                    _solNew[idx_cell][idx_sys] = _sol[idx_cell][idx_sys] -
                                                 ( _dT / _areas[idx_cell] ) * _solNew[idx_cell][idx_sys]   /* cell averaged flux */
                                                 - _dT * _sol[idx_cell][idx_sys] * _sigmaTAtEnergy[idx_l]; /* scattering */
                }
            }
            // Source Term
            // _solNew[idx_cell][0] += _dE * _Q[0][idx_cell][0];
        }
    }
}

void CSDPNSolver::PrepareVolumeOutput() {
    // std::cout << "Prepare Volume Output...";
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
            case MOMENTS:
                // As many entries as there are moments in the system
                _outputFields[idx_group].resize( _nSystem );
                _outputFieldNames[idx_group].resize( _nSystem );

                for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
                    for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                        _outputFields[idx_group][GlobalIndex( idx_l, idx_k )].resize( _nCells );

                        _outputFieldNames[idx_group][GlobalIndex( idx_l, idx_k )] =
                            std::string( "u_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
                    }
                }
                break;
            default: ErrorMessages::Error( "Volume Output Group not defined for CSD PN Solver!", CURRENT_FUNCTION ); break;
        }
    }
}

void CSDPNSolver::WriteVolumeOutput( unsigned idx_pseudoTime ) {
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();
    double maxDose;
    if( ( _settings->GetVolumeOutputFrequency() != 0 && idx_pseudoTime % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
        ( idx_pseudoTime == _nIter - 1 ) /* need sol at last iteration */ ) {

        for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
            switch( _settings->GetVolumeOutput()[idx_group] ) {
                case MINIMAL:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = _scalarFluxNew[idx_cell];
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

                default: ErrorMessages::Error( "Volume Output Group not defined for CSD PN Solver!", CURRENT_FUNCTION ); break;
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

// double CSDPNSolver::NormPDF( double x, double mu, double sigma ) {
//     return INV_SQRT_2PI / sigma * std::exp( -( ( x - mu ) * ( x - mu ) ) / ( 2.0 * sigma * sigma ) );
// }

// Vector CSDPNSolver::Time2Energy( const Vector& t, const double E_CutOff ) {
//     Interpolation interp( E_trans, E_tab );
//     Interpolation interp2( E_tab, E_trans );
//     return blaze::max( 0, interp( interp2( E_CutOff, 0 ) - t ) );
// }
//
// double CSDPNSolver::Time2Energy( const double t, const double E_CutOff ) {
//     Interpolation interp( E_trans, E_tab );
//     Interpolation interp2( E_tab, E_trans );
//     return std::fmax( 0.0, interp( E_CutOff - t ) );
// }
//
// Vector CSDPNSolver::Energy2Time( const Vector& E, const double E_CutOff ) {
//     Interpolation interp( E_tab, E_trans );
//     return blaze::max( 0, interp( E_CutOff - E ) );
// }
//
// double CSDPNSolver::Energy2Time( const double E, const double E_CutOff ) {
//     Interpolation interp( E_tab, E_trans );
//     return std::fmax( 0, interp( E_CutOff - E ) );
// }
