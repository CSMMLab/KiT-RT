#include "solvers/csdpnsolver.hpp"
#include "common/config.hpp"
#include "common/globalconstants.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "fluxes/numericalflux.hpp"
#include "kernels/scatteringkernelbase.hpp"
#include "problems/icru.hpp"
#include "problems/problembase.hpp"
#include "solvers/csdpn_starmap_constants.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"
// externals
#include "spdlog/spdlog.h"
#include <iostream>
#include <mpi.h>

#include "quadratures/quadraturebase.hpp"
#include "toolboxes/sphericalbase.hpp"

double normpdf( double x, double mu, double sigma ) {
    return INV_SQRT_2PI / sigma * std::exp( -( ( x - mu ) * ( x - mu ) ) / ( 2.0 * sigma * sigma ) );
}

Vector Time2Energy( const Vector& t, const double E_CutOff ) {
    Interpolation interp( E_trans, E_tab );
    Interpolation interp2( E_tab, E_trans );
    return blaze::max( 0, interp( interp2( E_CutOff, 0 ) - t ) );
}

double Time2Energy( const double t, const double E_CutOff ) {
    Interpolation interp( E_trans, E_tab );
    Interpolation interp2( E_tab, E_trans );
    return std::fmax( 0.0, interp( E_CutOff - t ) );
}

Vector Energy2Time( const Vector& E, const double E_CutOff ) {
    Interpolation interp( E_tab, E_trans );
    return blaze::max( 0, interp( E_CutOff - E ) );
}

double Energy2Time( const double E, const double E_CutOff ) {
    Interpolation interp( E_tab, E_trans );
    return std::fmax( 0, interp( E_CutOff - E ) );
}

CSDPNSolver::CSDPNSolver( Config* settings ) : PNSolver( settings ) {
    std::cout << "Start of constructor: E_ref = " << E_ref << std::endl;
    saveE_ref        = E_ref;
    _polyDegreeBasis = settings->GetMaxMomentDegree();

    // Limiter variables
    _solDx   = VectorVector( _nCells, Vector( _nSystem, 0.0 ) );
    _solDy   = VectorVector( _nCells, Vector( _nSystem, 0.0 ) );
    _limiter = VectorVector( _nCells, Vector( _nSystem, 0.0 ) );

    _basis = NULL;

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
    // std::cout << "eMax is " << maxE << " vs " << interpTrafoToE( eMaxTrafo ) << std::endl;
    // std::cout << "eMin is " << minE << " vs " << interpTrafoToE( eMinTrafo ) << std::endl;

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
    _s        = interpS( _energies );
    double dE = _eTrafo[2] - _eTrafo[1];
    Vector eTrafoMid( _nEnergies - 1 );
    for( unsigned n = 0; n < _nEnergies - 1; ++n ) {
        eTrafoMid[n] = _eTrafo[n] + dE / 2;
    }
    // compute Corresponding original energies at intermediate points
    Vector eMid( _nEnergies - 1 );
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        eMid[n] = interpTrafoToE( eMaxTrafo - eTrafoMid[n] );
    }
    _sMid = interpS( eMid );

    // write initial condition
    Vector pos_beam = Vector{ 0.5, 0.5 };
    _sol            = VectorVector( _nCells, Vector( _nSystem, 0.0 ) );
    double tmp      = 0.0;
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        double x            = _cellMidPoints[idx_cell][0];
        double y            = _cellMidPoints[idx_cell][1];
        const double stddev = .01;
        double f            = normpdf( x, pos_beam[0], stddev ) * normpdf( y, pos_beam[1], stddev );

        _sol[idx_cell][0] = f * StarMAPmoments[0];
        if( _sol[idx_cell][0] > tmp ) tmp = _sol[idx_cell][0];

        for( unsigned idx_sys = 1; idx_sys < _nSystem; idx_sys++ ) {

            //_sol[idx_cell][idx_sys] = f * StarMAPmoments[idx_sys];    // must be VectorVector
            _sol[idx_cell][idx_sys] = 0;    // isotropic
        }
    }

    _solNew = _sol;

    _dose = std::vector<double>( _settings->GetNCells(), 0.0 );

    _sigmaTAtEnergy = Vector( _polyDegreeBasis, 0.0 );

    std::cout << "End of constructor: E_ref = " << E_ref << std::endl;
}

CSDPNSolver::~CSDPNSolver() {
    if( _basis ) delete _basis;
}

void CSDPNSolver::IterPreprocessing( unsigned idx_iter ) {
    if( _reconsOrder > 1 ) {
        auto solDivRho = _sol;
        for( unsigned j = 0; j < _nCells; ++j ) {
            solDivRho[j] = _sol[j] / _density[j];
        }

        _mesh->ComputeSlopes( _nSystem, _solDx, _solDy, solDivRho );
        _mesh->ComputeLimiter( _nSystem, _solDx, _solDy, solDivRho, _limiter );
    }

    Vector sigmaSAtEnergy( _polyDegreeBasis );
    // compute scattering cross section at current energy
    for( unsigned idx_degree = 0; idx_degree < _polyDegreeBasis; ++idx_degree ) {
        // setup interpolation from E to sigma at degree idx_degree
        Interpolation interp( saveE_ref, blaze::column( sigma_ref, idx_degree ) );
        sigmaSAtEnergy[idx_degree] = interp( _energies[idx_iter + 1] );    // TODO: check if this is really idx_iter+1
    }
    for( unsigned idx_degree = 0; idx_degree < _polyDegreeBasis; ++idx_degree ) {
        _sigmaTAtEnergy[idx_degree] = ( sigmaSAtEnergy[0] - sigmaSAtEnergy[idx_degree] );
    }
}

void CSDPNSolver::SolverPreprocessing() {

    // cross sections do not need to be transformed to ETilde energy grid since e.g. TildeSigmaT(ETilde) = SigmaT(E(ETilde))
}

void CSDPNSolver::IterPostprocessing( unsigned idx_iter ) {
    // std::cout << "Iter Postprocessing...";
    // --- Update Solution ---
    _sol = _solNew;

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();

    unsigned n = idx_iter;
    // -- Compute Dose
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( n > 0 || n < _nEnergies - 1 ) {
            _dose[j] += _dE * ( _fluxNew[j] * _sMid[n] ) / _density[j];    // update dose with trapezoidal rule // diss Kerstin
            //_dose[j] += _dE * ( _fluxNew[j] * _s[_nEnergies - n - 1] ) / _density[j];    // update dose with explicit Euler rule // diss Kerstin
        }
        else {
            _dose[j] += 0.5 * _dE * ( _fluxNew[j] * _sMid[n] ) / _density[j];
            //_dose[j] += _dE * _fluxNew[j] * _s[_nEnergies - n - 1] / _density[j];
        }
    }
    // obj.dose.+= dE * ( uNew[:, 1] + uOUnc ) * obj.csd.SMid[n]./ obj.densityVec./ ( 1 + ( n == 1 || n == nEnergies ) );
    std::cout << "weight: " << _s[n] << " time: " << idx_iter * _dE << " energy: " << Time2Energy( idx_iter * _dE, _E_cutoff ) << " DONE."
              << std::endl;
}

void CSDPNSolver::FluxUpdate() {
    Vector solL( _nSystem, 0.0 );
    Vector solR( _nSystem, 0.0 );

    // Loop over all spatial cells
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {

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
                _solNew[idx_cell] += _g->Flux( _AxPlus,
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
                            _g->Flux( _AxPlus,
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
                            _g->Flux( _AxPlus, _AxMinus, _AyPlus, _AyMinus, _AzPlus, _AzMinus, solL, solR, _normals[idx_cell][idx_neighbor] );
                        break;
                    // default: first order solver
                    default:
                        _solNew[idx_cell] += _g->Flux( _AxPlus,
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
    _dE = fabs( _eTrafo[idx_energy + 1] - _eTrafo[idx_energy] );
// loop over all spatial cells
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
            for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                int idx_sys = GlobalIndex( idx_l, idx_k );
                if( implicitScattering ) {
                    _solNew[idx_cell][idx_sys] =
                        _sol[idx_cell][idx_sys] - ( _dE / _areas[idx_cell] ) * _solNew[idx_cell][idx_sys]; /* cell averaged flux */
                }
                else {
                    _solNew[idx_cell][idx_sys] = _sol[idx_cell][idx_sys] -
                                                 ( _dE / _areas[idx_cell] ) * _solNew[idx_cell][idx_sys]   /* cell averaged flux */
                                                 - _dE * _sol[idx_cell][idx_sys] * _sigmaTAtEnergy[idx_l]; /* scattering */
                }
            }
            // Source Term
            // _solNew[idx_cell][0] += _dE * _Q[0][idx_cell][0];
        }
    }
    // treat scattering implicitly
    if( implicitScattering ) {
        // loop over all spatial cells
        //#pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
            // Dirichlet cells stay at IC, farfield assumption
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
            for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
                for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                    int idx_sys                = GlobalIndex( idx_l, idx_k );
                    _solNew[idx_cell][idx_sys] = _solNew[idx_cell][idx_sys] / ( 1.0 + _dE * _sigmaTAtEnergy[idx_l] ); /* scattering */
                }
            }
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
            default: ErrorMessages::Error( "Volume Output Group not defined for PN Solver!", CURRENT_FUNCTION ); break;
        }
    }
    // std::cout << "DONE." << std::endl;
}

void CSDPNSolver::WriteVolumeOutput( unsigned idx_pseudoTime ) {
    // std::cout << "Write Volume Output...";
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
                case MOMENTS:
                    for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {
                        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                            _outputFields[idx_group][idx_sys][idx_cell] = _sol[idx_cell][idx_sys];
                        }
                    }
                    break;

                default: ErrorMessages::Error( "Volume Output Group not defined for CSD_PN_TRAFO Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
    // std::cout << "DONE." << std::endl;
}

Vector CSDPNSolver::ConstructFlux( unsigned ) {
    // for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
    //    flux += _moments[idx_quad] * ( _weights[idx_quad] * entropyFlux );
    //}
    return Vector( 1, 0.0 );
}
