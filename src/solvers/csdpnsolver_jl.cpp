#include "solvers/csdpnsolver_jl.hpp"
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

CSDPNSolver_JL::CSDPNSolver_JL( Config* settings ) : PNSolver( settings ) {
    //  std::cout << "Start of constructor: E_ref = " << E_ref << std::endl;
    saveE_ref        = E_ref;
    _polyDegreeBasis = settings->GetMaxMomentDegree();

    // Limiter variables
    _solDx   = VectorVector( _nCells, Vector( _nSystem, 0.0 ) );
    _solDy   = VectorVector( _nCells, Vector( _nSystem, 0.0 ) );
    _limiter = VectorVector( _nCells, Vector( _nSystem, 0.0 ) );
    _basis   = NULL;

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

    // write initial condition
    Vector pos_beam = Vector{ 0.5, 0.5 };
    _sol            = VectorVector( _nCells, Vector( _nSystem, 0.0 ) );

    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        double x            = _cellMidPoints[idx_cell][0];
        double y            = _cellMidPoints[idx_cell][1];
        const double stddev = .01;
        double f            = NormPDF( x, pos_beam[0], stddev ) * NormPDF( y, pos_beam[1], stddev );

        _sol[idx_cell][0] = f * StarMAPmoments[0];

        for( unsigned idx_sys = 1; idx_sys < _nSystem; idx_sys++ ) {
            _sol[idx_cell][idx_sys] = f * StarMAPmoments[idx_sys];    // must be VectorVector
        }
    }

    //// check normpdf
    // VectorVector testVec = VectorVector( 100, Vector( 2 + _nSystem, 0.0 ) );
    // double dx            = 1.0 / 100.0;
    // for( unsigned i = 0; i < 100; i++ ) {
    //     const double stddev = .01;
    //     testVec[i][0]       = dx * i;
    //     testVec[i][1]       = dx * i;
    //
    //    for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {
    //        testVec[i][idx_sys + 2] =
    //            normpdf( testVec[i][0], pos_beam[0], stddev ) * normpdf( testVec[i][1], pos_beam[1], stddev ) * StarMAPmoments[idx_sys];
    //    }
    //}
    // TextProcessingToolbox::PrintVectorVectorToFile( testVec, "IC_fullMoments.csv", 100, 2 + _nSystem );

    _solNew = _sol;

    _dose = std::vector<double>( _settings->GetNCells(), 0.0 );

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

    // // std::cout << "End of constructor: E_ref = " << E_ref << std::endl;
    // TextProcessingToolbox::PrintVectorToFile( _s, "stopping.csv", _nEnergies );
    // TextProcessingToolbox::PrintVectorToFile( _sMid, "_sMid.csv", _nEnergies - 1 );
    // TextProcessingToolbox::PrintVectorToFile( _energies, "energies.csv", _nEnergies );
    // TextProcessingToolbox::PrintVectorToFile( _eTrafo, "energiesTrafo.csv", _nEnergies );
    // TextProcessingToolbox::PrintMatrixToFile( _AxPlus, "AxPlus.csv", _nSystem );
    // TextProcessingToolbox::PrintMatrixToFile( _AxMinus, "AxMinus.csv", _nSystem );
    // TextProcessingToolbox::PrintMatrixToFile( _AyPlus, "AyPlus.csv", _nSystem );
    // TextProcessingToolbox::PrintMatrixToFile( _AyMinus, "AyMinus.csv", _nSystem );
}

CSDPNSolver_JL::~CSDPNSolver_JL() {
    if( _basis ) delete _basis;
}

void CSDPNSolver_JL::IterPreprocessing( unsigned idx_iter ) {
    if( _reconsOrder > 1 ) {
        auto solDivRho = _sol;
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
        Interpolation interp( saveE_ref, blaze::column( sigma_ref, idx_degree ) );
        sigmaSAtEnergy[idx_degree] = interp( _energies[idx_iter + 1] );
    }
    for( unsigned idx_degree = 0; idx_degree <= _polyDegreeBasis; ++idx_degree ) {
        _sigmaTAtEnergy[idx_degree] = ( sigmaSAtEnergy[0] - sigmaSAtEnergy[idx_degree] );
    }
}

void CSDPNSolver_JL::SolverPreprocessing() {
    // cross sections do not need to be transformed to ETilde energy grid since e.g. TildeSigmaT(ETilde) = SigmaT(E(ETilde))
}

void CSDPNSolver_JL::IterPostprocessing( unsigned idx_iter ) {
    // std::cout << "Iter Postprocessing...";
    // --- Update Solution ---
    _sol = _solNew;

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();

    unsigned n = idx_iter;
    // -- Compute Dose
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( n > 0 && n < _nEnergies - 1 ) {
            _dose[j] += _dE * ( _sol[j][0] * _sMid[n] ) / _density[j];    // update dose with trapezoidal rule // diss Kerstin
        }
        else {
            _dose[j] += 0.5 * _dE * ( _sol[j][0] * _sMid[n] ) / _density[j];
        }
    }
}

void CSDPNSolver_JL::FluxUpdate() {
    // Vector solL( _nSystem, 0.0 );
    // Vector solR( _nSystem, 0.0 );

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

void CSDPNSolver_JL::FVMUpdate( unsigned idx_energy ) {
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
#pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
            // Dirichlet cells stay at IC, farfield assumption
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
            for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
                for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                    int idx_sys = GlobalIndex( idx_l, idx_k );
                    // if( abs( _solNew[idx_cell][idx_sys] / ( 1.0 + _dE * _sigmaTAtEnergy[idx_l] ) ) > 1e-5 ) {
                    //     std::cout << "scatter at cell:" << idx_cell << " : " << _solNew[idx_cell][idx_sys] / ( 1.0 + _dE *
                    // _sigmaTAtEnergy[idx_l] )
                    //               << "\n";
                    // }
                    _solNew[idx_cell][idx_sys] = _solNew[idx_cell][idx_sys] / ( 1.0 + _dE * _sigmaTAtEnergy[idx_l] ); /* scattering */
                }
            }
        }
    }
}

void CSDPNSolver_JL::PrepareVolumeOutput() {
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
}

void CSDPNSolver_JL::WriteVolumeOutput( unsigned idx_pseudoTime ) {
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

                default: ErrorMessages::Error( "Volume Output Group not defined for CSD_PN_TRAFO Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
}

Vector CSDPNSolver_JL::ConstructFlux( unsigned ) {
    // for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
    //    flux += _moments[idx_quad] * ( _weights[idx_quad] * entropyFlux );
    //}
    return Vector( 1, 0.0 );
}

double CSDPNSolver_JL::NormPDF( double x, double mu, double sigma ) {
    return INV_SQRT_2PI / sigma * std::exp( -( ( x - mu ) * ( x - mu ) ) / ( 2.0 * sigma * sigma ) );
}

Vector CSDPNSolver_JL::Time2Energy( const Vector& t, const double E_CutOff ) {
    Interpolation interp( E_trans, E_tab );
    Interpolation interp2( E_tab, E_trans );
    return blaze::max( 0, interp( interp2( E_CutOff, 0 ) - t ) );
}

double CSDPNSolver_JL::Time2Energy( const double t, const double E_CutOff ) {
    Interpolation interp( E_trans, E_tab );
    Interpolation interp2( E_tab, E_trans );
    return std::fmax( 0.0, interp( E_CutOff - t ) );
}

Vector CSDPNSolver_JL::Energy2Time( const Vector& E, const double E_CutOff ) {
    Interpolation interp( E_tab, E_trans );
    return blaze::max( 0, interp( E_CutOff - E ) );
}

double CSDPNSolver_JL::Energy2Time( const double E, const double E_CutOff ) {
    Interpolation interp( E_tab, E_trans );
    return std::fmax( 0, interp( E_CutOff - E ) );
}
