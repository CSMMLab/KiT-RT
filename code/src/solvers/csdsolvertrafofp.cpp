#include "solvers/csdsolvertrafofp.h"
#include "common/config.h"
#include "common/io.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/problembase.h"
#include "quadratures/quadraturebase.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

CSDSolverTrafoFP::CSDSolverTrafoFP( Config* settings ) : SNSolver( settings ) {
    _dose = std::vector<double>( _settings->GetNCells(), 0.0 );

    // Set angle and energies
    _energies  = Vector( _nEnergies, 0.0 );    // equidistant
    _energyMin = 1e-4 * 0.511;
    _energyMax = 10e0;

    // write equidistant energy grid (false) or refined grid (true)
    GenerateEnergyGrid( false );

    // create 1D quadrature
    unsigned nq            = _settings->GetNQuadPoints();
    QuadratureBase* quad1D = QuadratureBase::Create( QUAD_GaussLegendre1D, nq );
    Vector w               = quad1D->GetWeights();
    VectorVector muVec     = quad1D->GetPoints();
    Vector mu( nq );
    for( unsigned k = 0; k < nq; ++k ) {
        mu[k] = muVec[k][0];
    }

    // setup Laplace Beltrami matrix L in slab geometry
    _L = Matrix( nq, nq, 0.0 );

    _FPMethod = 2;

    double DMinus = 0.0;
    double DPlus  = 0.0;
    for( unsigned k = 0; k < nq; ++k ) {
        DMinus = DPlus;
        DPlus  = DMinus - 2 * mu[k] * w[k];
        if( k > 0 ) {
            _L( k, k - 1 ) = DMinus / ( mu[k] - mu[k - 1] ) / w[k];
            _L( k, k )     = -DMinus / ( mu[k] - mu[k - 1] ) / w[k];
        }
        if( k < nq - 1 ) {
            _L( k, k + 1 ) = DPlus / ( mu[k + 1] - mu[k] ) / w[k];
            _L( k, k ) += -DPlus / ( mu[k + 1] - mu[k] ) / w[k];
        }
    }

    // Heney-Greenstein parameter
    double g = 0.8;

    // determine transport coefficients of Heney-Greenstein
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
        ICRU database( mu, _energies );
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

    _density = std::vector<double>( _nCells, 1.0 );
    // exit(EXIT_SUCCESS);
}

void CSDSolverTrafoFP::Solve() {
    auto log = spdlog::get( "event" );

    // save original energy field for boundary conditions
    auto energiesOrig = _energies;

    // setup incoming BC on left
    _sol = VectorVector( _density.size(), Vector( _settings->GetNQuadPoints(), 0.0 ) );    // hard coded IC, needs to be changed
    for( unsigned k = 0; k < _nq; ++k ) {
        if( _quadPoints[k][0] > 0 && !_RT ) _sol[0][k] = 1e5 * exp( -10.0 * pow( 1.0 - _quadPoints[k][0], 2 ) );
    }
    // hard coded boundary type for 1D testcases (otherwise cells will be NEUMANN)
    _boundaryCells[0]           = BOUNDARY_TYPE::DIRICHLET;
    _boundaryCells[_nCells - 1] = BOUNDARY_TYPE::DIRICHLET;

    // setup identity matrix for FP scattering
    Matrix identity( _nq, _nq, 0.0 );
    for( unsigned k = 0; k < _nq; ++k ) identity( k, k ) = 1.0;

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
    if( rank == 0 ) log->info( "{:10}   {:10}", "E", "dFlux" );

    // do substitution from psi to psiTildeHat (cf. Dissertation Kerstion Kuepper, Eq. 1.23)
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
    double densityMin = _density[0];
    for( unsigned j = 1; j < _nCells; ++j ) {
        if( densityMin > _density[j] ) densityMin = _density[j];
    }

    // cross sections do not need to be transformed to ETilde energy grid since e.g. TildeSigmaT(ETilde) = SigmaT(E(ETilde))

    // loop over energies (pseudo-time)
    for( unsigned n = 0; n < _nEnergies - 1; ++n ) {
        _dE = fabs( _energies[n + 1] - _energies[n] );    // is the sign correct here?

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

        _IL = identity - _beta * _L;

        // write BC for water phantom
        if( _RT ) {
            for( unsigned k = 0; k < _nq; ++k ) {
                if( _quadPoints[k][0] > 0 ) {
                    _sol[0][k] = 1e5 * exp( -200.0 * pow( 1.0 - _quadPoints[k][0], 2 ) ) *
                                 exp( -50.0 * pow( _energyMax - energiesOrig[_nEnergies - n - 1], 2 ) ) * _density[0] * _s[_nEnergies - n - 1];
                }
            }
        }

        // add FP scattering term implicitly
        for( unsigned j = 0; j < _nCells; ++j ) {
            if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
            //_sol[j] = blaze::solve( identity - _dE * _alpha2 * _L, psiNew[j] );
            _sol[j] = _IL * blaze::solve( _IL - _dE * _alpha * _L, _sol[j] );
        }

        // loop over all spatial cells
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
                        psiNew[j][i] += _g->Flux( _quadPoints[i],
                                                  _sol[j][i] / _density[j],
                                                  _sol[_neighbors[j][idx_neighbor]][i] / _density[_neighbors[j][idx_neighbor]],
                                                  _normals[j][idx_neighbor] ) /
                                        _areas[j];
                }
                // time update angular flux with numerical flux and total scattering cross section
                psiNew[j][i] = _sol[j][i] - _dE * psiNew[j][i];
            }
        }

        for( unsigned j = 0; j < _nCells; ++j ) {
            if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
            _sol[j] = psiNew[j];
        }

        for( unsigned j = 0; j < _nCells; ++j ) {
            fluxNew[j] = dot( _sol[j], _weights );
            if( n > 0 ) {
                _dose[j] += 0.5 * _dE * ( fluxNew[j] * _s[_nEnergies - n - 1] + fluxOld[j] * _s[_nEnergies - n] ) /
                            _density[j];    // update dose with trapezoidal rule
            }
            else {
                _dose[j] += _dE * fluxNew[j] * _s[_nEnergies - n - 1] / _density[j];
            }
            _solverOutput[j] = fluxNew[j];
        }

        // Save( n );
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 )
            log->info( "{:03.8f}  {:03.8f}  {:01.5e}  {:01.5e}", energiesOrig[_nEnergies - n - 1], _energies[n], _dE / densityMin, dFlux );
        if( std::isinf( dFlux ) || std::isnan( dFlux ) ) break;
    }
    Save( 1 );
}

void CSDSolverTrafoFP::Save() const {
    std::vector<std::string> fieldNames{ "dose", "normalized dose" };
    std::vector<std::vector<std::string>> fieldNamesWrapper{ fieldNames };

    std::vector<std::vector<double>> dose( 1, _dose );
    std::vector<std::vector<double>> normalizedDose( 1, _dose );
    double maxDose = *std::max_element( _dose.begin(), _dose.end() );
    for( unsigned i = 0; i < _dose.size(); ++i ) normalizedDose[0][i] /= maxDose;
    std::vector<std::vector<std::vector<double>>> results{ dose, normalizedDose };
    ExportVTK( _settings->GetOutputFile(), results, fieldNamesWrapper, _mesh );
}

void CSDSolverTrafoFP::Save( int currEnergy ) const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<std::vector<std::string>> fieldNamesWrapper{ fieldNames };

    std::vector<std::vector<double>> scalarField( 1, _solverOutput );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), results, fieldNamesWrapper, _mesh );
}

void CSDSolverTrafoFP::GenerateEnergyGrid( bool refinement ) {
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
