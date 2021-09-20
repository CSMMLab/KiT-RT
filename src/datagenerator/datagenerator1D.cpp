/*!
 * \file datagenerator1D.cpp
 * \brief Class to generate data for the neural entropy closure
 * \author S. Schotthoefer
 */

#include "datagenerator/datagenerator1D.h"
#include "common/config.h"
#include "quadratures/quadraturebase.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/sphericalbase.h"

#include <iostream>
#include <omp.h>

DataGenerator1D::DataGenerator1D( Config* settings ) : DataGeneratorRegression( settings ) {
    ComputeMoments();

    // Initialize Training Data
    if( _settings->GetAlphaSampling() )
        ComputeSetSizeAlpha();
    else
        ComputeSetSizeU();

    _uSol     = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
    _alpha    = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
    _hEntropy = std::vector<double>( _setSize, 0.0 );
}

DataGenerator1D::~DataGenerator1D() {}

void DataGenerator1D::ComputeMoments() {
    double my, phi;
    phi = 0;    // placeholder. will not be used

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        my                     = _quadPointsSphere[idx_quad][0];
        _momentBasis[idx_quad] = _basisGenerator->ComputeSphericalBasis( my, phi );
    }
}

void DataGenerator1D::SampleSolutionU() {
    // Use necessary conditions from Monreal, Dissertation, Chapter 3.2.1, Page 26

    // --- Determine stepsizes etc ---

    if( _maxPolyDegree == 0 ) {
        // --- sample u in order 0 ---
        // u_0 = <1*psi>
        double du = _settings->GetMaxValFirstMoment() / (double)_gridSize;

        for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
            _uSol[idx_set][0] = du * idx_set;
        }
    }
    else if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS && !_settings->GetNormalizedSampling() ) {
        if( _maxPolyDegree == 1 ) {
            unsigned c = 0;
            double du  = _settings->GetMaxValFirstMoment() / (double)_gridSize;

            double u0 = _settings->GetRealizableSetEpsilonU0();
            double u1;
            double N1;    // helper
                          // for( double u0 = _settings->GetRealizableSetEpsilonU0(); u0 < _settings->GetMaxValFirstMoment(); u0 += du ) {
            while( u0 < _settings->GetMaxValFirstMoment() ) {
                u1 = -u0 + u0 * _settings->GetRealizableSetEpsilonU0();
                N1 = u1 / u0;
                while( N1 < 1 - _settings->GetRealizableSetEpsilonU0() ) {
                    _uSol[c][0] = u0;    // u0  by Monreals notation
                    _uSol[c][1] = u1;    // u1  by Monreals notation
                    u1 += du;
                    N1 = u1 / u0;
                    c++;
                }
                u0 += du;
            }
        }
        else if( _maxPolyDegree == 2 ) {
            unsigned c = 0;
            double du  = _settings->GetMaxValFirstMoment() / (double)_gridSize;

            double u0 = _settings->GetRealizableSetEpsilonU0();
            double u1, u2;
            double N1, N2;    // helper
            while( u0 < _settings->GetMaxValFirstMoment() ) {
                u1 = -u0 + u0 * _settings->GetRealizableSetEpsilonU0();
                N1 = u1 / u0;
                while( N1 < 1 - _settings->GetRealizableSetEpsilonU0() ) {
                    u2 = u1 * u1 / u0 + u0 * _settings->GetRealizableSetEpsilonU0();
                    N2 = u2 / u0;
                    while( N2 < 1 - _settings->GetRealizableSetEpsilonU0() ) {
                        _uSol[c][0] = u0;    // u0 (normalized i.e. N0) by Monreals notation
                        _uSol[c][1] = u1;    // u1 (normalized i.e. N1) by Monreals notation
                        _uSol[c][2] = u2;    // u2 (normalized i.e. N2) by Monreals notation
                        u2 += du;
                        N2 = u2 / u0;
                        c++;
                    }
                    u1 += du;
                    N1 = u1 / u0;
                }
                u0 += du;
            }
        }
        else if( _maxPolyDegree == 3 ) {
            unsigned c = 0;
            double du  = _settings->GetMaxValFirstMoment() / (double)_gridSize;
            double u0  = _settings->GetRealizableSetEpsilonU0();
            double u1, u2, u3;
            double N1, N2, N3;    // helper
            while( u0 < _settings->GetMaxValFirstMoment() ) {
                u1 = -u0 + u0 * _settings->GetRealizableSetEpsilonU0();
                N1 = u1 / u0;
                while( N1 < 1 - _settings->GetRealizableSetEpsilonU0() ) {
                    u2 = u1 * u1 / u0 + u0 * _settings->GetRealizableSetEpsilonU0();
                    N2 = u2 / u0;
                    while( u2 < u0 - _settings->GetRealizableSetEpsilonU0() ) {
                        u3 = -u2 + u0 * ( N1 + N2 ) * ( N1 + N2 ) / ( 1 + N1 ) + u0 * _settings->GetRealizableSetEpsilonU1();
                        N3 = u3 / u0;
                        while( N3 < N2 - ( N1 - N2 ) * ( N1 - N2 ) / ( 1 - N1 ) - _settings->GetRealizableSetEpsilonU1() ) {
                            _uSol[c][0] = u0;    // u0  by Monreals notation
                            _uSol[c][1] = u1;    // u1  by Monreals notation
                            _uSol[c][2] = u2;    // u2  by Monreals notation
                            _uSol[c][3] = u3;    // u3  by Monreals notation
                            c++;
                            u3 += du;
                            N3 = u3 / u0;
                        }
                        u2 += du;
                        N2 = N3 / u0;
                    }
                    u1 += du;
                    N1 = u1 / u0;
                }
                u0 += du;
            }
        }
    }
    else if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS && _settings->GetNormalizedSampling() ) {
        if( _maxPolyDegree == 1 ) {
            // Carefull: This computes only normalized moments, i.e. sampling for u_0 = 1
            unsigned c = 0;
            double dN  = 2.0 / (double)_gridSize;
            double N1  = -1.0 + _settings->GetRealizableSetEpsilonU0();
            while( N1 < 1.0 - _settings->GetRealizableSetEpsilonU0() ) {
                _uSol[c][0] = 1;     // u0 (normalized i.e. N0) by Monreals notation
                _uSol[c][1] = N1;    // u1 (normalized i.e. N1) by Monreals notation
                N1 += dN;
                c++;
            }
        }
        if( _maxPolyDegree == 2 ) {
            // Carefull: This computes only normalized moments, i.e. sampling for u_0 = 1
            unsigned c = 0;
            double N1  = -1.0 + _settings->GetRealizableSetEpsilonU0();
            double N2;
            double dN = 2.0 / (double)_gridSize;
            while( N1 < 1.0 - _settings->GetRealizableSetEpsilonU0() ) {
                N2 = N1 * N1 + _settings->GetRealizableSetEpsilonU0();
                while( N2 < 1.0 - _settings->GetRealizableSetEpsilonU0() ) {
                    _uSol[c][0] = 1;     // u0 (normalized i.e. N0) by Monreals notation
                    _uSol[c][1] = N1;    // u1 (normalized i.e. N1) by Monreals notation
                    _uSol[c][2] = N2;    // u2 (normalized i.e. N2) by Monreals notation
                    N2 += dN;
                    c++;
                }
                N1 += dN;
            }
        }
        if( _maxPolyDegree == 3 ) {
            // Carefull: This computes only normalized moments, i.e. sampling for u_0 = 1, N1=0
            unsigned c = 0;
            double N1  = 0 + _settings->GetRealizableSetEpsilonU0();
            double N2, N3;
            double dN = 1.0 / (double)_gridSize;
            N2        = N1 * N1 + _settings->GetRealizableSetEpsilonU0();
            while( N2 < 1.0 - _settings->GetRealizableSetEpsilonU0() ) {
                N3 = -N2 + ( N1 + N2 ) * ( N1 + N2 ) / ( 1 + N1 ) + _settings->GetRealizableSetEpsilonU1();
                while( N3 < N2 - ( N1 - N2 ) * ( N1 - N2 ) / ( 1 - N1 ) - _settings->GetRealizableSetEpsilonU1() ) {
                    _uSol[c][0] = 1;     // u0  by Monreals notation
                    _uSol[c][1] = N1;    // u1  by Monreals notation
                    _uSol[c][2] = N2;    // u2  by Monreals notation
                    _uSol[c][3] = N3;    // u3  by Monreals notation
                    c++;
                    N3 += dN;
                }
                N2 += dN;
            }
        }
    }
    else {
        ErrorMessages::Error( "Sampling for this configuration is not yet supported", CURRENT_FUNCTION );
    }
}

void DataGenerator1D::CheckRealizability() {
    // Todo
}

void DataGenerator1D::ComputeSetSizeU() {
    if( _maxPolyDegree == 0 ) {
    }
    else if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS && !_settings->GetNormalizedSampling() ) {
        if( _maxPolyDegree == 1 ) {
            unsigned c = 0;
            double du  = _settings->GetMaxValFirstMoment() / (double)_gridSize;

            double u0 = _settings->GetRealizableSetEpsilonU0();
            double u1;
            double N1;    // helper
            while( u0 < _settings->GetMaxValFirstMoment() ) {
                u1 = -u0 + u0 * _settings->GetRealizableSetEpsilonU0();
                N1 = u1 / u0;
                while( N1 < 1 - _settings->GetRealizableSetEpsilonU0() ) {
                    u1 += du;
                    N1 = u1 / u0;
                    c++;
                }
                u0 += du;
            }
            _setSize = c;
        }
        else if( _maxPolyDegree == 2 ) {
            unsigned c = 0;
            double du  = _settings->GetMaxValFirstMoment() / (double)_gridSize;

            double u0 = _settings->GetRealizableSetEpsilonU0();
            double u1, u2;
            double N1, N2;    // helper
            while( u0 < _settings->GetMaxValFirstMoment() ) {
                u1 = -u0 + u0 * _settings->GetRealizableSetEpsilonU0();
                N1 = u1 / u0;
                if( N1 < 1 - _settings->GetRealizableSetEpsilonU0() ) {
                    while( N1 < 1 - _settings->GetRealizableSetEpsilonU0() ) {
                        u2 = u1 * u1 / u0 + u0 * _settings->GetRealizableSetEpsilonU0();
                        N2 = u2 / u0;
                        if( N2 < 1 - _settings->GetRealizableSetEpsilonU0() ) {
                            while( N2 < 1 - _settings->GetRealizableSetEpsilonU0() ) {
                                c++;
                                u2 += du;
                                N2 = u2 / u0;
                            }
                        }
                        u1 += du;
                        N1 = u1 / u0;
                    }
                }
                u0 += du;
            }
            _setSize = c;
        }
        else if( _maxPolyDegree == 3 ) {
            unsigned c = 0;
            double du  = _settings->GetMaxValFirstMoment() / (double)_gridSize;
            double u0  = _settings->GetRealizableSetEpsilonU0();
            double u1, u2, u3;
            double N1, N2, N3;    // helper
            while( u0 < _settings->GetMaxValFirstMoment() ) {
                u1 = -u0 + u0 * _settings->GetRealizableSetEpsilonU0();
                N1 = u1 / u0;
                while( N1 < 1 - _settings->GetRealizableSetEpsilonU0() ) {
                    u2 = u1 * u1 / u0 + u0 * _settings->GetRealizableSetEpsilonU0();
                    N2 = u2 / u0;
                    while( u2 < u0 - _settings->GetRealizableSetEpsilonU0() ) {
                        u3 = -u2 + u0 * ( N1 + N2 ) * ( N1 + N2 ) / ( 1 + N1 ) + u0 * _settings->GetRealizableSetEpsilonU1();
                        N3 = u3 / u0;
                        while( N3 < N2 - ( N1 - N2 ) * ( N1 - N2 ) / ( 1 - N1 ) - _settings->GetRealizableSetEpsilonU1() ) {
                            u3 += du;
                            N3 = u3 / u0;
                            c++;
                        }
                        u2 += du;
                        N2 = u2 / u0;
                    }
                    u1 += du;
                    N1 = u1 / u0;
                }
                u0 += du;
            }
            _setSize = c;
        }
        else {
            ErrorMessages::Error( "Sampling of training data of degree higher than 3 is not yet implemented.", CURRENT_FUNCTION );
        }
    }
    else if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS && _settings->GetNormalizedSampling() ) {
        if( _maxPolyDegree == 1 ) {
            // Carefull: This computes only normalized moments, i.e. sampling for u_0 = 1
            unsigned c = 0;
            double dN  = 2.0 / (double)_gridSize;
            double N1  = -1.0 + _settings->GetRealizableSetEpsilonU0();
            while( N1 < 1.0 - _settings->GetRealizableSetEpsilonU0() ) {
                N1 += dN;
                c++;
            }
            _setSize = c;
        }
        else if( _maxPolyDegree == 2 ) {
            // Carefull: This computes only normalized moments, i.e. sampling for u_0 = 1
            unsigned c = 0;
            double N1  = -1.0 + _settings->GetRealizableSetEpsilonU0();
            double N2;
            double dN = 2.0 / (double)_gridSize;
            while( N1 < 1.0 - _settings->GetRealizableSetEpsilonU0() ) {
                N2 = N1 * N1 + _settings->GetRealizableSetEpsilonU0();
                while( N2 < 1.0 - _settings->GetRealizableSetEpsilonU0() ) {
                    c++;
                    N2 += dN;
                }
                N1 += dN;
            }
            _setSize = c;
        }
        else if( _maxPolyDegree == 3 ) {
            // Carefull: This computes only normalized moments, i.e. sampling for u_0 = 1, N1=0
            unsigned c = 0;
            double N1  = 0 + _settings->GetRealizableSetEpsilonU0();
            double N2, N3;
            double dN = 1.0 / (double)_gridSize;
            N2        = N1 * N1 + _settings->GetRealizableSetEpsilonU0();
            while( N2 < 1.0 - _settings->GetRealizableSetEpsilonU0() ) {
                N3 = -N2 + ( N1 + N2 ) * ( N1 + N2 ) / ( 1 + N1 ) + _settings->GetRealizableSetEpsilonU1();
                while( N3 < N2 - ( N1 - N2 ) * ( N1 - N2 ) / ( 1 - N1 ) - _settings->GetRealizableSetEpsilonU1() ) {
                    c++;
                    N3 += dN;
                }
                N2 += dN;
            }
            _setSize = c;
        }
        else {
            ErrorMessages::Error( "Sampling of training data of degree higher than 3 is not yet implemented.", CURRENT_FUNCTION );
        }
    }
}
