/*!
 * \file datagenerator3D.cpp
 * \brief Class to generate data for the neural entropy closure
 * \author S. Schotthoefer
 */

#include "datagenerator/datagenerator2D.h"
#include "common/config.h"
#include "quadratures/quadraturebase.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/sphericalbase.h"

#include <iostream>
#include <omp.h>

DataGenerator2D::DataGenerator2D( Config* settings ) : DataGeneratorRegression( settings ) {
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

DataGenerator2D::~DataGenerator2D() {}

void DataGenerator2D::ComputeMoments() {
    double my, phi;

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        my                     = _quadPointsSphere[idx_quad][0];
        phi                    = _quadPointsSphere[idx_quad][1];
        _momentBasis[idx_quad] = _basisGenerator->ComputeSphericalBasis( my, phi );
    }
}

void DataGenerator2D::SampleSolutionU() {
    // Use necessary conditions from Monreal, Dissertation, Chapter 3.2.1, Page 26

    // different processes for different
    if( _maxPolyDegree == 0 ) {
        // --- sample u in order 0 ---
        // u_0 = <1*psi>

        // --- Determine stepsizes etc ---
        double du0 = _settings->GetMaxValFirstMoment() / (double)_gridSize;

        for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
            _uSol[idx_set][0] = du0 * idx_set;
        }
    }
    else if( _settings->GetNormalizedSampling() ) {
        if( _maxPolyDegree == 1 && _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
            // Sample points on unit sphere.
            // Use MC sampling: x randUniform ==> r = sqrt(x)*u_0Max
            //                  y radnUniform ==> a = 2*pi*y,

            // Create generator
            std::default_random_engine generator;
            std::uniform_real_distribution<double> distribution( 0.0, 1.0 );

#pragma omp parallel for schedule( guided )
            for( unsigned long idx_set = 0; idx_set < _setSize; idx_set++ ) {
                double mu  = std::sqrt( distribution( generator ) );
                double phi = 2 * M_PI * distribution( generator );

                if( std::sqrt( 1 - mu * mu ) > 1 - _settings->GetRealizableSetEpsilonU0() ) {
                    idx_set--;
                    continue;
                }
                else {
                    _uSol[idx_set][0] = 1.0;
                    _uSol[idx_set][1] = std::sqrt( 1 - mu * mu ) * std::cos( phi );
                    _uSol[idx_set][2] = std::sqrt( 1 - mu * mu ) * std::sin( phi );
                }
            }
        }
        else {
            ErrorMessages::Error( "Sampling for this configuration is not yet supported", CURRENT_FUNCTION );
        }
    }
    else if( !_settings->GetNormalizedSampling() ) {
        if( _maxPolyDegree == 1 && _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
            // Sample points on unit sphere.
            // Use MC sampling: x randUniform ==> r = sqrt(x)*u_0Max
            //                  y radnUniform ==> a = 2*pi*y,

            // Create generator
            std::default_random_engine generator;
            std::uniform_real_distribution<double> distribution( 0.0, 1.0 );

            double du     = 0.001;
            long gridSize = (long)( (double)_settings->GetMaxValFirstMoment() / du );
            long c        = 0;
            for( long idx_set = 0; idx_set < gridSize; idx_set++ ) {

                double radiusU0 = du * ( idx_set + 1 );
                // Boundary correction
                if( radiusU0 < _settings->GetRealizableSetEpsilonU0() ) {
                    radiusU0 = _settings->GetRealizableSetEpsilonU0();    // Boundary close to 0
                }
                // number of points for unit circle should scale with (u_0)^2,
                long currSetSize = _gridSize;    //(long)( ( radiusU0 * radiusU0 ) * (double)_gridSize );

                for( long idx_currSet = 0; idx_currSet < currSetSize; idx_currSet++ ) {
                    double mu  = std::sqrt( distribution( generator ) );
                    double phi = 2 * M_PI * distribution( generator );
                    if( std::sqrt( 1 - mu * mu ) < _settings->GetRealizableSetEpsilonU1() &&
                        std::sqrt( 1 - mu * mu ) > radiusU0 - _settings->GetRealizableSetEpsilonU1() ) {
                        ErrorMessages::Error( "sampling set is empty. Boundaries overlap", CURRENT_FUNCTION );    // empty set
                    }
                    if( std::sqrt( 1 - mu * mu ) > 1 - _settings->GetRealizableSetEpsilonU1() * radiusU0 ) {
                        idx_currSet--;
                        continue;
                    }
                    else if( std::sqrt( 1 - mu * mu ) < _settings->GetRealizableSetEpsilonU1() * radiusU0 ) {
                        idx_currSet--;
                        continue;
                    }
                    else {
                        _uSol[c][0] = radiusU0;
                        _uSol[c][1] = radiusU0 * std::sqrt( 1 - mu * mu ) * std::cos( phi );
                        _uSol[c][2] = radiusU0 * std::sqrt( 1 - mu * mu ) * std::sin( phi );
                        c++;
                    }
                }
            }
        }
        else {
            ErrorMessages::Error( "Sampling for this configuration is not yet supported", CURRENT_FUNCTION );
        }
    }
}

void DataGenerator2D::CheckRealizability() {
    double epsilon = _settings->GetRealizableSetEpsilonU0();
    if( _maxPolyDegree == 1 ) {
        //#pragma omp parallel for schedule( guided )
        for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
            double normU1 = 0.0;
            Vector u1( 3, 0.0 );
            if( _uSol[idx_set][0] < epsilon ) {
                if( std::abs( _uSol[idx_set][1] ) > 0 || std::abs( _uSol[idx_set][2] ) > 0 ) {
                    ErrorMessages::Error( "Moment not realizable [code 0]. Values: (" + std::to_string( _uSol[idx_set][0] ) + "|" +
                                              std::to_string( _uSol[idx_set][1] ) + "|" + std::to_string( _uSol[idx_set][2] ) + ")",
                                          CURRENT_FUNCTION );
                }
            }
            else {
                u1     = { _uSol[idx_set][1], _uSol[idx_set][2] };
                normU1 = norm( u1 );
                if( normU1 / _uSol[idx_set][0] > _settings->GetRealizableSetEpsilonU1() ) {    // Danger Hardcoded
                    std::cout << "normU1 / _uSol[" << idx_set << "][0]: " << normU1 / _uSol[idx_set][0] << "\n";
                    std::cout << "normU1: " << normU1 << " | _uSol[idx_set][0] " << _uSol[idx_set][0] << "\n";
                    ErrorMessages::Error( "Moment to close to boundary of realizable set [code 1].\nBoundary ratio: " +
                                              std::to_string( normU1 / _uSol[idx_set][0] ),
                                          CURRENT_FUNCTION );
                }
                if( normU1 / _uSol[idx_set][0] <= 0 /*+ 0.5 * epsilon*/ ) {
                    // std::cout << "_uSol" << _uSol[idx_set][1] << " | " << _uSol[idx_set][2] << " | " << _uSol[idx_set][3] << " \n";
                    // std::cout << "normU1 / _uSol[" << idx_set << "][0]: " << normU1 / _uSol[idx_set][0] << "\n";
                    // std::cout << "normU1: " << normU1 << " | _uSol[idx_set][0] " << _uSol[idx_set][0] << "\n";
                    ErrorMessages::Error( "Moment to close to boundary of realizable set [code 2].\nBoundary ratio: " +
                                              std::to_string( normU1 / _uSol[idx_set][0] ),
                                          CURRENT_FUNCTION );
                }
            }
        }
    }
}

void DataGenerator2D::ComputeSetSizeU() {
    if( _maxPolyDegree == 0 ) {
    }
    else if( _settings->GetNormalizedSampling() ) {
        if( _maxPolyDegree == 1 && _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
            // Do nothing. Setsize is already given.
        }
        else {
            ErrorMessages::Error( "Sampling for this configuration is not yet supported", CURRENT_FUNCTION );
        }
    }
    else if( !_settings->GetNormalizedSampling() ) {
        if( _maxPolyDegree == 1 && _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {

            double du     = 0.001;
            long gridSize = (long)( (double)_settings->GetMaxValFirstMoment() / du );    // Hardcoded
            long c        = 0;
            for( long idx_set = 0; idx_set < gridSize; idx_set++ ) {

                double radiusU0 = du * ( idx_set + 1 );
                // Boundary correction
                if( radiusU0 < _settings->GetRealizableSetEpsilonU0() ) {
                    radiusU0 = _settings->GetRealizableSetEpsilonU0();    // Boundary close to 0
                }
                // number of points for unit circle should scale with (u_0)^2,
                long currSetSize = _gridSize;    // (long)( ( radiusU0 * radiusU0 ) * (double)_gridSize );

                for( long idx_currSet = 0; idx_currSet < currSetSize; idx_currSet++ ) {
                    c++;
                }
            }
            _setSize = c;
        }
        else {
            ErrorMessages::Error( "Sampling for this configuration is not yet supported", CURRENT_FUNCTION );
        }
    }
}
