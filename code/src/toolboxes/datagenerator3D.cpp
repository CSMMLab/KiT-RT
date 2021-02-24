/*!
 * \file datagenerator3D.cpp
 * \brief Class to generate data for the neural entropy closure
 * \author S. Schotthoefer
 */

#include "toolboxes/datagenerator3D.h"
#include "common/config.h"
#include "quadratures/quadraturebase.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/sphericalbase.h"

#include <iostream>
#include <omp.h>

DataGenerator3D::DataGenerator3D( Config* settings ) : DataGeneratorBase( settings ) {
    ComputeMoments();

    // Initialize Training Data
    ComputeSetSize();

    _uSol     = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
    _alpha    = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
    _hEntropy = std::vector<double>( _setSize, 0.0 );
}

DataGenerator3D::~DataGenerator3D() {}

void DataGenerator3D::ComputeMoments() {
    double my, phi;

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        my                 = _quadPointsSphere[idx_quad][0];
        phi                = _quadPointsSphere[idx_quad][1];
        _moments[idx_quad] = _basis->ComputeSphericalBasis( my, phi );
    }
}

void DataGenerator3D::SampleSolutionU() {
    // Use necessary conditions from Monreal, Dissertation, Chapter 3.2.1, Page 26

    // --- Determine stepsizes etc ---
    double du0 = _settings->GetMaxValFirstMoment() / (double)_gridSize;

    // different processes for different
    if( _LMaxDegree == 0 ) {
        // --- sample u in order 0 ---
        // u_0 = <1*psi>

        for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
            _uSol[idx_set][0] = du0 * idx_set;
        }
    }
    else if( _LMaxDegree == 1 && _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
        // Sample points on unit sphere.
        QuadratureBase* quad = QuadratureBase::Create( _settings );
        VectorVector qpoints = quad->GetPoints();    // carthesian coordinates.
        unsigned long nq     = (unsigned long)quad->GetNq();

        // --- sample u in order 0 ---
        // u_0 = <1*psi>

#pragma omp parallel for schedule( guided )
        for( unsigned long idx_set = 0; idx_set < _gridSize; idx_set++ ) {
            unsigned long outerIdx = ( idx_set - 1 ) * ( idx_set ) / 2;    // sum over all radii up to current
            outerIdx *= nq;                                                // in each radius step, use all quad points

            //  --- sample u in order 1 ---
            /* order 1 has 3 elements. (omega_x, omega_y,  omega_z) = omega, let u_1 = (u_x, u_y, u_z) = <omega*psi>
             * Condition u_0 >= norm(u_1)   */
            double radiusU0 = du0 * ( idx_set + 1 );
            // Boundary correction
            if( radiusU0 < _settings->GetRealizableSetEpsilonU0() ) {
                radiusU0 = _settings->GetRealizableSetEpsilonU0();    // Boundary close to 0
            }

            unsigned long localIdx = 0;
            double radius          = 0.0;
            unsigned long innerIdx = 0;
            // loop over all radii
            for( unsigned long idx_subset = 0; idx_subset < idx_set; idx_subset++ ) {
                radius = du0 * ( idx_subset + 1 );    // dont use radius 0 ==> shift by one
                // Boundary correction
                if( radius < _settings->GetRealizableSetEpsilonU0() ) {
                    radius = _settings->GetRealizableSetEpsilonU0();    // Boundary close to 0
                }
                if( radius / radiusU0 > 0.95 * _settings->GetRealizableSetEpsilonU1() ) {
                    // small offset to take care of rounding errors in quadrature
                    radius = radiusU0 * 0.95 * _settings->GetRealizableSetEpsilonU1();    // "outer" boundary
                }

                localIdx = outerIdx + idx_subset * nq;

                for( unsigned long quad_idx = 0; quad_idx < nq; quad_idx++ ) {
                    innerIdx = localIdx + quad_idx;    // gives the global index

                    _uSol[innerIdx][0] = radiusU0;    // Prevent 1st order moments to be too close to the boundary
                    // scale quadpoints with radius
                    _uSol[innerIdx][1] = radius * qpoints[quad_idx][0];
                    _uSol[innerIdx][2] = radius * qpoints[quad_idx][1];
                    _uSol[innerIdx][3] = radius * qpoints[quad_idx][2];
                }
            }
        }
    }
    else {
        ErrorMessages::Error( "Sampling for this configuration is not yet supported", CURRENT_FUNCTION );
    }
}

void DataGenerator3D::CheckRealizability() {
    double epsilon = _settings->GetRealizableSetEpsilonU0();
    if( _LMaxDegree == 1 ) {
#pragma omp parallel for schedule( guided )
        for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
            double normU1 = 0.0;
            Vector u1( 3, 0.0 );
            if( _uSol[idx_set][0] < epsilon ) {
                if( std::abs( _uSol[idx_set][1] ) > 0 || std::abs( _uSol[idx_set][2] ) > 0 || std::abs( _uSol[idx_set][3] ) > 0 ) {
                    ErrorMessages::Error( "Moment not realizable [code 0]. Values: (" + std::to_string( _uSol[idx_set][0] ) + "|" +
                                              std::to_string( _uSol[idx_set][1] ) + "|" + std::to_string( _uSol[idx_set][2] ) + "|" +
                                              std::to_string( _uSol[idx_set][3] ) + ")",
                                          CURRENT_FUNCTION );
                }
            }
            else {
                u1     = { _uSol[idx_set][1], _uSol[idx_set][2], _uSol[idx_set][3] };
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

void DataGenerator3D::ComputeSetSize() {
    if( _LMaxDegree == 0 ) {
    }
    else if( _LMaxDegree == 1 && _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
        // Sample points on unit sphere.
        QuadratureBase* quad = QuadratureBase::Create( _settings );
        unsigned long nq     = (unsigned long)quad->GetNq();

        // Allocate memory.
        _setSize = nq * _gridSize * ( _gridSize - 1 ) / 2;

        delete quad;
    }
    else {
        ErrorMessages::Error( "Sampling of training data of degree higher than 1 is not yet implemented.", CURRENT_FUNCTION );
    }
}
