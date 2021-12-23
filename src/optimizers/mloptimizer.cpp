/*!
 * @file newtonoptimizer.cpp
 * @brief class for solving the minimal entropy optimization problem using a neural network
 * @author S. Schotthöfer
 */

#include "optimizers/mloptimizer.hpp"
#include "common/config.hpp"
#include "toolboxes/errormessages.hpp"

// Only build optimizer, if tensorflow backend is enabled
#ifdef BUILD_ML
#include "entropies/entropybase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/sphericalbase.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"

#include <iostream>

MLOptimizer::MLOptimizer( Config* settings ) : OptimizerBase( settings ) {

    _quadrature = QuadratureBase::Create( settings );
    _nq         = _quadrature->GetNq();
    _weights    = _quadrature->GetWeights();

    // construct input tensor
    SphericalBase* tempBase = SphericalBase::Create( _settings );
    _nSystem                = tempBase->GetBasisSize();
    delete tempBase;

    std::string modelFolder = TENSORFLOW_MODEL_PATH;

    // Choose the right model depending on spherical basis, basis degree and spatial dimension

    std::string polyDegreeStr = std::to_string( _settings->GetMaxMomentDegree() );
    std::string dimStr        = std::to_string( _settings->GetDim() );
    std::string modelMkStr    = std::to_string( _settings->GetModelMK() );
    std::string basisTypeStr;
    switch( _settings->GetSphericalBasisName() ) {
        case SPHERICAL_HARMONICS: basisTypeStr = "Harmonic"; break;
        case SPHERICAL_MONOMIALS: basisTypeStr = "Monomial"; break;
    }

    std::string tfModelPath = modelFolder + "/" + basisTypeStr + "_Mk" + modelMkStr + "_M" + polyDegreeStr + "_" + dimStr + "D";
    // std::cout << "Load Tensorflow model from:\n ";
    // std::cout << tfModelPath << "\n";

    // Load model
    // std::cout << _settings->GetNCells() << "\n";

    _tfModel = new cppflow::model( tfModelPath );

    _modelServingVectorU.resize( _settings->GetNCells() * ( _nSystem - 1 ) );
    //_modelInput = cppflow::fill( { int( _settings->GetNCells() ), int( _nSystem - 1 ) }, 1.0f );    //{ _settings->GetNCells(), _nSystem - 1 }

    // VectorVector uTest      = VectorVector( _settings->GetNCells(), Vector( _nSystem, 1.0 ) );
    // VectorVector alphaTest  = VectorVector( _settings->GetNCells(), Vector( _nSystem, 1.0 ) );
    // VectorVector momentTest = VectorVector( _nq, Vector( _nSystem, 1.0 ) );
    //
    // SolveMultiCell( uTest, alphaTest, momentTest );

    // ErrorMessages::Error( "ML build not configured. Please activate cmake flag BUILD_ML.", CURRENT_FUNCTION );
}

MLOptimizer::~MLOptimizer() {}

void MLOptimizer::Solve( Vector& alpha, Vector& u, const VectorVector& /*moments*/, unsigned /*idx_cell*/ ) {}

void MLOptimizer::SolveMultiCell( VectorVector& alpha, VectorVector& u, const VectorVector& moments ) {

    // Only for debugging.... Needs to go to constructor
    VectorVector momentsRed = VectorVector( _nq, Vector( _nSystem - 1, 0.0 ) );
#pragma omp parallel for
    for( unsigned idx_nq = 0; idx_nq < _nq; idx_nq++ ) {    // copy (reduced) moments
        for( unsigned idx_sys = 1; idx_sys < _nSystem; idx_sys++ ) {
            momentsRed[idx_nq][idx_sys - 1] = moments[idx_nq][idx_sys];
        }
    }

    // Transform VectorVector to flattened vector<float> and normalize data
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
        if( u[idx_cell][0] > 0 ) {
            //_modelServingVectorU[idx_cell * ( _nSystem - 1 ) + 0] = 0.0;
            //_modelServingVectorU[idx_cell * ( _nSystem - 1 ) + 1] = 0.5;
            for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                _modelServingVectorU[idx_cell * ( _nSystem - 1 ) + idx_sys] = (float)( u[idx_cell][idx_sys + 1] / u[idx_cell][0] );
            }
        }
        else {
            ErrorMessages::Error( "Particle Density is zero causing divide by zero error.", CURRENT_FUNCTION );
        }
    }
    // Create tensor from flattened vector
    _modelInput = cppflow::tensor( _modelServingVectorU, { _settings->GetNCells(), _nSystem - 1 } );

    // Call Model (change call depending on model mk) // Specific for MK11 2D now
    std::vector<cppflow::tensor> output = _tfModel->operator()(
        { { "serving_default_input_1:0", _modelInput } }, { "StatefulPartitionedCall:0", "StatefulPartitionedCall:1", "StatefulPartitionedCall:2" } );

    // reform output to vector<float>
    _modelServingVectorAlpha = output[1].get_data<float>();
    // std::cout << output[1] << "\n";
    //  Postprocessing
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
        Vector alphaRed = Vector( _nSystem - 1, 0.0 );    // local reduced alpha
        for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
            alphaRed[idx_sys]            = (double)_modelServingVectorAlpha[idx_cell * ( _nSystem - 1 ) + idx_sys];
            alpha[idx_cell][idx_sys + 1] = alphaRed[idx_sys];
        }
        // Restore alpha_0
        double integral = 0.0;
        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            integral += _entropy->EntropyPrimeDual( dot( alphaRed, momentsRed[idx_quad] ) ) * _weights[idx_quad];
        }
        alpha[idx_cell][0] = -log( integral );    // log trafo
        // rescale alpha_0
        alpha[idx_cell][0] += log( u[idx_cell][0] );
    }
    // postprocessing depending on model mk
    // TextProcessingToolbox::PrintVectorVector( alpha );
}

void MLOptimizer::ReconstructMoments( Vector& sol, const Vector& alpha, const VectorVector& moments ) {
    double entropyReconstruction = 0.0;
    for( unsigned idx_sys = 0; idx_sys < sol.size(); idx_sys++ ) {
        sol[idx_sys] = 0.0;
    }
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        // Make entropyReconstruction a member vector, s.t. it does not have to be re-evaluated in ConstructFlux
        entropyReconstruction = _entropy->EntropyPrimeDual( blaze::dot( alpha, moments[idx_quad] ) );
        sol += moments[idx_quad] * ( _weights[idx_quad] * entropyReconstruction );
    }
}

#else
MLOptimizer::MLOptimizer( Config* settings ) : OptimizerBase( settings ) {
    ErrorMessages::Error( "ML build not configured. Please activate cmake flage BUILD_ML.", CURRENT_FUNCTION );
}
MLOptimizer::~MLOptimizer() {}
#endif
