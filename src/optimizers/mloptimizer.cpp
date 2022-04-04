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

    // construct support structures
    SphericalBase* tempBase  = SphericalBase::Create( _settings );
    _nSystem                 = tempBase->GetBasisSize();
    VectorVector momentBasis = VectorVector( _nq, Vector( _nSystem, 0.0 ) );
    double my, phi;
    VectorVector quadPointsSphere = _quadrature->GetPointsSphere();
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        my                    = quadPointsSphere[idx_quad][0];
        phi                   = quadPointsSphere[idx_quad][1];
        momentBasis[idx_quad] = tempBase->ComputeSphericalBasis( my, phi );
    }
    _reducedMomentBasis = VectorVector( _nq, Vector( _nSystem - 1, 0.0 ) );
    for( unsigned idx_nq = 0; idx_nq < _nq; idx_nq++ ) {    // copy (reduced) moments
        for( unsigned idx_sys = 1; idx_sys < _nSystem; idx_sys++ ) {
            _reducedMomentBasis[idx_nq][idx_sys - 1] = momentBasis[idx_nq][idx_sys];
        }
    }
    delete tempBase;

    // Choose the right model depending on spherical basis, basis degree and spatial dimension
    std::string polyDegreeStr = std::to_string( _settings->GetMaxMomentDegree() );
    std::string dimStr        = std::to_string( _settings->GetDim() );
    std::string modelMkStr    = std::to_string( _settings->GetModelMK() );
    std::string basisTypeStr;
    switch( _settings->GetSphericalBasisName() ) {
        case SPHERICAL_HARMONICS: basisTypeStr = "Harmonic"; break;
        case SPHERICAL_MONOMIALS: basisTypeStr = "Monomial"; break;
    }
    // std::string modelFolder = TENSORFLOW_MODEL_PATH;
    std::string tfModelName = basisTypeStr + "_Mk" + modelMkStr + "_M" + polyDegreeStr + "_" + dimStr + "D";
    if( _settings->GetNeuralModelGamma() > 0 ) {
        tfModelName += +"_gamma" + std::to_string( _settings->GetNeuralModelGamma() );
    }
    std::string tfModelPath = TENSORFLOW_MODEL_PATH;
    tfModelPath += "/" + tfModelName;
    std::cout << "Load Tensorflow model from:\n ";
    std::cout << tfModelPath << "\n";

    // Load model
    _tfModel = new cppflow::model( tfModelPath );                                // load model
    _modelServingVectorU.resize( _settings->GetNCells() * ( _nSystem - 1 ) );    // reserve size for model servitor

    // Specify model input name
    // Call Model (change call depending on model mk) (Seems to be randomly assigned by tensorflow)
    _tfModelInputName = "";
    // std::cout << _tfModelName << "\n";
    if( tfModelName.compare( "Monomial_Mk11_M1_2D" ) == 0 || tfModelName.compare( "Monomial_Mk12_M1_2D" ) == 0 ||
        tfModelName.compare( "Monomial_Mk12_M1_2D_gamma3" ) == 0 || tfModelName.compare( "Monomial_Mk13_M3_2D" ) == 0 ) {
        _tfModelInputName = "serving_default_input_1:0";
    }
    else if( tfModelName.compare( "Monomial_Mk12_M1_2D_gamma1" ) == 0 || tfModelName.compare( "Monomial_Mk12_M1_2D_gamma2" ) == 0 ) {
        _tfModelInputName = "serving_default_x:0";
    }
    else {
        ErrorMessages::Error( "Model input name unknown. Use Tensorflow CLI to determine input name and add it to source code", CURRENT_FUNCTION );
    }
}

MLOptimizer::~MLOptimizer() {}

void MLOptimizer::Solve( Vector& alpha, Vector& u, const VectorVector& /*moments*/, unsigned /*idx_cell*/ ) {}

void MLOptimizer::SolveMultiCell( VectorVector& alpha, VectorVector& u, const VectorVector& /*moments*/ ) {

#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
        if( abs( u[idx_cell][0] ) > 1e-9 ) {
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

    // Call model
    std::vector<cppflow::tensor> output = _tfModel->operator()(
        { { _tfModelInputName, _modelInput } }, { "StatefulPartitionedCall:0", "StatefulPartitionedCall:1", "StatefulPartitionedCall:2" } );

    // reform output to vector<float>
    _modelServingVectorAlpha = output[1].get_data<float>();
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
            integral += _entropy->EntropyPrimeDual( dot( alphaRed, _reducedMomentBasis[idx_quad] ) ) * _weights[idx_quad];
        }
        alpha[idx_cell][0] = -log( integral );    // log trafo
        // rescale alpha_0 ==> done by normalized MN Solver
        // alpha[idx_cell][0] += log( u[idx_cell][0] );
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
