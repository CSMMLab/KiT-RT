/*!
 * @file newtonoptimizer.cpp
 * @brief class for solving the minimal entropy optimization problem using a neural network
 * @author S. Schotthöfer
 */

#include "optimizers/neuralnetworkoptimizer.hpp"
#include "common/config.hpp"
#include "toolboxes/errormessages.hpp"

// Only build optimizer, if tensorflow backend is enabled
#ifdef BUILD_ML
#include "entropies/entropybase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"
#include "velocitybasis/sphericalbase.hpp"

#include <iostream>

NeuralNetworkOptimizer::NeuralNetworkOptimizer( Config* settings ) : OptimizerBase( settings ) {

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
    tfModelPath += "/" + tfModelName + "/best_model";
    auto log = spdlog::get( "event" );
    log->info( "Load Tensorflow model from:\n " + tfModelPath + "\n" );

    // Load model
    _tfModel             = new cppflow::model( tfModelPath );    // load model
    unsigned servingSize = _settings->GetNCells();
    if( _settings->GetEnforceNeuralRotationalSymmetry() ) {
        if( _settings->GetMaxMomentDegree() > 2 ) {
            ErrorMessages::Error( "This postprocessing step is currently only for M1 and M2 models available.", CURRENT_FUNCTION );
        }
        servingSize *= 2;    // Double number of vectors, since we mirror the rotated vector
        _rotationMats.resize( _settings->GetNCells() );
        _rotationMatsT.resize( _settings->GetNCells() );
    }
    _modelServingVectorU.resize( servingSize * ( _nSystem - 1 ) );    // reserve size for model servitor

    // Specify model input name
    // Call Model (change call depending on model mk) (Seems to be randomly assigned by tensorflow)
    _tfModelInputName                            = "";
    bool model_found                             = false;
    std::vector<std::string> input_name_1_models = {
        "Monomial_Mk11_M1_2D",        "Monomial_Mk12_M1_2D",        "Monomial_Mk12_M1_2D_gamma3", "Monomial_Mk11_M2_2D",
        "Monomial_Mk11_M2_2D_gamma1", "Monomial_Mk11_M2_2D_gamma2", "Monomial_Mk11_M2_2D_gamma3", "Monomial_Mk12_M2_2D",
        "Monomial_Mk12_M2_2D_gamma1", "Monomial_Mk12_M2_2D_gamma2", "Monomial_Mk12_M2_2D_gamma3", "Monomial_Mk11_M3_2D",
        "Monomial_Mk11_M3_2D_gamma1", "Monomial_Mk11_M3_2D_gamma2", "Monomial_Mk11_M3_2D_gamma3", "Monomial_Mk12_M3_2D",
        "Monomial_Mk12_M3_2D_gamma1", "Monomial_Mk12_M3_2D_gamma2", "Monomial_Mk12_M3_2D_gamma3", "Monomial_Mk13_M3_2D",
        "Monomial_Mk13_M3_2D_gamma1", "Monomial_Mk13_M3_2D_gamma2", "Monomial_Mk13_M3_2D_gamma3", "Monomial_Mk11_M4_2D",
        "Monomial_Mk11_M4_2D_gamma1", "Monomial_Mk11_M4_2D_gamma2", "Monomial_Mk11_M4_2D_gamma3", "Monomial_Mk12_M4_2D",
        "Monomial_Mk12_M4_2D_gamma1", "Monomial_Mk12_M4_2D_gamma2", "Monomial_Mk12_M4_2D_gamma3", "Monomial_Mk11_M5_2D",
        "Monomial_Mk11_M5_2D_gamma1", "Monomial_Mk11_M5_2D_gamma2", "Monomial_Mk11_M5_2D_gamma3", "Monomial_Mk12_M5_2D",
        "Monomial_Mk12_M5_2D_gamma1", "Monomial_Mk12_M5_2D_gamma2", "Monomial_Mk12_M5_2D_gamma3" };
    for( std::vector<std::string>::iterator it = input_name_1_models.begin(); it != input_name_1_models.end(); ++it ) {
        if( tfModelName.compare( *it ) == 0 ) {
            model_found       = true;
            _tfModelInputName = "serving_default_input_1:0";
            break;
        }
    }
    std::vector<std::string> input_name_2_models = { "Monomial_Mk12_M1_2D_gamma1", "Monomial_Mk12_M1_2D_gamma2" };
    for( std::vector<std::string>::iterator it = input_name_2_models.begin(); it != input_name_2_models.end(); ++it ) {
        if( tfModelName.compare( *it ) == 0 ) {
            model_found       = true;
            _tfModelInputName = "serving_default_x:0";
            break;
        }
    }
    if( !model_found ) {
        ErrorMessages::Error( "Model input name <" + tfModelName + "> unknown. Use Tensorflow CLI to determine input name and add it to source code",
                              CURRENT_FUNCTION );
    }
}

NeuralNetworkOptimizer::~NeuralNetworkOptimizer() {}

void NeuralNetworkOptimizer::Solve( Vector& alpha, Vector& u, const VectorVector& /*moments*/, unsigned /*idx_cell*/ ) {}

Matrix NeuralNetworkOptimizer::CreateRotator( const Vector& uFirstMoment ) {
    double a = uFirstMoment[0];
    double b = uFirstMoment[1];
    double c, s, r;

    r = norm( uFirstMoment );    // sqrt( a * a + b * b );
    c = a / r;
    s = -b / r;

    return Matrix{ { c, -s }, { s, c } };    // Rotation Matrix
}

Vector NeuralNetworkOptimizer::RotateM1( Vector& vec, Matrix& R ) { return R * vec; }

Matrix NeuralNetworkOptimizer::RotateM2( Matrix& vec, Matrix& R, Matrix& Rt ) { return R * vec * Rt; }

void NeuralNetworkOptimizer::SolveMultiCell( VectorVector& alpha, VectorVector& u, const VectorVector& /*moments*/, Vector& alpha_norms ) {

    unsigned servingSize = _settings->GetNCells();
    Matrix rot180{ { -1.0, 0.0 }, { 0.0, -1.0 } };

    if( _settings->GetEnforceNeuralRotationalSymmetry() ) {    // Rotation Preprocessing
        if( _settings->GetMaxMomentDegree() == 1 ) {
#pragma omp parallel for
            for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
                Vector u1{ u[idx_cell][1], u[idx_cell][2] };
                _rotationMats[idx_cell]  = CreateRotator( u1 );
                _rotationMatsT[idx_cell] = blaze::trans( _rotationMats[idx_cell] );

                u1 = RotateM1( u1, _rotationMats[idx_cell] );    //, _rotationMatsT[idx_cell] );
                // Save rotated moment
                _modelServingVectorU[idx_cell * ( _nSystem - 1 )]     = (float)( u1[0] );
                _modelServingVectorU[idx_cell * ( _nSystem - 1 ) + 1] = (float)( u1[1] );    // should be zero

                // Rotate Moment by 180 degrees and save mirrored moment
                RotateM1( u1, rot180 );
                _modelServingVectorU[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 )] = (float)( u1[0] );    // Only first moment is mirrored
                _modelServingVectorU[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + 1] = (float)( u1[1] );    // should be zero
            }
        }
        else {
            // Rotate everything to x-Axis of first moment tensor
            //std::cout << "rotation_active\n";
#pragma omp parallel for
            for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
                Vector u1{ u[idx_cell][1], u[idx_cell][2] };
                Matrix u2{ { u[idx_cell][3], u[idx_cell][4] }, { u[idx_cell][4], u[idx_cell][5] } };

                _rotationMats[idx_cell]  = CreateRotator( u1 );
                _rotationMatsT[idx_cell] = blaze::trans( _rotationMats[idx_cell] );

                u1 = RotateM1( u1, _rotationMats[idx_cell] );
                u2 = RotateM2( u2, _rotationMats[idx_cell], _rotationMatsT[idx_cell] );

                _modelServingVectorU[idx_cell * ( _nSystem - 1 )]     = (float)( u1[0] );
                _modelServingVectorU[idx_cell * ( _nSystem - 1 ) + 1] = (float)( u1[1] );    // should be zero
                _modelServingVectorU[idx_cell * ( _nSystem - 1 ) + 2] = (float)( u2( 0, 0 ) );
                _modelServingVectorU[idx_cell * ( _nSystem - 1 ) + 3] = (float)( u2( 0, 1 ) );
                _modelServingVectorU[idx_cell * ( _nSystem - 1 ) + 4] = (float)( u2( 1, 1 ) );

                // Rotate Moment by 180 degrees and save mirrored moment
                u1 = RotateM1( u1, rot180 );
                u2 = RotateM2( u2, rot180, rot180 );                                                                  // mirror matrix is symmetri
                _modelServingVectorU[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 )] = (float)( u1[0] );    // Only first moment is mirrored
                _modelServingVectorU[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + 1] = (float)( u1[1] );    // should be zero

                // Even Moments cancel rotation
                _modelServingVectorU[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + 2] = (float)( u2( 0, 0 ) );
                _modelServingVectorU[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + 3] = (float)( u2( 0, 1 ) );
                _modelServingVectorU[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + 4] = (float)( u2( 1, 1 ) );
            }
        }
        servingSize *= 2;
    }
    else {    // No Postprocessing
#pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
            for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                _modelServingVectorU[idx_cell * ( _nSystem - 1 ) + idx_sys] = (float)( u[idx_cell][idx_sys + 1] );
            }
        }
    }

    // Create tensor from flattened vector
    _modelInput = cppflow::tensor( _modelServingVectorU, { servingSize, _nSystem - 1 } );

    // Call model
    std::vector<cppflow::tensor> output = _tfModel->operator()(
        { { _tfModelInputName, _modelInput } }, { "StatefulPartitionedCall:0", "StatefulPartitionedCall:1", "StatefulPartitionedCall:2" } );

    // reform output to vector<float>
    _modelServingVectorAlpha = output[1].get_data<float>();

    // Postprocessing
    if( _settings->GetEnforceNeuralRotationalSymmetry() ) {    // Rotational postprocessing
        if( _settings->GetMaxMomentDegree() == 1 ) {
#pragma omp parallel for
            for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
                Vector alphaRed       = Vector( _nSystem - 1, 0.0 );    // local reduced alpha
                Vector alphaRedMirror = Vector( _nSystem - 1, 0.0 );    // local reduced mirrored alpha
                for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                    alphaRed[idx_sys]       = (double)_modelServingVectorAlpha[idx_cell * ( _nSystem - 1 ) + idx_sys];
                    alphaRedMirror[idx_sys] = (double)_modelServingVectorAlpha[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + idx_sys];
                    // Mirror order 1 Moments
                    alphaRedMirror[idx_sys] =
                        -1 * (double)_modelServingVectorAlpha[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + idx_sys];

                    alphaRed[idx_sys] = ( alphaRed[idx_sys] + alphaRedMirror[idx_sys] ) / 2;    // average (and store in alphaRed)
                }

                // Rotate Back
                alphaRed = RotateM1( alphaRed, _rotationMatsT[idx_cell] );

                // Restore alpha_0
                double integral = 0.0;
                for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                    integral += _entropy->EntropyPrimeDual( dot( alphaRed, _reducedMomentBasis[idx_quad] ) ) * _weights[idx_quad];
                }
                alpha[idx_cell][0] = -log( integral );    // log trafo

                // Store output
                for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                    alpha[idx_cell][idx_sys + 1] = alphaRed[idx_sys];
                }
            }
        }
        else {    // Degree 2
#pragma omp parallel for
            for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
                Vector alphaRed       = Vector( _nSystem - 1, 0.0 );    // local reduced alpha
                Vector alphaRedMirror = Vector( _nSystem - 1, 0.0 );    // local reduced mirrored alpha
                for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                    alphaRed[idx_sys]       = (double)_modelServingVectorAlpha[idx_cell * ( _nSystem - 1 ) + idx_sys];
                    alphaRedMirror[idx_sys] = (double)_modelServingVectorAlpha[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + idx_sys];
                    if( idx_sys < 2 ) {    // Miror order 1 moments back
                        alphaRedMirror[idx_sys] =
                            -1 * (double)_modelServingVectorAlpha[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + idx_sys];
                    }
                    alphaRed[idx_sys] = ( alphaRed[idx_sys] + alphaRedMirror[idx_sys] ) / 2;    // average (and store in alphaRed)
                }
                alpha_norms[idx_cell] = norm( alphaRed ) * norm( alphaRed );

                // Rotate Back
                Vector alpha1{ alphaRed[0], alphaRed[1] };
                Matrix alpha2{ { alphaRed[2], alphaRed[3] * 0.5 }, { alphaRed[3] * 0.5, alphaRed[4] } };

                alpha1 = RotateM1( alpha1, _rotationMatsT[idx_cell] );                             // Rotate Back
                alpha2 = RotateM2( alpha2, _rotationMatsT[idx_cell], _rotationMats[idx_cell] );    // Rotate Back
                // Store back-rotated alpha
                alphaRed[0] = alpha1[0];
                alphaRed[1] = alpha1[1];
                alphaRed[2] = alpha2( 0, 0 );
                alphaRed[3] = 2 * alpha2( 1, 0 );
                alphaRed[4] = alpha2( 1, 1 );

                // Restore alpha_0
                double integral = 0.0;
                for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                    integral += _entropy->EntropyPrimeDual( dot( alphaRed, _reducedMomentBasis[idx_quad] ) ) * _weights[idx_quad];
                }
                alpha[idx_cell][0] = -log( integral );    // log trafo

                // Store output
                for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                    alpha[idx_cell][idx_sys + 1] = alphaRed[idx_sys];
                }
            }
        }
    }
    else {
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
        }
    }
    _modelServingVectorAlpha.clear();
}

void NeuralNetworkOptimizer::ReconstructMoments( Vector& sol, const Vector& alpha, const VectorVector& moments ) {
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
NeuralNetworkOptimizer::NeuralNetworkOptimizer( Config* settings ) : OptimizerBase( settings ) {
    ErrorMessages::Error( "ML build not configured. Please activate cmake flage BUILD_ML.", CURRENT_FUNCTION );
}
NeuralNetworkOptimizer::~NeuralNetworkOptimizer() {}
#endif
