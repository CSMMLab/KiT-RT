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
#include <math.h>

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
        "Harmonic_Mk11_M1_2D",        "Harmonic_Mk11_M1_2D_gamma1", "Harmonic_Mk11_M1_2D_gamma2", "Harmonic_Mk11_M1_3D_gamma3",
        "Harmonic_Mk12_M1_3D",        "Harmonic_Mk12_M1_3D_gamma1", "Harmonic_Mk12_M1_3D_gamma2", "Harmonic_Mk12_M1_3D_gamma3",
        "Harmonic_Mk11_M1_3D",        "Harmonic_Mk11_M1_3D_gamma1", "Harmonic_Mk11_M1_3D_gamma2", "Harmonic_Mk11_M1_3D_gamma3",
        "Monomial_Mk11_M1_2D",        "Monomial_Mk11_M1_2D_gamma1", "Monomial_Mk11_M1_2D_gamma2", "Monomial_Mk11_M1_2D_gamma3",
        "Monomial_Mk12_M1_2D",        "Monomial_Mk12_M1_2D_gamma3", "Monomial_Mk11_M2_2D",        "Monomial_Mk11_M2_2D_gamma1",
        "Monomial_Mk11_M2_2D_gamma2", "Monomial_Mk11_M2_2D_gamma3", "Monomial_Mk12_M2_2D",        "Monomial_Mk12_M2_2D_gamma1",
        "Monomial_Mk12_M2_2D_gamma2", "Monomial_Mk12_M2_2D_gamma3", "Monomial_Mk11_M3_2D",        "Monomial_Mk11_M3_2D_gamma1",
        "Monomial_Mk11_M3_2D_gamma2", "Monomial_Mk11_M3_2D_gamma3", "Monomial_Mk12_M3_2D",        "Monomial_Mk12_M3_2D_gamma1",
        "Monomial_Mk12_M3_2D_gamma2", "Monomial_Mk12_M3_2D_gamma3", "Monomial_Mk13_M3_2D",        "Monomial_Mk13_M3_2D_gamma1",
        "Monomial_Mk13_M3_2D_gamma2", "Monomial_Mk13_M3_2D_gamma3", "Monomial_Mk11_M4_2D",        "Monomial_Mk11_M4_2D_gamma1",
        "Monomial_Mk11_M4_2D_gamma2", "Monomial_Mk11_M4_2D_gamma3", "Monomial_Mk12_M4_2D",        "Monomial_Mk12_M4_2D_gamma1",
        "Monomial_Mk12_M4_2D_gamma2", "Monomial_Mk12_M4_2D_gamma3", "Monomial_Mk11_M5_2D",        "Monomial_Mk11_M5_2D_gamma1",
        "Monomial_Mk11_M5_2D_gamma2", "Monomial_Mk11_M5_2D_gamma3", "Monomial_Mk12_M5_2D",        "Monomial_Mk12_M5_2D_gamma1",
        "Monomial_Mk12_M5_2D_gamma2", "Monomial_Mk12_M5_2D_gamma3" };
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

void NeuralNetworkOptimizer::Solve( Vector& /*alpha*/, Vector& /* u*/, const VectorVector& /*moments*/, unsigned /*idx_cell*/ ) {}

Matrix NeuralNetworkOptimizer::CreateRotator( const Vector& uFirstMoment ) {
    double a = uFirstMoment[0];
    double b = uFirstMoment[1];
    double c, s, r;

    r = norm( uFirstMoment );    // sqrt( a * a + b * b );
    c = a / r;
    s = -b / r;

    return Matrix{ { c, -s }, { s, c } };    // Rotation Matrix
}

Matrix NeuralNetworkOptimizer::CreateRotatorSphericalHarmonics(double theta, double x, double y ) {
    // Assumes that spherical harmonics degree is > 1

    double r = norm( Vector{x,y} ); // sqrt( a * a + b * b );
    double c = x / r;
    double s = -y / r;
    if (r <  std::numeric_limits<double>::epsilon()  ){
        s = 0.0;
        c = 1.0;
    }
    if (isnan(c) || isnan(s)){
        std::cout << x << " | " << y << "\n";
        ErrorMessages::Error("nan detected" ,CURRENT_FUNCTION);
    }

    //double c = cos(theta);
    //double s = sin(theta);

    double c2 = c*c - s*s;
    double s2 = 2*s*c;

    Matrix R = Matrix(_nSystem,_nSystem,0.0);

    // Build R by going through submatrices
    R(0,0) =1.0;
    if ( _settings->GetMaxMomentDegree()>=1){
        R(1,1) = c;
        R(3,1) = s;
        R(2,2) = 1.0;
        R(1,3) = -s;
        R(3,3) = c;
    }
    if ( _settings->GetMaxMomentDegree()>=2){
        R(4,4) = c2;
        R(4,8) = -s2;
        R(5,5) = c;
        R(5,7) = -s;
        R(6,6) = 1;
        R(7,5) = s;
        R(7,7) = c;
        R(8,4) = s2;
        R(8,8) = c2;
    }
    if ( _settings->GetMaxMomentDegree()>=3){
        ErrorMessages::Error("Rotation Matrix for spherical harmonics with degree >2 not yet implementd.",CURRENT_FUNCTION);
    }
    return R;    // Rotation Matrix
}

Matrix NeuralNetworkOptimizer::CreateRotatorSphericalHarmonics2D( double theta, double x, double y ) {
    // Assumes that spherical harmonics degree is > 1

    theta = std::atan2( x, y );

    double c  = cos( theta );
    double s  = sin( theta );
    double c2 = cos( 2 * theta );
    double s2 = sin( 2 * theta );
    double c3 = cos( 3 * theta );
    double s3 = sin( 3 * theta );
    double c4 = cos( 4 * theta );
    double s4 = sin( 4 * theta );

    Matrix R = Matrix( _nSystem, _nSystem, 0.0 );

    // Build R by going through submatrices
    R( 0, 0 ) = 1.0;
    if( _settings->GetMaxMomentDegree() >= 1 ) {
        R( 1, 1 ) = c;
        R( 2, 2 ) = c;
        R( 1, 2 ) = -s;
        R( 2, 1 ) = s;
    }
    if( _settings->GetMaxMomentDegree() >= 2 ) {
        R( 3, 3 ) = c2;
        R( 4, 4 ) = 1;
        R( 5, 5 ) = c2;
        // off-diagonal
        R( 3, 5 ) = -s2;
        R( 5, 3 ) = s2;
    }
    if( _settings->GetMaxMomentDegree() >= 3 ) {
        // diagonal
        R( 6, 6 ) = c3;
        R( 7, 7 ) = c;
        R( 8, 8 ) = c;
        R( 9, 9 ) = c3;
        // off-diagonal
        R( 6, 9 ) = -s3;
        R( 7, 8 ) = -s;
        R( 8, 7 ) = s;
        R( 9, 6 ) = s3;
    }
    if( _settings->GetMaxMomentDegree() >= 4 ) {
        // diagonal
        R( 10, 10 ) = c4;
        R( 11, 11 ) = c2;
        R( 12, 12 ) = 1;
        R( 13, 13 ) = c2;
        R( 14, 14 ) = c4;
        // off-diagonal
        R( 10, 14 ) = -s4;
        R( 11, 13 ) = -s2;
        R( 13, 11 ) = s2;
        R( 14, 10 ) = s4;
    }    // Further entries are sN and cN in above notation (see https://3dvar.com/Green2003Spherical.pdf)
    if( _settings->GetMaxMomentDegree() >= 5 ) {
        ErrorMessages::Error( "Rotation Matrix for spherical harmonics with degree >5 not yet implementd.", CURRENT_FUNCTION );
    }
    return R;
}

Vector NeuralNetworkOptimizer::RotateM1( Vector& vec, Matrix& R ) { return R * vec; }

Matrix NeuralNetworkOptimizer::RotateM2( Matrix& vec, Matrix& R, Matrix& Rt ) { return R * vec * Rt; }

void NeuralNetworkOptimizer::InferenceMonomial( VectorVector& alpha, VectorVector& u,const VectorVector& moments, Vector& alpha_norms ) {
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
            // std::cout << "rotation_active\n";
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
    else {    // No Preprocessing
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
                Vector alphaCorr = Vector(_nSystem,0.0);
                Vector uCorr =  Vector(_nSystem,0.0);

                for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                    alphaRed[idx_sys]       = (double)_modelServingVectorAlpha[idx_cell * ( _nSystem - 1 ) + idx_sys];
                    alphaRedMirror[idx_sys] = (double)_modelServingVectorAlpha[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + idx_sys];
                    // Mirror order 1 Moments
                    alphaRedMirror[idx_sys] =
                        -1 * (double)_modelServingVectorAlpha[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + idx_sys];
                    alphaRed[idx_sys] = ( alphaRed[idx_sys] + alphaRedMirror[idx_sys] ) / 2;    // average (and store in alphaRed)
                    alphaCorr[idx_sys+1] = alphaRed[idx_sys];
                }

                // Correct rotation offset made by nn approximation
                {
                // Restore alpha_0
                double integral = 0.0;
                for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                    integral += _entropy->EntropyPrimeDual( dot( alphaRed, _reducedMomentBasis[idx_quad] ) ) * _weights[idx_quad];
                }
                alphaCorr[0] =  -log( integral );    // log trafo
                for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                    uCorr += moments[idx_quad]*_entropy->EntropyPrimeDual( dot( alphaCorr, moments[idx_quad] ) ) * _weights[idx_quad];
                }
                //Rotate uCorr onto x-axis
                Vector u1{ uCorr[1], uCorr[2] };
                Vector u1_star(2,0.0);

                Matrix rotationMatsCorr  = CreateRotator( u1 );

                u1_star = RotateM1( u1, rotationMatsCorr );    //, _rotationMatsT[idx_cell] );
                if (u1_star[1]> 0.000001){
                    ErrorMessages::Error("Somethings wrong" + std::to_string(u1[1]) + " " + std::to_string(u1_star[1]),CURRENT_FUNCTION);
                }
                alphaRed = RotateM1( alphaRed, rotationMatsCorr ); // rotate corresponding alpha
                }
                // Rotate Back
                alphaRed = RotateM1( alphaRed, _rotationMatsT[idx_cell] );
                alpha_norms[idx_cell] = norm( alphaRed ) * norm( alphaRed );

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

                alpha_norms[idx_cell] = norm( alphaRed ) * norm( alphaRed );


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

void NeuralNetworkOptimizer::InferenceSphericalHarmonics2D( VectorVector& alpha, VectorVector& u,const VectorVector& moments, Vector& alpha_norms ) {
    unsigned servingSize = _settings->GetNCells();

    Matrix rot180 = CreateRotatorSphericalHarmonics2D( M_PI, -1.0,0.0 );

    if( _settings->GetEnforceNeuralRotationalSymmetry() ) {    // Rotation Preprocessing
#pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
            _rotationMats[idx_cell]  = CreateRotatorSphericalHarmonics2D( 0.0, u[idx_cell][1], u[idx_cell][2] );
            _rotationMatsT[idx_cell] = blaze::trans( _rotationMats[idx_cell] );

            u[idx_cell] = _rotationMats[idx_cell] * u[idx_cell];
            // Save rotated moment
            for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                _modelServingVectorU[idx_cell * ( _nSystem - 1 ) + idx_sys] = (float)( u[idx_cell][idx_sys + 1] );
            }

            // Rotate Moment by 180 degrees and save mirrored moment
            u[idx_cell] = rot180 * u[idx_cell];

            for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                _modelServingVectorU[( _settings->GetNCells() + idx_cell ) + idx_sys] = (float)( u[idx_cell][idx_sys + 1] );
            }
        }
        servingSize *= 2;
    }
    else {    // No Preprocessing
#pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
            u[idx_cell][2] = 0.0;  //manually enforce slab geometry
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
#pragma omp parallel for
            for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
                Vector alphaRed       = Vector( _nSystem, 0.0 );    // local reduced alpha (with dummy entry at 0)
                Vector alphaRedMirror = Vector( _nSystem, 0.0 );    // local reduced mirrored alpha (with dummy entry at 0)
                for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                    alphaRed[idx_sys + 1] = (double)_modelServingVectorAlpha[idx_cell * ( _nSystem - 1 ) + idx_sys];
                    alphaRedMirror[idx_sys + 1] =
                        (double)_modelServingVectorAlpha[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + idx_sys];
                }
                // Mirror back
                alphaRedMirror = rot180 * alphaRedMirror;
                // Average
                alphaRed = ( alphaRed + alphaRedMirror ) / 2.0;
                // Rotate Back
                alphaRed = _rotationMatsT[idx_cell] * alphaRed;
                // Compute norms for scaling
                alpha_norms[idx_cell] = norm( alphaRed ) * norm( alphaRed );

                // Restore alpha_0
                double integral = 0.0;
                for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                    integral += _entropy->EntropyPrimeDual( dot( alphaRed, _reducedMomentBasis[idx_quad] ) ) * _weights[idx_quad];
                }
                alpha[idx_cell][0] = -(log( integral )+ log(moments[0][0]) )/moments[0][0];     // log trafo

                // Store output
                for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                    alpha[idx_cell][idx_sys + 1] = alphaRed[idx_sys];
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
            alpha[idx_cell][0] = -(log( integral )+ log(moments[0][0]) )/moments[0][0];     // log trafo
        }
    }
    _modelServingVectorAlpha.clear();
}

void NeuralNetworkOptimizer::InferenceSphericalHarmonics( VectorVector& alpha, VectorVector& u,const VectorVector& moments, Vector& alpha_norms ) {
    unsigned servingSize = _settings->GetNCells();

    Matrix rot180 = CreateRotatorSphericalHarmonics( M_PI, -1.0,0.0 );

    if( _settings->GetEnforceNeuralRotationalSymmetry() ) {    // Rotation Preprocessing
//#pragma omp parallel for
            for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {

                u[idx_cell][2] = 0.0;  //manually enforce slab geometry

                double theta = atan2( u[idx_cell][1], u[idx_cell][3]); // Calculate angle using atan2 function

                _rotationMats[idx_cell]  = CreateRotatorSphericalHarmonics( theta, u[idx_cell][1], u[idx_cell][3] );
                _rotationMatsT[idx_cell] = blaze::trans( _rotationMats[idx_cell] );
                Vector Ru =_rotationMats[idx_cell] * u[idx_cell];
                Vector RRu =  _rotationMatsT[idx_cell]*Ru;
                if (norm(RRu-u[idx_cell])>1e-6){
                    TextProcessingToolbox::PrintMatrix(  _rotationMats[idx_cell]);
                    TextProcessingToolbox::PrintMatrix(  _rotationMatsT[idx_cell]);
                    std::cout << u[idx_cell]<<" | "<< Ru << " | "<<RRu <<" " << norm(RRu-u[idx_cell]) <<"\n";
                    ErrorMessages::Error("Rotations not proper" ,CURRENT_FUNCTION);
                }

                //std::cout << u[idx_cell]<<" | "<< a  <<" | "<< b <<" | "<< r <<" | " << theta<< "\n";
                //exit(1);
                //std::cout << u[idx_cell]<<" \n| "<< _rotationMats[idx_cell] * u[idx_cell]  <<"\n";

                //if (normu1_slab > 1e-15){
                //    std::cout << "ux,uy " << u[idx_cell][1]<< u[idx_cell][3] << "\n ";

                //    std::cout << "theta " << theta << "\n ";

                //    std::cout << "u orig: " << u[idx_cell]<<" \n| u_rot: "<< Ru  <<"\n";
                //    Vector t  = _rotationMatsT[idx_cell]* Ru;
                //    std::cout << "u backrot: " << t << "\n --------------\n";
                //}

                //if (Ru[3]/u[idx_cell][3] >1e-10){
                //    std::cout << "ux,uy " << u[idx_cell][1]<< u[idx_cell][3] << "\n ";

                //    std::cout << "theta " << theta << "\n ";

                //    std::cout << "u orig: " << u[idx_cell]<<" \n| u_rot: "<< Ru  <<"\n";
                //    Vector t  = _rotationMatsT[idx_cell]* Ru;
                //    std::cout << "u backrot: " << t << "\n --------------\n";

                //    std::cout << "alert \n --------------\n";
                //    std::cout << _rotationMats[idx_cell] << "\n"<<_rotationMatsT[idx_cell] ;
                //    exit(1);
                //}
                Vector Ru_minus = rot180 * Ru;

                for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                     _modelServingVectorU[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 )  + idx_sys] =  (float)( Ru[idx_sys + 1] );
                     _modelServingVectorU[idx_cell * ( _nSystem - 1 ) + idx_sys] = (float)(Ru_minus[idx_sys + 1] );
                }
            }
        servingSize *= 2;
    }
    else {    // No Preprocessing
#pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
            u[idx_cell][2] = 0.0;  //manually enforce slab geometry
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
#pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
            Vector alphaRed = Vector( _nSystem - 1, 0.0 );    // local reduced alpha
            Vector alphaP = Vector( _nSystem , 0.0 );    // local  alpha
            Vector alpha_P_Mirror = Vector( _nSystem , 0.0 );    // local mirrored alpha


            for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                alphaP[idx_sys + 1]         = (double)_modelServingVectorAlpha[idx_cell * ( _nSystem - 1 ) + idx_sys];
                alpha_P_Mirror[idx_sys+1] = (double)_modelServingVectorAlpha[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + idx_sys];
            }
            alpha_P_Mirror = rot180*alpha_P_Mirror;
            alpha[idx_cell]= ( alphaP + alpha_P_Mirror ) / 2;    // average (and store in alpha)


            alpha[idx_cell] = _rotationMatsT[idx_cell] * alpha[idx_cell]; // Rotate back
            //alpha[idx_cell][2] = 0.0;    //manually enforce slab geometry
            // compute norm squared for dynamic ansatz  alpha[idx_cell][0] = 0 by construction
            alpha_norms[idx_cell] = norm(  alpha[idx_cell] ) * norm(  alpha[idx_cell] );


            // Restore alpha_0
            double integral = 0.0;
            for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                integral += _entropy->EntropyPrimeDual( dot( alphaRed, _reducedMomentBasis[idx_quad] ) ) * _weights[idx_quad];
            }
            alpha[idx_cell][0] = -(log( integral )+ log(moments[0][0]) )/moments[0][0];     // log trafo
        }


        //for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {

        //        Vector alpha_P       = Vector( _nSystem , 0.0 );    // local reduced alpha
        //        Vector alpha_P_Mirror = Vector( _nSystem , 0.0 );    // local reduced mirrored alpha
        //        Vector alpha_P_Red       = Vector( _nSystem-1 , 0.0 );    // local reduced alpha


        //        for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) { // Load tf vectors
        //            alpha_P[idx_sys+1]       = (double)_modelServingVectorAlpha[idx_cell * ( _nSystem - 1 ) + idx_sys];
        //            alpha_P_Mirror[idx_sys+1] = (double)_modelServingVectorAlpha[( _settings->GetNCells() + idx_cell ) * ( _nSystem - 1 ) + idx_sys];
        //        }
        //        //alpha_P[2]=0.0;
        //        //alpha_P_Mirror[2]=0.0;

        //        // Rotate back mirrored alpha
        //        //alpha_P_Mirror = rot180*alpha_P_Mirror;
        //        alpha[idx_cell]= alpha_P;//( alpha_P + alpha_P_Mirror ) / 2;    // average (and store in alpha)
        //        alpha_norms[idx_cell] = 0.0;// compute norm squared for dynamic ansatz.
        //        // Rotated back to position of original u
        //        alpha_P = _rotationMatsT[idx_cell] * alpha_P;
        //        for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) { // Load tf vectors
        //            alpha_P_Red[idx_sys]   = alpha_P[1 + idx_sys];
        //            alpha_norms[idx_cell] += alpha_P_Red[idx_sys]* alpha_P_Red[idx_sys];
        //        }

        //        // Restore alpha_0
        //        double integral = 0.0;
        //        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        //            integral += _entropy->EntropyPrimeDual( dot( alpha_P_Red, _reducedMomentBasis[idx_quad] ) ) * _weights[idx_quad];
        //        }
        //        alpha[idx_cell][0] = -( log( integral ) + log( moments[0][0] ) ) / moments[0][0];    // log trafo
        //    }
   }

    else {
#pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _settings->GetNCells(); idx_cell++ ) {
            Vector alphaRed = Vector( _nSystem - 1, 0.0 );    // local reduced alpha
            for( unsigned idx_sys = 0; idx_sys < _nSystem - 1; idx_sys++ ) {
                alphaRed[idx_sys]            = (double)_modelServingVectorAlpha[idx_cell * ( _nSystem - 1 ) + idx_sys];
                alpha[idx_cell][idx_sys + 1] = alphaRed[idx_sys];
            }
            alphaRed[1] = 0.0; //manually enforce slab geometry
            alpha[idx_cell][2] = 0.0;    //manually enforce slab geometry
            alpha_norms[idx_cell] = norm( alphaRed ) * norm( alphaRed ); // compute norm squared for dynamic ansatz.

            // Restore alpha_0
            double integral = 0.0;
            for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                integral += _entropy->EntropyPrimeDual( dot( alphaRed, _reducedMomentBasis[idx_quad] ) ) * _weights[idx_quad];
            }
            alpha[idx_cell][0] = -(log( integral )+ log(moments[0][0]) )/moments[0][0];     // log trafo
            Vector u_recons = Vector(_nSystem, 0.0);

            for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                u_recons += _entropy->EntropyPrimeDual( dot( alpha[idx_cell], moments[idx_quad] ) ) * _weights[idx_quad]* moments[idx_quad];
            }
            //if (norm(u_recons-u[idx_cell])/norm(u[idx_cell])>1e-2){
            //    std::cout << u[idx_cell]<<" \n\n "<< u_recons << "\n\n >> " << norm(u_recons-u[idx_cell]) <<"\n";
            //    ErrorMessages::Error("Reconstruction failed" ,CURRENT_FUNCTION);
            //}
        }
    }
    _modelServingVectorAlpha.clear();
}

void NeuralNetworkOptimizer::SolveMultiCell( VectorVector& alpha, VectorVector& u, const VectorVector& moments, Vector& alpha_norms ) {
    if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
        if (_settings->GetDim() == 2){
            InferenceSphericalHarmonics2D( alpha, u,moments, alpha_norms );
        }
        else{
            InferenceSphericalHarmonics( alpha, u,moments, alpha_norms );
        }
    }
    else {
        InferenceMonomial( alpha, u,moments, alpha_norms );
    }
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
