#include "fluxes/upwindflux.h"

#include "blaze/math/expressions/DMatDVecMultExpr.h"      // for DVecScalarMul...
#include "blaze/math/expressions/DMatScalarMultExpr.h"    // for operator*
#include "blaze/math/expressions/DVecDVecAddExpr.h"       // for DVecDVecAddExpr
#include "blaze/math/expressions/DVecDVecSubExpr.h"       // for DVecDVecSubExpr
#include "blaze/math/expressions/DVecScalarMultExpr.h"    // for operator*
#include "blaze/math/expressions/DenseMatrix.h"           // for DenseMatrix
#include "blaze/math/expressions/DenseVector.h"           // for DenseVector
#include "blaze/math/expressions/MatScalarMultExpr.h"     // for operator*
#include "blaze/math/expressions/Vector.h"                // for isSame
#include "blaze/math/simd/Add.h"                          // for operator+
#include "blaze/math/simd/BasicTypes.h"                   // for operator+=
#include "blaze/math/simd/Mult.h"                         // for operator*
#include "blaze/math/simd/Set.h"                          // for set
#include "blaze/math/simd/Sub.h"                          // for operator-
#include "blaze/math/simd/Sum.h"                          // for sum
#include "blaze/math/smp/default/DenseVector.h"           // for smpAddAssign
#include "fluxes/numericalflux.h"                         // for NumericalFlux
#include <emmintrin.h>                                    // for _mm_mul_pd

UpwindFlux::UpwindFlux() : NumericalFlux() {}

double UpwindFlux::Flux( const Vector& Omega, double psiL, double psiR, const Vector& n ) const {
    double inner = Omega[0] * n[0] + Omega[1] * n[1];
    if( inner > 0 ) {
        return inner * psiL;
    }
    else {
        return inner * psiR;
    }
}

Vector UpwindFlux::Flux( const Matrix AxPlus,
                         const Matrix AxMinus,
                         const Matrix AyPlus,
                         const Matrix AyMinus,
                         const Matrix /*AzPlus*/,
                         const Matrix /*AzMinus*/,
                         const Vector psiL,
                         const Vector psiR,
                         const Vector n ) const {
    // 2d only atm!!!

    Vector resultFlux( psiR.size(), 0 );
    // x dir
    if( n[0] > 0 ) {
        resultFlux += n[0] * AxPlus * psiL + n[0] * AxMinus * psiR;
    }
    else {
        resultFlux += n[0] * AxPlus * psiR + n[0] * AxMinus * psiL;
    }
    // y dir
    if( n[1] > 0 ) {
        resultFlux += n[1] * AyPlus * psiL + n[1] * AyMinus * psiR;
    }
    else {
        resultFlux += n[1] * AyPlus * psiR + n[1] * AyMinus * psiL;
    }

    return resultFlux;
}

void UpwindFlux::FluxVanLeer( const Matrix& Ax,
                              const Matrix& AxAbs,
                              const Matrix& /*Ay*/,
                              const Matrix& /*AyAbs*/,
                              const Matrix& Az,
                              const Matrix& AzAbs,
                              const Vector& psiL,
                              const Vector& psiR,
                              const Vector& n,
                              Vector& resultFlux ) const {

    // resultFlux += n[0] * ( AxPlus * psiL + AxMinus * psiR ) + n[1] * ( AyPlus * psiL + AyMinus * psiR );
    resultFlux += 0.5 * ( n[0] * ( Ax * ( psiL + psiR ) - AxAbs * ( psiR - psiL ) ) + n[1] * ( Az * ( psiL + psiR ) - AzAbs * ( psiR - psiL ) ) );
}
