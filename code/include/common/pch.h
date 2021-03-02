#ifndef PCH_H
#define PCH_H

#include <algorithm>                 // for find
#include <assert.h>                  // for assert
#include <bits/exception.h>          // for exception
#include <bits/types/struct_tm.h>    // for tm
#include <chrono>                    // for seconds
#include <cmath>                     // for sqrt, acos, cos
#include <emmintrin.h>               // for _mm_mul_pd
#include <ext/alloc_traits.h>        // for __alloc_trait...
#include <filesystem>                // for path
#include <fstream>
#include <iostream>    // for ifstream, basic_istream
#include <iterator>    // for begin, end
#include <limits>      // for numeric_limits
#include <map>         // for map
#include <math.h>      // for sqrt, acos, cos
#include <memory>      // for allocator_traits<>::value_type
#include <mpi.h>
#include <omp.h>
#include <stdlib.h>    // for abs
#include <string.h>    // for strncpy
#include <string>      // for string
#include <time.h>      // for localtime, strftime, time
#include <utility>     // for pair
#include <vector>      // for vector

#include "blaze/math/CompressedMatrix.h"
#include "blaze/math/DynamicMatrix.h"
#include "blaze/math/DynamicVector.h"
#include "blaze/math/Vector.h"                            // for dot
#include "blaze/math/dense/DenseVector.h"                 // for operator*=
#include "blaze/math/expressions/DVecDVecAddExpr.h"       // for DVecDVecAddExpr
#include "blaze/math/expressions/DVecDVecInnerExpr.h"     // for operator*
#include "blaze/math/expressions/DVecDVecSubExpr.h"       // for DVecDVecSubExpr
#include "blaze/math/expressions/DVecNormExpr.h"          // for l2Norm, norm
#include "blaze/math/expressions/DVecReduceExpr.h"        // for min, sum
#include "blaze/math/expressions/DVecScalarDivExpr.h"     // for operator/
#include "blaze/math/expressions/DVecScalarMultExpr.h"    // for DVecScalarMul...
#include "blaze/math/expressions/DVecTransExpr.h"         // for trans
#include "blaze/math/expressions/DenseVector.h"           // for DenseVector
#include "blaze/math/expressions/SVecSVecInnerExpr.h"     // for operator*
#include "blaze/math/expressions/SVecTransExpr.h"         // for trans
#include "blaze/math/expressions/Vector.h"                // for derestrict
#include "blaze/math/simd/Add.h"                          // for operator+
#include "blaze/math/simd/BasicTypes.h"                   // for operator+=
#include "blaze/math/simd/Mult.h"                         // for operator*
#include "blaze/math/simd/Reduce.h"                       // for reduce
#include "blaze/math/simd/Set.h"                          // for set
#include "blaze/math/simd/Sub.h"                          // for operator-
#include "blaze/math/simd/Sum.h"                          // for sum
#include "blaze/math/smp/default/DenseVector.h"           // for smpAssign
#include "blaze/math/views/Row.h"                         // for row
#include "blaze/math/views/row/Sparse.h"                  // for Row
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_sinks.h"
#include "spdlog/spdlog.h"

#include "common/config.h"
#include "common/globalconstants.h"
#include "common/mesh.h"
#include "common/optionstructure.h"
#include "common/typedef.h"
#include "entropies/entropybase.h"
#include "fluxes/numericalflux.h"
#include "io.h"
#include "kernels/scatteringkernelbase.h"
#include "optimizers/optimizerbase.h"
#include "problems/problembase.h"
#include "quadratures/quadraturebase.h"
#include "solvers/solverbase.h"
#include "toolboxes/errormessages.h"

#endif    // PCH_H
