#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>

typedef blaze::DynamicMatrix<double> Matrix;
typedef std::vector<blaze::DynamicVector<double>> VecVec;
typedef std::vector<blaze::DynamicVector<unsigned>> VecVecU;
typedef blaze::DynamicVector<double> Vector;
typedef blaze::DynamicVector<unsigned> VectorU;

#endif    // TYPEDEFS_H
