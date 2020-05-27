#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>

typedef blaze::DynamicMatrix<double> Matrix;
typedef blaze::CompressedMatrix<double> SparseMatrix;
typedef std::vector<blaze::DynamicVector<double>> VectorVector;
typedef std::vector<blaze::DynamicVector<unsigned>> VectorVectorU;
typedef blaze::DynamicVector<double> Vector;
typedef blaze::DynamicVector<unsigned> VectorU;

#endif    // TYPEDEFS_H
