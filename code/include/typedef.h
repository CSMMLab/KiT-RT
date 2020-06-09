#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>

typedef blaze::DynamicMatrix<double> Matrix;    // Row Major by default!
typedef blaze::DynamicMatrix<double, blaze::columnMajor> MatrixCol;
typedef blaze::SymmetricMatrix<blaze::DynamicMatrix<double, blaze::rowMajor>> SymMatrix;    // Row Major by default!
typedef blaze::SymmetricMatrix<blaze::DynamicMatrix<double, blaze::columnMajor>> SymMatrixCol;
typedef blaze::DiagonalMatrix<blaze::DynamicMatrix<double, blaze::rowMajor>> DiagMatrix;
typedef blaze::CompressedMatrix<double> SparseMatrix;
typedef std::vector<blaze::DynamicVector<double>> VectorVector;
typedef std::vector<blaze::DynamicVector<unsigned>> VectorVectorU;
typedef blaze::DynamicVector<double> Vector;       // Column Vector by default!
typedef blaze::DynamicVector<unsigned> VectorU;    // Column Vector by default!

#endif    // TYPEDEFS_H
