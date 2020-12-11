/*!
 * \file datagenerator.h
 * \brief Class to generate data for the neural entropy closure
 * \author S. Schotthoefer
 */

#ifndef DATAGENERATOR_H
#define DATAGENERATOR_H

#include "common/typedef.h"
#include <vector>

class SphericalBase;
class QuadratureBase;
class Config;
class NewtonOptimizer;
class EntropyBase;

class nnDataGenerator
{
  public:
    /*! @brief: Class constructor. Generates training data for neural network approaches using
     *          spherical harmonics and an entropy functional and the quadrature specified by
     *          the options file.
     *   @param: setSize: number of elements in training set
     *           basisSize: length of spherical harmonics basis (maybe redundant)*/
    nnDataGenerator( Config* settings );
    ~nnDataGenerator() {}

    /*! @brief: computes the training data set.
     *          Realizable set is sampled uniformly.
     *          Prototype: 1D, u\in[0,100] */
    void computeTrainingData();

    /*!   @brief: Writes the training data to file
     *            Filename encryption: [TODO] */
    void writeTrainingDataToCSV();

  private:
    Config* _settings; /*! @brief config class for global information */

    VectorVector _uSol;            /*! @brief: vector with moments. Size: (setSize,basisSize)*/
    VectorVector _alpha;           /*! @brief: vector with Lagrange multipliers. Size: (setSize,basisSize)*/
    std::vector<double> _hEntropy; /*! @brief: vector with entropy values. Size: (setSize) */

    unsigned long _setSize;
    unsigned short _LMaxDegree; /*! @brief: Max Order of Spherical Harmonics */
    unsigned _nTotalEntries;    /*! @brief: Total number of equations in the system */

    QuadratureBase* _quadrature;    /*! @brief quadrature to create members below */
    unsigned _nq;                   /*! @brief number of quadrature points */
    VectorVector _quadPoints;       /*!  @brief quadrature points, dim(_quadPoints) = (_nq,spatialDim) */
    Vector _weights;                /*!  @brief quadrature weights, dim(_weights) = (_nq) */
    VectorVector _quadPointsSphere; /*!  @brief (my,phi), dim(_quadPoints) = (_nq,2) */

    SphericalBase* _basis; /*! @brief: Class to compute and store current spherical harmonics basis */
    VectorVector _moments; /*! @brief: Moment Vector pre-computed at each quadrature point: dim= _nq x _nTotalEntries */

    NewtonOptimizer* _optimizer; /*! @brief: Class to solve minimal entropy problem */
    EntropyBase* _entropy;       /*! @brief: Class to handle entropy functional evaluations */

    // Main methods
    void SampleSolutionU();        /*! @brief: Samples solution vectors u */
    void ComputeEntropyH_dual();   /*! @brief: Compute the entropy functional at (u,alpha) in dual formulation */
    void ComputeEntropyH_primal(); /*! @brief:  Compute the entropy functional at (u,alpha) in primal formulation */

    // IO routines
    void PrintTrainingData(); /*! @brief : Print computed training data to csv file and screen */

    // Helper functions
    void ComputeMoments(); /*! @brief: Pre-Compute Moments at all quadrature points. */
};
#endif    // DATAGENERATOR_H
