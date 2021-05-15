/*!
 * \file datageneratorbase.h
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

class DataGeneratorBase
{
  public:
    /*! @brief Class constructor. Generates training data for neural network approaches using
     *          spherical harmonics and an entropy functional and the quadrature specified by
     *          the options file.
     *   @param settings config class with global information*/
    DataGeneratorBase( Config* settings );
    virtual ~DataGeneratorBase();

    /*! @brief Create a datagenerator (1D or 3D)
     *  @param settings Pointer to the config file
     *  @returns Pointer to the createt basis class */
    static DataGeneratorBase* Create( Config* settings );

    /*! @brief computes the training data set.
     *          Realizable set is sampled uniformly.
     *          Prototype: 1D, u in [0,100] */
    void ComputeTrainingData();

    inline VectorVector GetuSol() { return _uSol; }                /*! @brief Get the computed solution vector uSol */
    inline VectorVector GetAlpha() { return _alpha; }              /*! @brief Get the computed vector alpha */
    inline std::vector<double> GethEntropy() { return _hEntropy; } /*! @brief Get the computed entropy value h */

  protected:
    Config* _settings; /*!< @brief config class for global information */

    VectorVector _uSol;            /*!< @brief vector with moments. Size: (setSize,basisSize)*/
    VectorVector _alpha;           /*!< @brief vector with Lagrange multipliers. Size: (setSize,basisSize)*/
    std::vector<double> _hEntropy; /*!< @brief vector with entropy values. Size: (setSize) */

    unsigned long _setSize; /*!< @brief Size of the whole training Set */
    unsigned long
        _gridSize; /*!< @brief Size of the grid discretizing moment U0 for higher order sampling (has different uses for different samplers)*/

    unsigned short _maxPolyDegree; /*!< @brief Max Order of Spherical Harmonics */
    unsigned _nTotalEntries;       /*!< @brief Total number of equations in the system */

    QuadratureBase* _quadrature;    /*!< @brief quadrature to create members below */
    unsigned _nq;                   /*!< @brief number of quadrature points */
    VectorVector _quadPoints;       /*!<  @brief quadrature points, dim(_quadPoints) = (_nq,spatialDim) */
    Vector _weights;                /*!<  @brief quadrature weights, dim(_weights) = (_nq) */
    VectorVector _quadPointsSphere; /*!<  @brief (my,phi), dim(_quadPoints) = (_nq,2) */

    SphericalBase* _basisGenerator; /*!< @brief Class to compute and store current spherical harmonics basis */
    VectorVector _momentBasis; /*!< @brief Moment Vector pre-computed at each quadrature point: dim= _nq x _nTotalEntries */

    NewtonOptimizer* _optimizer; /*!< @brief Class to solve minimal entropy problem */
    EntropyBase* _entropy;       /*!< @brief Class to handle entropy functional evaluations */

    // Main methods
    virtual void SampleSolutionU() = 0; /*!< @brief Samples solution vectors u */
    void SampleMultiplierAlpha();       /*!< @brief Sample Lagrange multipliers alpha */
    void ComputeEntropyH_dual();        /*!< @brief Compute the entropy functional at (u,alpha) in dual formulation */
    void ComputeEntropyH_primal();      /*!< @brief  Compute the entropy functional at (u,alpha) in primal formulation */
    void ComputeRealizableSolution();   /*!< @brief make u the realizable moment to alpha, since Newton has roundoff errors. */

    // IO routines
    void PrintTrainingData(); /*!< @brief : Print computed training data to csv file and screen */
    void PrintLoadScreen();   /*!< @brief Print screen IO*/

    // Helper functions
    virtual void ComputeMoments()     = 0; /*!< @brief Pre-Compute Moments at all quadrature points. */
    virtual void CheckRealizability() = 0; /*!< @brief Debugging helper. Will be removed */
    virtual void ComputeSetSizeU()    = 0; /*!< @brief Computes the size of the training set, depending on the chosen settings.*/
    void ComputeSetSizeAlpha();            /*!< @brief Computes the seSize for alphasampling */
    void AdaptBasisSize();                 /*!< @brief In case of normal sampling, deletes zero order basis for dimension reduction. */
};
#endif    // DATAGENERATOR_H
