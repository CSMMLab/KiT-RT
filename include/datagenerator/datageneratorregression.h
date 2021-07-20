#ifndef DATAGENERATORREGRESSION_H
#define DATAGENERATORREGRESSION_H

#include "datageneratorbase.h"

class DataGeneratorRegression : public DataGeneratorBase
{
  public:
    /*! @brief Class constructor. Generates training data for neural network approaches using
     *          spherical harmonics and an entropy functional and the quadrature specified by
     *          the options file.
     *   @param settings config class with global information*/
    DataGeneratorRegression( Config* settings );
    virtual ~DataGeneratorRegression();

    void ComputeTrainingData() override;

    inline VectorVector GetuSol() { return _uSol; }                /*! @brief Get the computed solution vector uSol */
    inline VectorVector GetAlpha() { return _alpha; }              /*! @brief Get the computed vector alpha */
    inline std::vector<double> GethEntropy() { return _hEntropy; } /*! @brief Get the computed entropy value h */

  protected:
    VectorVector _uSol;            /*!< @brief vector with moments. Size: (setSize,basisSize)*/
    std::vector<double> _hEntropy; /*!< @brief vector with entropy values. Size: (setSize) */
    unsigned long
        _gridSize; /*!< @brief Size of the grid discretizing moment U0 for higher order sampling (has different uses for different samplers)*/

    // Main methods
    virtual void SampleSolutionU() = 0; /*!< @brief Samples solution vectors u */
    void ComputeEntropyH_dual();        /*!< @brief Compute the entropy functional at (u,alpha) in dual formulation */
    void ComputeEntropyH_primal();      /*!< @brief  Compute the entropy functional at (u,alpha) in primal formulation */
    void ComputeRealizableSolution();   /*!< @brief make u the realizable moment to alpha, since Newton has roundoff errors. */

    // IO routines
    void PrintTrainingData() override; /*!< @brief : Print computed training data to csv file and screen */

    // Helper functions
    virtual void ComputeMoments()     = 0; /*!< @brief Pre-Compute Moments at all quadrature points. */
    virtual void CheckRealizability() = 0; /*!< @brief Debugging helper. Will be removed */
    virtual void ComputeSetSizeU()    = 0; /*!< @brief Computes the size of the training set, depending on the chosen settings.*/
    void ComputeSetSizeAlpha();            /*!< @brief Computes the seSize for alphasampling */
};

#endif    // DATAGENERATORREGRESSION_H
