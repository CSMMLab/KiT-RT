/*!
 * \file datagenerator1D.h
 * \brief Class to generate data for the neural entropy closure in 1D spatial dimensions
 * \author S. Schotthoefer
 */

#ifndef DATAGENERATOR1D_H
#define DATAGENERATOR1D_H

#include "toolboxes/datageneratorbase.h"

class DataGenerator1D : public DataGeneratorBase
{
  public:
    /*! @brief Class constructor. Generates training data for neural network approaches using
     *          spherical harmonics and an entropy functional and the quadrature specified by
     *          the options file.
     *   @param setSize: number of elements in training set
     *           basisSize: length of spherical harmonics basis (maybe redundant)*/
    DataGenerator1D( Config* settings );
    ~DataGenerator1D();

  private:
    // Main methods
    void SampleSolutionU() override; /*!< @brief Samples solution vectors u */

    // Helper functions
    void ComputeMoments() override; /*!< @brief Pre-Compute Moments at all quadrature points. */
    void ComputeSetSize() override; /*!< @brief Computes the size of the training set, depending on the chosen settings.*/
    void CheckRealizability() override;    // Debugging helper
};
#endif    // DATAGENERATOR1D_H
