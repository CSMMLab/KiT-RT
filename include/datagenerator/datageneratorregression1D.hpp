/*!
 * \file datageneratorregression1D.h
 * \brief Class to generate data for the neural entropy closure in 1D spatial dimensions
 * \author S. Schotthoefer
 */

#ifndef DATAGENERATORREGRESSION1D_H
#define DATAGENERATORREGRESSION1D_H

#include "datageneratorregression.hpp"

class DataGeneratorRegression1D : public DataGeneratorRegression
{
  public:
    /*! @brief Class constructor. Generates training data for neural network approaches using
     *          spherical harmonics and an entropy functional and the quadrature specified by
     *          the options file.
     *   @param settings config class with global information*/
    DataGeneratorRegression1D( Config* settings );
    ~DataGeneratorRegression1D();

  private:
    // Main methods
    void SampleSolutionU() override; /*!< @brief Samples solution vectors u */

    // Helper functions
    void ComputeMoments() override;        /*!< @brief Pre-Compute Moments at all quadrature points. */
    void ComputeSetSizeU() override;       /*!< @brief Computes the size of the training set, depending on the chosen settings.*/
    void CheckRealizability() override;    // Debugging helper
};
#endif    // DATAGENERATORREGRESSION1D_H
