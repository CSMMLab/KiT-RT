/*!
 * \file datageneratorregression3D.h
 * \brief Class to generate data for the neural entropy closure in 3D spatial dimensions
 * \author S. Schotthoefer
 */

#ifndef DATAGENERATORREGRESSION3D_H
#define DATAGENERATORREGRESSION3D_H

#include "datageneratorregression.hpp"

class DataGeneratorRegression3D : public DataGeneratorRegression
{
  public:
    /*! @brief Class constructor. Generates training data for neural network approaches using
     *          spherical harmonics and an entropy functional and the quadrature specified by
     *          the options file.
     *   @param settings config class with global information*/
    DataGeneratorRegression3D( Config* settings );
    ~DataGeneratorRegression3D();

  private:
    // Main methods
    void SampleSolutionU() override; /*!< @brief Samples solution vectors u */

    // Helper functions
    void ComputeMoments() override;        /*!< @brief Pre-Compute Moments at all quadrature points. */
    void ComputeSetSizeU() override;       /*!< @brief Computes the size of the training set, depending on the chosen settings.*/
    void CheckRealizability() override;    // Debugging helper
};
#endif    // DATAGENERATORREGRESSION3D_H
