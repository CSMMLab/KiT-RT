/*!
 * \file datagenerator3D.h
 * \brief Class to generate data for the neural entropy closure in 3D spatial dimensions
 * \author S. Schotthoefer
 */

#ifndef DATAGENERATOR3D_H
#define DATAGENERATOR3D_H

#include "datageneratorregression.h"

class DataGenerator3D : public DataGeneratorRegression
{
  public:
    /*! @brief Class constructor. Generates training data for neural network approaches using
     *          spherical harmonics and an entropy functional and the quadrature specified by
     *          the options file.
     *   @param settings config class with global information*/
    DataGenerator3D( Config* settings );
    ~DataGenerator3D();

  private:
    // Main methods
    void SampleSolutionU() override; /*!< @brief Samples solution vectors u */

    // Helper functions
    void ComputeMoments() override;        /*!< @brief Pre-Compute Moments at all quadrature points. */
    void ComputeSetSizeU() override;       /*!< @brief Computes the size of the training set, depending on the chosen settings.*/
    void CheckRealizability() override;    // Debugging helper
};
#endif    // DATAGENERATOR3D_H
