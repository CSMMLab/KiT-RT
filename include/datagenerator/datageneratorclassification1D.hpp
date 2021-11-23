#ifndef DATAGENERATORCLASSIFICATION1D_H
#define DATAGENERATORCLASSIFICATION1D_H

#include "datageneratorclassification.hpp"

class DataGeneratorClassification1D : public DataGeneratorClassification
{
  public:
    /*! @brief Class constructor. Generates training data for neural network approaches using
     *          spherical harmonics and an entropy functional and the quadrature specified by
     *          the options file.
     *   @param settings config class with global information*/
    DataGeneratorClassification1D( Config* settings );
    ~DataGeneratorClassification1D();


  protected:
    void ComputeMoments() final;        /*!< @brief Pre-Compute Moments at all quadrature points. */
    void SampleMultiplierAlpha() final; /*!< @brief Sample Lagrange multipliers alpha, with mean values corresponding to a maxwellian distribution */
    void PrintTrainingData() final;     /*!< @brief : Print computed training data to csv file and screen */
};

#endif    // DATAGENERATORCLASSIFICATION1D_H
