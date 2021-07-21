#ifndef DATAGENERATORCLASSIFICATION_H
#define DATAGENERATORCLASSIFICATION_H

#include "datageneratorbase.h"

class DataGeneratorClassification : public DataGeneratorBase
{
  public:
    /*! @brief Class constructor. Generates training data for neural network approaches using
     *          spherical harmonics and an entropy functional and the quadrature specified by
     *          the options file.
     *   @param settings config class with global information*/
    DataGeneratorClassification( Config* settings );
    ~DataGeneratorClassification();

    void ComputeTrainingData() override;

  protected:
    Vector _pdfClassification; /*!< @brief One-hot vector with classification, if kinetic pdf is within or outside KL Divergence threshold*/
    Vector _maxwellian;        /*!< @brief Maxwellian pdf evaluated at the quadrature points */
    double _leftBound;
    double _rightBound;

    // IO routines
    void PrintTrainingData() override; /*!< @brief : Print computed training data to csv file and screen */

    // Helper functions
    void ComputeMoments() override; /*!< @brief Pre-Compute Moments at all quadrature points. */
    void ClassifyDensity();         /*!< @brief Checks, if the pdf of each Lagrange multiplier is within the KL distance of the maxwellian */
    /*!< @brief Computes the Kullback Leibler Divergence of the pdfs f1 and f2, both pfds are evaluated at their quadrature points
                  @param: f1,f2. Evaluation of the pdf at their quadrature points. length of vector must be _nq.
             */
    double ComputeKLDivergence( Vector& f1, Vector& f2 );
    /*!< @brief Evalutes the maxwellian at the quadrature points of _quadrature.
         @param rho: density
         @param u: bulk velocity
         @param T: Temperature*/
    Vector ComputeMaxwellian( double rho, double u, double T );
};

#endif    // DATAGENERATORCLASSIFICATION_H
