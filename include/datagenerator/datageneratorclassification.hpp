#ifndef DATAGENERATORCLASSIFICATION_H
#define DATAGENERATORCLASSIFICATION_H

#include "datageneratorbase.hpp"

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
    Vector _pdfClassification;    /*!< @brief One-hot vector with classification, if kinetic pdf is within or outside KL Divergence threshold*/
    Vector _maxwellian;           /*!< @brief Maxwellian pdf evaluated at the quadrature points */
    double _maxVelocity;          /*!< @brief Bound of the velocity domain (1D) */
    VectorVector _kineticDensity; /*!< @brief vector if sampled kinetic densities, evaluated at quadrature points */

    // IO routines
    void PrintTrainingData() override; /*!< @brief : Print computed training data to csv file and screen */

    // Helper functions
    virtual void ComputeMoments() override = 0; /*!< @brief Pre-Compute Moments at all quadrature points. */
    void ClassifyDensity(); /*!< @brief Checks, if the pdf of each Lagrange multiplier is within the KL distance of the maxwellian */
    /*!< @brief Computes the Kullback Leibler Divergence of the pdfs f1 and f2, both pfds are evaluated at their quadrature points
                  @param: f1,f2. Evaluation of the pdf at their quadrature points. length of vector must be _nq.
             */
    double ComputeKLDivergence( Vector& f1, Vector& f2 );
    /*!< @brief Evalutes the maxwellian at the quadrature points of _quadrature.
         @param rho: density
         @param u: bulk velocity
         @param T: Temperature*/
    Vector ComputeMaxwellian( double rho, double u, double T );

    /*!< @brief Computes the kinetic density from the given Lagrange multipliers alpha at the quadrature points.
     *          f = exp(alpha*m) and stores it in _kineticDensity             */
    void ReconstructKineticDensity();

    /*!< @brief Sample Lagrange multipliers alpha, with mean values corresponding to a maxwellian distribution */
    virtual void SampleMultiplierAlpha() override = 0;
};

#endif    // DATAGENERATORCLASSIFICATION_H
