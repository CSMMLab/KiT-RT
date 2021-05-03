#ifndef ELECTRONRT_H
#define ELECTRONRT_H

#include "epics.h"
#include "problembase.h"

class ElectronRT : public ProblemBase
{
  private:
    /**
     * @brief LoadXSH2O loads (total) scattering cross sections for water from file and saves them in _xsH2O and _totalxsH2O
     * @param fileSigmaS is name of scattering cross section file
     * @param fileSigmaT is name of total cross section file
     */
    void LoadXSH20( std::string fileSigmaS, std::string fileSigmaT );

    ElectronRT() = delete;

  public:
    ElectronRT( Config* settings, Mesh* mesh );
    virtual ~ElectronRT();

    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
    std::vector<Matrix> GetScatteringXSE( const Vector& energies, const Matrix& angles ) override;
    /**
     * @brief GetTotalXSE gives back vector of total cross sections for
     *        energies in vector energy
     * @param energies is the energy the cross section is queried for
     */
    Vector GetTotalXSE( const Vector& energies ) override;
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    VectorVector SetupIC() override;
    std::vector<double> GetDensity( const VectorVector& cellMidPoints ) override;
};

#endif    // ELECTRONRT_H
