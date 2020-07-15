#ifndef ELECTRONRT_H
#define ELECTRONRT_H

#include "physics.h"
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

    virtual VectorVector GetScatteringXS( const std::vector<double>& energies );
    virtual VectorVector GetTotalXS( const std::vector<double>& energies );
    virtual VectorVector GetScatteringXSE( const Vector& energies, const Vector& angles );
    virtual Vector GetTotalXSE( const Vector& energies );
    virtual std::vector<VectorVector> GetExternalSource( const std::vector<double>& energies );
    virtual std::vector<double> GetStoppingPower( const std::vector<double>& energies );
    virtual VectorVector SetupIC();
    std::vector<double> GetDensity( const VectorVector& cellMidPoints );
};

#endif    // ELECTRONRT_H
