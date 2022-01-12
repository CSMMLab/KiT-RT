#ifndef CSDPNSOLVER_JL_H
#define CSDPNSOLVER_JL_H

#include "solvers/pnsolver.hpp"

class SphericalBase;

class CSDPNSolver_JL : public PNSolver
{
  private:
    std::vector<double> _dose; /*!< @brief TODO */
    double _E_cutoff;
    Vector saveE_ref;

    Vector _eTrafo;
    Vector _sigmaTAtEnergy;

    // Physics acess
    Vector _angle; /*!< @brief angles */

    std::vector<Matrix> _sigmaSE; /*!<  @brief scattering cross section for all energies*/
    Vector _sigmaTE;              /*!<  @brief total cross section for all energies*/

    VectorVector _basisAtQuad; /*!<  @brief spherical harmonics basis at quadrature points*/

    // --- Private member variables ---
    unsigned short _polyDegreeBasis; /*!< @brief Max polynomial degree of the basis */

    // Moment basis
    SphericalBase* _basis; /*!< @brief Class to compute and store current spherical harmonics basis */

    // Quadrature related members
    VectorVector _quadPoints;       /*!<  @brief quadrature points, dim(_quadPoints) = (_nq,spatialDim) */
    Vector _weights;                /*!<  @brief quadrature weights, dim(_weights) = (_nq) */
    VectorVector _quadPointsSphere; /*!<  @brief (my,phi), dim(_quadPoints) = (_nq,2) */

  public:
    /**
     * @brief CSDPNSolver constructor
     * @param settings stores all needed information
     */
    CSDPNSolver_JL( Config* settings );

    virtual ~CSDPNSolver_JL();

    // virtual Solve() override;

  private:
    void SolverPreprocessing() override;

    void IterPreprocessing( unsigned /*idx_iter*/ ) override;
    void IterPostprocessing( unsigned /*idx_iter*/ ) override;

    void FluxUpdate() override;
    void FVMUpdate( unsigned idx_energy ) override;

    void PrepareVolumeOutput() override;
    void WriteVolumeOutput( unsigned idx_pseudoTime ) override;

    Vector ConstructFlux( unsigned idx_cell );

    // Helper Functions for CSD
    double NormPDF( double x, double mu, double sigma );
    Vector Time2Energy( const Vector& t, const double E_CutOff );
    double Time2Energy( const double t, const double E_CutOff );
    Vector Energy2Time( const Vector& E, const double E_CutOff );
    double Energy2Time( const double E, const double E_CutOff );
};

#endif    // CSDPNSOLVER_JL_H
