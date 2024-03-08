#ifndef HALFLATTICE_H
#define HALFLATTICE_H

#include "problembase.hpp"

class SphericalBase;

class HalfLattice_SN : public ProblemBase
{
  private:
    HalfLattice_SN() = delete;

    Vector _sigmaS; /*!< @brief Vector of scattering crosssections */
    Vector _sigmaT; /*!< @brief Vector of total crosssections */

    std::map<unsigned, bool> _ghostCellsReflectingY;     /*!< map that indicates if a ghostcell has a fixed value or is a mirroring boundary */
    std::map<unsigned, unsigned> _quadratureYReflection; /*!< map that gives a Reflection against the y axis for the velocity ordinates */

    // Lattice QOIS
    double _curAbsorptionLattice;    /*!< @brief Absorption of particles at Lattice checkerboard regions at current time step */
    double _totalAbsorptionLattice;  /*!< @brief Absorption of particles at Lattice checkerboard regions integrated until current time step */
    double _curMaxAbsorptionLattice; /*!< @brief Maximum pointwise absorption of particles at Lattice checkerboard regions  until current time step */

    bool IsAbsorption( const Vector& pos ) const; /*!< @return True if pos is in absorption region, False otherwise */
    bool IsSource( const Vector& pos ) const;     /*!< @return True if pos is in source region, False otherwise */

  protected:
    void SetGhostCells() override;                  /*!< @brief Sets vector of ghost cells for boundary conditions */
    unsigned GetBlockID( const Vector& pos ) const; /*!< @brief Returns checkerboard field id (0-48, row major) of the Lattice test case */

  public:
    HalfLattice_SN( Config* settings, Mesh* mesh, QuadratureBase* quad );
    virtual ~HalfLattice_SN();

    virtual VectorVector GetScatteringXS( const Vector& energies ) override;
    virtual VectorVector GetTotalXS( const Vector& energies ) override;
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;

    const Vector& GetGhostCellValue( int idx_cell, const Vector& cell_sol ) override final;

    void ComputeCurrentOutflow( const VectorVector& solution ) override;         /*!<   @brief Computes Problemspeagnostic Scalar QOI */
    void ComputeMaxOrdinatewiseOutflow( const VectorVector& solution ) override; /*!<   @brief Computes Problemspeagnostic Scalar QOI */

    double GetCurAbsorptionLattice() override final;
    double GetTotalAbsorptionLattice() override final;
    double GetMaxAbsorptionLattice() override final;
    void ComputeTotalAbsorptionLattice( double dT ) override;
    void ComputeCurrentAbsorptionLattice( const Vector& scalarFlux ) override;
    void ComputeMaxAbsorptionLattice( const Vector& scalarFlux ) override;
};

class HalfLattice_Moment : public ProblemBase
{
  private:
    Vector _sigmaS; /*!< @brief Vector of scattering crosssections len: numCells  */
    Vector _sigmaT; /*!< @brief Vector of total crosssections. len: numCells*/

    HalfLattice_Moment() = delete;

    bool isAbsorption( const Vector& pos ) const; /*!< @return True if pos is in absorption region, False otherwise */
    bool isSource( const Vector& pos ) const;     /*!< @return True if pos is in source region, False otherwise */

    /*!
     * @brief Gets the global index for given order l of Legendre polynomials and given
     *        order k of Legendre functions.
     *        Note: This is code doubling from PNSolver::GlobalIndex
     * @param l  order of Legendre polynomial
     * @param k  order of Legendre function
     * @returns globalIndex
     */
    int GlobalIndex( int l, int k ) const;

  public:
    HalfLattice_Moment( Config* settings, Mesh* mesh, QuadratureBase* quad );
    virtual ~HalfLattice_Moment();

    virtual VectorVector GetScatteringXS( const Vector& energies ) override;
    virtual VectorVector GetTotalXS( const Vector& energies ) override;
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
};

#endif    // HALFLATTICE_H
