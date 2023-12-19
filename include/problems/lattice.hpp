#ifndef LATTICE_H
#define LATTICE_H

#include "problembase.hpp"

class SphericalBase;

class Lattice_SN : public ProblemBase
{
  private:
    Vector _scatteringXS; /*!< @brief Vector of scattering crosssections */
    Vector _totalXS;      /*!< @brief Vector of total crosssections */

    Lattice_SN() = delete;

    bool IsAbsorption( const Vector& pos ) const; /*!< @return True if pos is in absorption region, False otherwise */
    bool IsSource( const Vector& pos ) const;     /*!< @return True if pos is in source region, False otherwise */

  protected:
    void SetGhostCells() override;                  /*!< @brief Sets vector of ghost cells for boundary conditions */
    unsigned GetBlockID( const Vector& pos ) const; /*!< @brief Returns checkerboard field id (0-48, row major) of the Lattice test case */

  public:
    Lattice_SN( Config* settings, Mesh* mesh );
    virtual ~Lattice_SN();

    virtual VectorVector GetScatteringXS( const Vector& energies ) override;
    virtual VectorVector GetTotalXS( const Vector& energies ) override;
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;

    const Vector& GetGhostCellValue( int idx_cell, const Vector& cell_sol ) override final;
};

class Lattice_Moment : public ProblemBase
{
  private:
    Vector _scatteringXS; /*!< @brief Vector of scattering crosssections len: numCells  */
    Vector _totalXS;      /*!< @brief Vector of total crosssections. len: numCells*/

    Lattice_Moment() = delete;

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
    Lattice_Moment( Config* settings, Mesh* mesh );
    virtual ~Lattice_Moment();

    virtual VectorVector GetScatteringXS( const Vector& energies ) override;
    virtual VectorVector GetTotalXS( const Vector& energies ) override;
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
};

#endif    // LATTICE_H
