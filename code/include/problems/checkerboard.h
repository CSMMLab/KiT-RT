#ifndef CHECKERBOARD_H
#define CHECKERBOARD_H

#include "common/typedef.h"    // for Vector, VectorVector
#include "problembase.h"
#include <vector>    // for vector
class Config;
class Mesh;

class Checkerboard_SN : public ProblemBase
{
  private:
    Vector _scatteringXS; /*!< @brief Vector of scattering crosssections */
    Vector _totalXS;      /*!< @brief Vector of total crosssections */

    Checkerboard_SN() = delete;

    bool isAbsorption( const Vector& pos ) const; /*!< @return True if pos is in absorption region, False otherwise */
    bool isSource( const Vector& pos ) const;     /*!< @return True if pos is in source region, False otherwise */

  public:
    Checkerboard_SN( Config* settings, Mesh* mesh );
    virtual ~Checkerboard_SN();

    virtual VectorVector GetScatteringXS( const Vector& energies ) override;
    virtual VectorVector GetTotalXS( const Vector& energies ) override;
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
};

class Checkerboard_PN : public ProblemBase
{
  private:
    Vector _scatteringXS; /*!< @brief Vector of scattering crosssections len: numCells  */
    Vector _totalXS;      /*!< @brief Vector of total crosssections. len: numCells*/

    Checkerboard_PN() = delete;

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
    Checkerboard_PN( Config* settings, Mesh* mesh );
    virtual ~Checkerboard_PN();

    virtual VectorVector GetScatteringXS( const Vector& energies ) override;
    virtual VectorVector GetTotalXS( const Vector& energies ) override;
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
};

#endif    // CHECKERBOARD_H
