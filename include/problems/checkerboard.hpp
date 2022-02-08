#ifndef CHECKERBOARD_H
#define CHECKERBOARD_H

#include "problembase.hpp"

class SphericalBase;

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

class Checkerboard_Moment : public ProblemBase
{
  private:
    Vector _scatteringXS; /*!< @brief Vector of scattering crosssections len: numCells  */
    Vector _totalXS;      /*!< @brief Vector of total crosssections. len: numCells*/

    Checkerboard_Moment() = delete;

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
    Checkerboard_Moment( Config* settings, Mesh* mesh );
    virtual ~Checkerboard_Moment();

    virtual VectorVector GetScatteringXS( const Vector& energies ) override;
    virtual VectorVector GetTotalXS( const Vector& energies ) override;
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
};

class Checkerboard_SN_1D : public ProblemBase
{
  private:
    Vector _scatteringXS; /*!< @brief Vector of scattering crosssections */
    Vector _totalXS;      /*!< @brief Vector of total crosssections */

    Checkerboard_SN_1D() = delete;

    bool isAbsorption( const Vector& pos ) const; /*!< @return True if pos is in absorption region, False otherwise */
    bool isSource( const Vector& pos ) const;     /*!< @return True if pos is in source region, False otherwise */

  public:
    Checkerboard_SN_1D( Config* settings, Mesh* mesh );
    virtual ~Checkerboard_SN_1D();

    virtual VectorVector GetScatteringXS( const Vector& energies ) override;
    virtual VectorVector GetTotalXS( const Vector& energies ) override;
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
};

class Checkerboard_Moment_1D : public ProblemBase
{
  private:
    Vector _scatteringXS; /*!< @brief Vector of scattering crosssections len: numCells  */
    Vector _totalXS;      /*!< @brief Vector of total crosssections. len: numCells*/

    Checkerboard_Moment_1D() = delete;

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
    Checkerboard_Moment_1D( Config* settings, Mesh* mesh );
    virtual ~Checkerboard_Moment_1D();

    virtual VectorVector GetScatteringXS( const Vector& energies ) override;
    virtual VectorVector GetTotalXS( const Vector& energies ) override;
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
};

#endif    // CHECKERBOARD_H
