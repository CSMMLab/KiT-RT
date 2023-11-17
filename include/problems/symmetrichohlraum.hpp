#ifndef SYMMETRIC_HOHLRAUM_H
#define SYMMETRIC_HOHLRAUM_H

#include "problems/problembase.hpp"

class SymmetricHohlraum : public ProblemBase
{
  private:
    SymmetricHohlraum() = delete;
    Vector _scatteringXS; /*!< @brief Vector of scattering crosssections */
    Vector _totalXS;      /*!< @brief Vector of total crosssections */

  protected:
    void SetGhostCells() override ; /*!< @brief Sets vector of ghost cells for boundary conditions */

  public:
    SymmetricHohlraum( Config* settings, Mesh* mesh );
    virtual ~SymmetricHohlraum();
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
    VectorVector GetScatteringXS( const Vector& energies ) override final;
    VectorVector GetTotalXS( const Vector& energies ) override final;
    const Vector& GetGhostCellValue(int idx_cell,  const Vector& cell_sol) override final;
};

class SymmetricHohlraum_Moment : public SymmetricHohlraum
{
  private:
    SymmetricHohlraum_Moment() = delete;
       // TODO void SetGhostCells() override final; /*!< @brief Sets vector of ghost cells for boundary conditions */

  public:
    SymmetricHohlraum_Moment( Config* settings, Mesh* mesh );
    ~SymmetricHohlraum_Moment();
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override final;
    VectorVector SetupIC() override final;
};

#endif    // SYMMETRIC_HOHLRAUM_H
