#ifndef SYMMETRIC_HOHLRAUM_H
#define SYMMETRIC_HOHLRAUM_H

#include "problems/problembase.hpp"

class SymmetricHohlraum : public ProblemBase
{
  private:
    SymmetricHohlraum() = delete;
    Vector _sigmaS; /*!< @brief Vector of scattering crosssections */
    Vector _sigmaT; /*!< @brief Vector of total crosssections */

  protected:
    void SetGhostCells() override; /*!< @brief Sets vector of ghost cells for boundary conditions */

    double _curAbsorptionHohlraumCenter;       /*!< @brief Absorption of particles at Hohlraum center at current time step  */
    double _curAbsorptionHohlraumVertical;     /*!< @brief Absorption of particles at Hohlraum vertical walls at current time step */
    double _curAbsorptionHohlraumHorizontal;   /*!< @brief Absorption of particles at Hohlraum horizontal walls at current time step */
    double _totalAbsorptionHohlraumCenter;     /*!< @brief Absorption of particles at Hohlraum center integrated until current time step  */
    double _totalAbsorptionHohlraumVertical;   /*!< @brief Absorption of particles at Hohlraum vertical walls integrated until current time step */
    double _totalAbsorptionHohlraumHorizontal; /*!< @brief Absorption of particles at Hohlraum horizontal walls integrated until current time step */
    double _varAbsorptionHohlraumGreen;        /*!< @brief Absorption of particles at Hohlraum green center cells integrated at current time step */

  public:
    SymmetricHohlraum( Config* settings, Mesh* mesh, QuadratureBase* quad );
    virtual ~SymmetricHohlraum();
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
    VectorVector GetScatteringXS( const Vector& energies ) override final;
    VectorVector GetTotalXS( const Vector& energies ) override final;
    const Vector& GetGhostCellValue( int idx_cell, const Vector& cell_sol ) override final;

    void ComputeCurrentAbsorptionHohlraum( const Vector& scalarFlux ) override; /*!<   @brief Computes Problemspecific Scalar QOI */
    void ComputeTotalAbsorptionHohlraum( double dT ) override;                  /*!<   @brief Computes Problemspecific Scalar QOI */
    // void ComputeCurrentProbeMoment() override;                                  /*!<   @brief Computes Problemspecific Scalar QOI */
    void ComputeVarAbsorptionGreen( const Vector& scalarFlux ) override; /*!<   @brief Computes Problemspecific Scalar QOI */
    double GetCurAbsorptionHohlraumCenter() { return _curAbsorptionHohlraumCenter; };
    double GetCurAbsorptionHohlraumVertical() { return _curAbsorptionHohlraumVertical; };
    double GetCurAbsorptionHohlraumHorizontal() { return _curAbsorptionHohlraumHorizontal; };
    double GetTotalAbsorptionHohlraumCenter() { return _totalAbsorptionHohlraumCenter; };
    double GetTotalAbsorptionHohlraumVertical() { return _totalAbsorptionHohlraumVertical; };
    double GetTotalAbsorptionHohlraumHorizontal() { return _totalAbsorptionHohlraumHorizontal; };
    double GetVarAbsorptionHohlraumGreen() { return _varAbsorptionHohlraumGreen; };
};

class SymmetricHohlraum_Moment : public SymmetricHohlraum
{
  private:
    SymmetricHohlraum_Moment() = delete;
    // TODO void SetGhostCells() override final; /*!< @brief Sets vector of ghost cells for boundary conditions */

  public:
    SymmetricHohlraum_Moment( Config* settings, Mesh* mesh, QuadratureBase* quad );
    ~SymmetricHohlraum_Moment();
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override final;
    VectorVector SetupIC() override final;
};

#endif    // SYMMETRIC_HOHLRAUM_H
