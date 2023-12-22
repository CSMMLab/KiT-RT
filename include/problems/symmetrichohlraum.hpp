#ifndef SYMMETRIC_HOHLRAUM_H
#define SYMMETRIC_HOHLRAUM_H

#include "problems/problembase.hpp"

class SymmetricHohlraum : public ProblemBase
{
  private:
    SymmetricHohlraum() = delete;

    std::vector<unsigned> linspace2D( const std::vector<double>& start, const std::vector<double>& end, unsigned num_points );

  protected:
    Vector _sigmaS; /*!< @brief Vector of scattering crosssections */
    Vector _sigmaT; /*!< @brief Vector of total crosssections */

    std::vector<double> _cornerUpperLeftGreen;  /*!< @brief Coord of corner of the green area (minus thickness/2 of it) */
    std::vector<double> _cornerLowerLeftGreen;  /*!< @brief Coord of corner of the green area (minus thickness/2 of it) */
    std::vector<double> _cornerUpperRightGreen; /*!< @brief Coord of corner of the green area (minus thickness/2 of it) */
    std::vector<double> _cornerLowerRightGreen; /*!< @brief Coord of corner of the green area (minus thickness/2 of it) */
    double _thicknessGreen;                     /*!< @brief thickness of the green area */

    void SetGhostCells() override; /*!< @brief Sets vector of ghost cells for boundary conditions */

    void SetProbingCellsLineGreen(); /*!< @brief Sets cell ids for probing cells on the green line of the hohlraum */

    double _curAbsorptionHohlraumCenter;       /*!< @brief Absorption of particles at Hohlraum center at current time step  */
    double _curAbsorptionHohlraumVertical;     /*!< @brief Absorption of particles at Hohlraum vertical walls at current time step */
    double _curAbsorptionHohlraumHorizontal;   /*!< @brief Absorption of particles at Hohlraum horizontal walls at current time step */
    double _totalAbsorptionHohlraumCenter;     /*!< @brief Absorption of particles at Hohlraum center integrated until current time step  */
    double _totalAbsorptionHohlraumVertical;   /*!< @brief Absorption of particles at Hohlraum vertical walls integrated until current time step */
    double _totalAbsorptionHohlraumHorizontal; /*!< @brief Absorption of particles at Hohlraum horizontal walls integrated until current time step */
    double _varAbsorptionHohlraumGreen;        /*!< @brief Absorption of particles at Hohlraum green center cells integrated at current time step */
    std::vector<unsigned> _probingCells;       /*!< @brief Indices of cells that contain a probing sensor */
    VectorVector _probingMoments;              /*!< @brief Solution Momnets at the probing cells that contain a probing sensor */
    unsigned _nProbingCellsLineGreen;          /*!< @brief Number of sampling cells that contain a probing sensor for the sliding window */
    std::vector<unsigned> _probingCellsLineGreen;     /*!< @brief Indices of cells that contain a probing sensor for the sliding window */
    std::vector<double> _absorptionValsIntegrated;    /*!< @brief Avg Absorption value at the sampleing points of lineGreen */
    std::vector<double> _varAbsorptionValsIntegrated; /*!< @brief Var in Avg Absorption value at the sampleing points of lineGreen */

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
    void ComputeCurrentProbeMoment( const VectorVector& solution ) override;    /*!<   @brief Computes Problemspecific Scalar QOI */
    void ComputeVarAbsorptionGreen( const Vector& scalarFlux ) override;        /*!<   @brief Computes Problemspecific Scalar QOI */
    void ComputeQOIsGreenProbingLine( const Vector& scalarFlux ) override;

    double GetCurAbsorptionHohlraumCenter() override { return _curAbsorptionHohlraumCenter; };
    double GetCurAbsorptionHohlraumVertical() override { return _curAbsorptionHohlraumVertical; };
    double GetCurAbsorptionHohlraumHorizontal() override { return _curAbsorptionHohlraumHorizontal; };
    double GetTotalAbsorptionHohlraumCenter() override { return _totalAbsorptionHohlraumCenter; };
    double GetTotalAbsorptionHohlraumVertical() override { return _totalAbsorptionHohlraumVertical; };
    double GetTotalAbsorptionHohlraumHorizontal() override { return _totalAbsorptionHohlraumHorizontal; };
    double GetVarAbsorptionHohlraumGreen() override { return _varAbsorptionHohlraumGreen; };
    const VectorVector& GetCurrentProbeMoment() const override { return _probingMoments; };
    virtual const std::vector<double>& GetCurrentProbeValuesGreenLine() const override { return _absorptionValsIntegrated; };
    virtual const std::vector<double>& GetCurrentVarProbeValuesGreenLine() const override { return _varAbsorptionValsIntegrated; };
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
    void ComputeCurrentProbeMoment( const VectorVector& solution ) override final; /*!<   @brief Computes Problemspecific Scalar QOI */
};

#endif    // SYMMETRIC_HOHLRAUM_H
