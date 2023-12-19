#ifndef HOHLRAUM_H
#define HOHLRAUM_H

#include "problems/problembase.hpp"

class Hohlraum : public ProblemBase
{
  private:
    Hohlraum() = delete;
    Vector _scatteringXS; /*!< @brief Vector of scattering crosssections */
    Vector _totalXS;      /*!< @brief Vector of total crosssections */

  public:
    Hohlraum( Config* settings, Mesh* mesh, QuadratureBase* quad );
    virtual ~Hohlraum();
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
    VectorVector GetScatteringXS( const Vector& energies ) override final;
    VectorVector GetTotalXS( const Vector& energies ) override final;
};

class Hohlraum_Moment : public Hohlraum
{
  private:
    Hohlraum_Moment() = delete;

  public:
    Hohlraum_Moment( Config* settings, Mesh* mesh, QuadratureBase* quad );
    ~Hohlraum_Moment();
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override final;
    VectorVector SetupIC() override final;
};

#endif    // HOHLRAUM_H
