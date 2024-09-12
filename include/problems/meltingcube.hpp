#ifndef MELTINGCUBE_H
#define MELTINGCUBE_H

#include "problembase.hpp"

class MeltingCube : public ProblemBase
{
  private:
    MeltingCube() = delete;

  protected:
    double _sigmaS; /*!< @brief Scattering coefficient */

  public:
    MeltingCube( Config* settings, Mesh* mesh, QuadratureBase* quad );
    ~MeltingCube();
    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
};

class MeltingCube_SN : public MeltingCube
{
  private:
    MeltingCube_SN() = delete;

  public:
    MeltingCube_SN( Config* settings, Mesh* mesh, QuadratureBase* quad );
    ~MeltingCube_SN();
    VectorVector SetupIC() override;
};

class MeltingCube_Moment : public MeltingCube
{
  private:
    MeltingCube_Moment() = delete;

  public:
    MeltingCube_Moment( Config* settings, Mesh* mesh, QuadratureBase* quad );
    ~MeltingCube_Moment();
    VectorVector SetupIC() override;
};

class MeltingCube_SN_1D : public MeltingCube
{
  private:
    MeltingCube_SN_1D() = delete;

  public:
    MeltingCube_SN_1D( Config* settings, Mesh* mesh, QuadratureBase* quad );
    ~MeltingCube_SN_1D();

    VectorVector SetupIC() override;
};

class MeltingCube_Moment_1D : public MeltingCube
{
  private:
    MeltingCube_Moment_1D() = delete;

  public:
    MeltingCube_Moment_1D( Config* settings, Mesh* mesh, QuadratureBase* quad );
    ~MeltingCube_Moment_1D();

    VectorVector SetupIC() override;
};

#endif    // MELTINGCUBE_H
