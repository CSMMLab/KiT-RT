#ifndef WATERPHANTOM_H
#define WATERPHANTOM_H

#include "problems/problembase.hpp"

class Waterphantom : public ProblemBase
{
  private:
    Waterphantom() = delete;

  public:
    Waterphantom( Config* settings, Mesh* mesh );
    ~Waterphantom();
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    VectorVector SetupIC() override;
    std::vector<double> GetDensity( const VectorVector& cellMidPoints ) override final;
    VectorVector GetScatteringXS( const Vector& energies ) override final;
    VectorVector GetTotalXS( const Vector& energies ) override final;
};

class Waterphantom_Moment : public Waterphantom
{
  private:
    Waterphantom_Moment() = delete;

  public:
    Waterphantom_Moment( Config* settings, Mesh* mesh );
    ~Waterphantom_Moment();
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override final;
    VectorVector SetupIC() override final;
};

#endif    // WATERPHANTOM_H
