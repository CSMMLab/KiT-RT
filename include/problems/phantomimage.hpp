#ifndef PHANTOMIMAGE_H
#define PHANTOMIMAGE_H

#include "problems/problembase.hpp"

class PhantomImage : public ProblemBase
{
  private:
    PhantomImage() = delete;
    double _sigmaS; /*!< @brief Scattering coefficient */
  public:
    PhantomImage( Config* settings, Mesh* mesh, QuadratureBase* quad );
    virtual ~PhantomImage();
    VectorVector SetupIC() override;
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    std::vector<double> GetDensity( const VectorVector& cellMidPoints ) override;
    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
};

#endif    // PHANTOMIMAGE_H
