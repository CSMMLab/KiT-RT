#ifndef QSPHERICALTRIANGLE_H
#define QSPHERICALTRIANGLE_H

#include "quadraturebase.hpp"

class QSphericalTessalation : public QuadratureBase
{
    // Implementation is done accordingly to Kendall Atkinson 1981, Australian Matematical Society.

  public:
    QSphericalTessalation( Config* settings );

    virtual ~QSphericalTessalation() {}

  protected:
    virtual bool CheckOrder();
    virtual inline void SetName() override { _name = "Tesselated Spherical Triangle Quadrature"; }
    void SetNq() override;
    void SetPointsAndWeights() override;
    void SetConnectivity() override;

  private:
    std::array<double, 3> compute_centroid( const std::array<std::array<double, 3>, 3>& triangle );
    std::array<double, 3> map_to_unit_sphere( const std::array<double, 3>& point );
    double dot_product( const std::array<double, 3>& v1, const std::array<double, 3>& v2 );
    double angle_between_vectors( const std::array<double, 3>& a, const std::array<double, 3>& b, const std::array<double, 3>& c );
    double spherical_triangle_area( const std::array<double, 3>& a, const std::array<double, 3>& b, const std::array<double, 3>& c );
    std::array<double, 3> midpoint( const std::array<double, 3>& p1, const std::array<double, 3>& p2 );
    void reflect_and_permute( const std::vector<std::array<double, 3>>& points,
                              const std::vector<double>& weights,
                              std::vector<std::array<double, 3>>& full_points,
                              std::vector<double>& full_weights );

    std::vector<std::array<std::array<double, 3>, 3>> generate_tessellation( const std::array<std::array<double, 3>, 3>& triangle, int order );
};

#endif    // QSPHERICALTRIANGLE_H
