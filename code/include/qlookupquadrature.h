#ifndef QLOOKUPQUADRATURE_H
#define QLOOKUPQUADRATURE_H

#include "quadrature.h"

class QLookupQuadrature : public Quadrature
{
  public:
    QLookupQuadrature( unsigned order );
    virtual ~QLookupQuadrature() {}

    //helper
    inline std::vector<unsigned> getAvailOrders() const {return _availableOrders;}    /*! @returns: Vector with avalable orders for level symmetric quadrature */

  protected:
    bool CheckOrder();        /*! @brief checks if given order is available for this quadrature rule. */
    void SetNq() override;    /*! @brief: Assumes, that _order is available in lookup table */
    void SetPointsAndWeights() override; /*! @brief reads in n_points gridpoints and -weights from given filename. */

    virtual void SetAvailOrders() = 0; /*! @brief Sets vector with avaialbe Orders and corresponding vector with Nq. */
    virtual void SetDataInfo() = 0;    /*! @brief Sets data path for lookup table */

    std::vector<unsigned> _availableOrders; /*! @brief: Vector with available orders for lookup table */
    std::vector<unsigned> _nqByOrder; /*! @brief: Vector with number of quadrature points of each listed order */
    std::string _dataFiles; /*! path to the lookup table */
    std::string _dataFileSuffix; /*! suffix to data files of this quadrature */

};
#endif // QLOOKUPQUADRATURE_H
