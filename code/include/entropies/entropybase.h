#ifndef ENTROPYBASE_H
#define ENTROPYBASE_H

#include "settings/config.h"

class EntropyBase
{
  public:
    inline EntropyBase() {}

    virtual inline ~EntropyBase() {}

    static EntropyBase* Create( Config* settings );

    /*! @brief: computes the entropy functional
     *  @param: z = point where the functional should be evaluated.
     *          z must be in domain of the functional
     *  @returns: value of entropy functional at z */
    virtual double Entropy( double z ) = 0;

    /*! @brief: computes the derivative of the entropy functional
     *  @param: z = point where the derivative should be evaluated.
     *          z must be in domain of the functional
     *  @returns: value of derivative at z */
    virtual double EntropyPrime( double z ) = 0;

    /*! @brief: computes the dual of the  entropy functional
     *  @param: z = point where the dual of the functional should be evaluated.
     *          z must be in domain of the duaö functional
     *  @returns: value of entropy functional at z */
    virtual double EntropyDual( double y ) = 0;

    /*! @brief: computes the  dual of the derivative of the entropy functional
     *  @param: z = point where the dual of the  derivative should be evaluated.
     *          z must be in domain of the dual functional
     *  @returns: value of dual of the derivative at z */
    virtual double EntropyPrimeDual( double y ) = 0;

    /*! @brief: checks, if value is in domain of entropy */
    virtual bool CheckDomain( double z ) = 0;
};

#endif    // ENTROPYBASE_H
