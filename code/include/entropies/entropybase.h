#ifndef ENTROPYBASE_H
#define ENTROPYBASE_H

// Foward declaration
class Config;

class EntropyBase
{
  public:
    inline EntropyBase() {}

    virtual inline ~EntropyBase() {}

    static EntropyBase* Create( Config* settings );

    /*! @brief computes the entropy functional
     *  @param z = point where the functional should be evaluated.
     *          z must be in domain of the functional
     *  @returns: value of entropy functional at z */
    virtual double Entropy( double z ) = 0;

    /*! @brief computes eta'(z).z must be in domain of the functional.
     *  @param z = point where the derivative should be evaluated. */
    virtual double EntropyPrime( double z ) = 0;

    /*! @brief computes the dual of the  entropy functional
     *  @param z = point where the dual of the functional should be evaluated.
     *          z must be in domain of the duaö functional
     *  @returns: value of entropy functional at z */
    virtual double EntropyDual( double y ) = 0;

    // /*! @brief computes m*eta_*'(alpha*m). alpha*m must be in domain of the functional
    //  *          alpha, m and grad must be of same size.
    //  *  @param alpha = point where the derivative should be evaluated.
    //  *  @param m = moment vector
    //  *  @param grad: vector in which the resulting gradient is stored */
    // virtual void EntropyPrimeDual( Vector& alpha, Vector& m, Vector& grad ) = 0;

    /*! @brief computes eta_*'(y).
     *  @param y = point where the derivative should be evaluated.
     *  @returns: value of the derivative of the  entropy functional at y */
    virtual double EntropyPrimeDual( double z ) = 0;

    /*! @brief computes the hessian of the dual entropy functional
     *  @param y = point where the hessian should be evaluated;
     *  @returns: value of the hessian at y */
    virtual double EntropyHessianDual( double y ) = 0;

    /*! @brief checks, if value is in domain of entropy */
    virtual bool CheckDomain( double z ) = 0;
};

#endif    // ENTROPYBASE_H
