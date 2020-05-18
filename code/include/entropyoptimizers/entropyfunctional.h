#ifndef ENTROPYFUNCTIONAL_H
#define ENTROPYFUNCTIONAL_H

class EntropyFunctionalBase
{
  public:
    EntropyFunctionalBase() {}

    static EntropyFunctionalBase* CreateEntropyFunctional();

    /*! \brief: Computes the value of the specific entropy functional. */
    double virtual ComputeEntropy( double x ) = 0;

    /*! \brief: Computes the sensitivity value of the specific entropy functional. */
    double virtual ComputeEntropySens( double x ) = 0;

    /*! \brief: Computes the value of the Legendre dual of the specific entropy functional. */
    double virtual ComputeLegendreDual( double x ) = 0;

    /*! \brief: Computes the sensitivity value ofthe Legendre dual of the specific entropy functional. */
    double virtual ComputeLegendreDualSens( double x ) = 0;

  protected:
};

class QuadraticEntropy : EntropyFunctionalBase
{
  public:
    QuadraticEntropy() {}

    /*! \brief: Computes the value of the specific entropy functional. */
    double ComputeEntropy( double x ) override { return 0.5 * x * x; }

    /*! \brief: Computes the sensitivity value of the specific entropy functional. */
    double virtual ComputeEntropySens( double x ) override { return x; }

    /*! \brief: Computes the value of the Legendre dual of the specific entropy functional. */
    double virtual ComputeLegendreDual( double x ) override { return 0.5 * x * x; }

    /*! \brief: Computes the sensitivity value ofthe Legendre dual of the specific entropy functional. */
    double virtual ComputeLegendreDualSens( double x ) override { return x; }

  protected:
};

#endif    // ENTROPYFUNCTIONAL_H
