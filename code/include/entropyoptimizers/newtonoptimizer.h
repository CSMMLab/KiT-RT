//#ifndef NEWTONOPTIMIZER_H
//#define NEWTONOPTIMIZER_H
//
//#include "entropyfunctional.h"
//#include "settings/config.h"
//#include <vector>
//
// class NewtonOptimizer
//{
//  public:
//    NewtonOptimizer( Config* settings );
//
//    /*! \brief: Solves the optimization problem for given moments*/
//    void Solve( std::vector<std::vector<double>>* moments );
//
//    /*! \brief Get pointer to solution vector for each cell */
//    std::vector<std::vector<double>>* GetSolution() { return &_lambda; }
//
//  protected:
//    Config* _settings; /*! \brief: Settings file */
//
//    // EntropyFunctionalBase _entropyFunctional; /*! \brief: entropy functional */
//    std::vector<std::vector<double>> _lambda; /*! \brief: Solution of the optimization problem */
//
//    int _nCells;    /*! \brief: Number of cells, i.e. outer length of vectorvector _lambda */
//    int _nMomentEq; /*! \brief: Number of moment equations, i.e. length of vector _lambda per cell */
//};
//
//#endif    // NEWTONOPTIMIZER_H
