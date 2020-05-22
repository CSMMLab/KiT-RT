///*!
// * \file entropyclosuredriver.h
// * \brief The main subroutines for driving entropyclosure routines.
// * \author S. Schotthöfer
// *
// */
//
//#ifndef ENTROPYCLOSUREDRIVER_H
//#define ENTROPYCLOSUREDRIVER_H
//
//#include "entropyoptimizers/newtonoptimizer.h"
//#include "settings/config.h"
//#include "solvers/advectionsolver.h"
//
// class EntropyClosureDriver
//{
//  public:
//    EntropyClosureDriver( Config* settings );
//
//    /*!
//     * \brief Launch the computation.
//     */
//    void StartSolver();
//
//  protected:
//    Config* _driverSettings;
//    double _startTime;       /*! \brief: Wall time at start of computation */
//    double _stopTime;        /*! \brief: Wall time at end of computation */
//    unsigned long _timeIter; /*! \brief: Physical time iteration */
//    double _nInnerIter;      /*! \brief: Maximum amount of inner Iterations */
//
//    std::vector<double>* _kineticDensity;       /*! \brief: pointer to kinetic density (for each cell) */
//    std::vector<std::vector<double>>* _moments; /*! \brief: pointer to moment vector (for each cell), i.e solution of _solver */
//    std::vector<std::vector<double>>* _lambda;  /*! \brief: pointer to vector of lagrange multipliers (for each cell), i.e. solution of _optimizer
//    */ NewtonOptimizer* _optimizer;                /*! \brief optimizer to reconstruct the entropy density */ AdvectionSolver* _solver; /*! \brief
//    solver to solve the entropy moment system */
//
//    /*!
//     * \brief Preprocess the iteration
//     */
//    void Preprocess( unsigned long timeIter );
//
//    /*!
//     * \brief Run the physical time iteration of the solver.
//     */
//    void Run();
//
//    /*!
//     * \brief Postprocess the iteration.
//     */
//    void Postprocess();
//
//    /*!
//     * \brief Output the solution in solution file.
//     */
//    void Output( unsigned long timeIter );
//
//    /*!
//     * \brief Monitor Convergence of the solution.
//     */
//    bool Monitor( unsigned long timeIter );
//};
//
//#endif    // ENTROPYCLOSUREDRIVER_H
