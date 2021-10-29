//
// Created by chinsp on 29/10/21.
//

#ifndef REFINESNSOLVER_H
#define REFINESNSOLVER_H

#include "solvers/snsolver.h"

class RefineSNSolver : public SNSolver
{
public:
    RefineSNSolver( Config* settings );
    ~RefineSNSolver() {}

    void Solve() override;
    void IterPreprocessing( unsigned iter ) override;

private:
    Matrix GetScatteringKernel( unsigned nq, Vector weights);   /*! Have to use new function because in original we cannot give _nq and _w */
    void CalculateInterpolation( unsigned quadIter);            /*! @brief Calculates Interpolation between 3 nearest neighbors */
    void InterpolateFlux( );
    void AddArtificialScattering( unsigned iter );

    clock_t endbuild;

    unsigned _nqF;              /*! @brief Nr of points in Fine Quadrature */
    unsigned _nqC;              /*! @brief Nr of points in Coarse Quadrature */
    VectorVector _quadPointsF;  /*! @brief Given points in Fine Quadrature */
    VectorVector _quadPointsC;  /*! @brief Given points in Coarse Quadrature */
    Vector _weightsF;           /*! @brief Given weights of points in Fine Quadrature */
    Vector _weightsC;           /*! @brief Given weights of points in Coarsee Quadrature */
    // VectorVectorU _neighF;
    VectorVectorU _neighC;      /*! @brief Neighbors of Quadrature points in coarse Grid */

    VectorVector _solF;

    VectorVectorU _quadC2F; /*! @brief VectorVector with fine quadrature points for each coarse quadrature point ( dim = [_nqC][nr of corresp. fine points]) */
    // VectorU _quadF2C;

    VectorU _refine;        /*! @brief Vector which tells in which direction we have to refine (dim = _nqC) */
    VectorU _refineIter;    /*! @brief Vector which tells how many Iterations are done with refinement (reset after given nr of Iterations) */
    VectorU _refineIterOld;

    VectorU _quadIDsOld;    /*! @brief _quadIDs from last Iteration*/
    VectorU _quadIDs;       /*! @brief Vector with IDs of quadrature points in original quadratures ( dim = _nq ) */

    blaze::CompressedMatrix<double> _solInter;       /*! @brief interpolation matrix (dim = (_nq, _nqOld)) */

    std::vector<VectorVector> _QextF;    /*!  @brief  external source term in Fine quadrature */
    std::vector<VectorVector> _QextC;    /*!  @brief  external source term in Coarse quadrature */

    unsigned _nqOld;                /*! @brief Nr of quadrature points in last Iteration */
    VectorVector _quadPointsOld;    /*! @brief Given quadrature points in last Iteration */

    // have to put this 3  in solverbase
    double _sigmaAS;
    double _beta;
    Matrix _ASKernel;       /*! @brief Kernel of Artifial Scattering ( dim = (_nq, _nq)) */

    // These are just debugging helpers
    clock_t start, startiter, t, tt, t1, t2, t21, t22, t23, t24, t3 , t4;
    FILE* Times;
    FILE* radFlux;
    FILE*  Refine;
    FILE* Info;
    FILE* ASFile;
};

#endif //REFINESNSOLVER_H
