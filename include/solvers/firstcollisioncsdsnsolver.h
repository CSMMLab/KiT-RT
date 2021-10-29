//
// Created by chinsp on 29/10/21.
//

#ifndef FIRSTCOLLISIONCSDSNSOLVER_H
#define FIRSTCOLLISIONCSDSNSOLVER_H

#include "solverbase.h"

class QuadratureBase;
class ScatteringKernel;

class FirstCollisionCSDSNSolver : public SolverBase( settings )
{

public:
    /*! @brief Constructor of FirstCollisonSNSolver */
    FirstCollisionCSDSNSolver( Config * settings );
    /*! @brief Destructor of FirstCollisionSNSolver */
    ~FirstCollisionCSDSNSolver() {}

    /*! @brief: Solver function */
    void Solve() override;

private:
    /*! @brief Performs preprocessing for the solver (CSD transformations) */
    void SolverPreprocessing( ) override;
    /*! @brief Performs preprocessing for the current solver iteration */
    void IterPreprocessing( unsigned idx_pseudotime ) override;
    /*! @brief Performs postprocessing for the current solver iteration */
    void IterPostprocessing( ) override; // not used
    void IterPostprocessing( unsigned iter );
    /*! @brief: Computes the flux of the solution to check conservation properties */
    void ComputeRadFlux() override;
    void ComputeRadFlux( unsigned iter );

    void FluxUpdate() override;
    /*! @brief Constructs the flux update for the first collision source of the current iteration and stores it in _solNewFC*/
    void FluxUpdateUncollided( );
    /*! @brief Constructs  the flux update for the current iteration and stores it in _solNew*/
    void FluxUpdateCollided( );

    void FVMUpdate( unsigned iter ) override;
    /*! @brief Computes the finite Volume update step for the first collision source of the current iteration */
    void FVMUpdateUncollided( unsigned iter );
    /*! @brief Computes the finite Volume update step of the current iteration */
    void FVMUpdateCollided( unsigned iter );

    /*! @brief Computes the first collision source _FCSource of the current iteration */
    void ComputeFirstCollisionSource( unsigned iter );
    /*! @brief: Adds the artificial scattering term to the updated first collision source */
    // void AddArtificialScattering( unsigned iter );
    void GetASTransportCoefficients();

    /*! @brief: Computation of the refined version of first collision _solFC */
    void IterRefinement( unsigned iter );
    /*! @brief: Refinement of the Quadrature depending on external source */
    void RefineQuadrature();
    /*! @brief: Interpolation of between solution on quadPoints in last iteartion and quadPoints in this iteration */
    void InterpolateSolution();
    /*! @brief: Determines if refinenment is necessary (per iteration step) */
    void DetermineRefinement( unsigned iter );

    /*! @brief: CConstructor of the angular flux */
    Vector ConstructFluxSN( unsigned idx_cell, bool FC );

    void Refinement();
    Matrix SetupLaplaceBeltrami( VectorVector p, Vector w );

    void PrepareVolumeOutput() override;
    void WriteVolumeOutput( unsigned idx_pseudoTime ) override;


    // Variables

private:
    // CSD Variables
    bool _RT;
    std::vector<double> _dose;
    Vector _energies;
    Vector _energiesOrig;
    Vector _angle;
    double _energyMax;

    std::vector<Matrix> _sigmaSE;
    std::vector<Matrix> _sigmaSECoarse;     /*!  @brief scattering cross section for all energies*/
    Vector _sigmaTE;                  /*!  @brief total cross section for all energies*/
    std::vector<Matrix> _sigmaSEFine;       /*!  @brief scattering cross section for all energies*/
    Vector _sigmaTEFine;                    /*!  @brief total cross section for all energies*/
    Vector _sigmaTECoarse;

    // FP Variables
    unsigned _FPMethod;
    blaze::DiagonalMatrix<Matrix> _identityFC;
    blaze::DiagonalMatrix<Matrix> _identity;
    Matrix _LC;
    Matrix _LFC;
    Matrix _xi;
    Matrix _xiF;
    Matrix _MatInterp;
    Matrix ILFC;
    Matrix IL;
    double alpha;
    double beta;
    double alpha2;
    double _sigmaAS;
    double _betaAS;
    Vector _xiAS;

    // general
    blaze::DiagonalMatrix<Matrix> _scatteringKernel;       /*! @brief Scattering Kernel of coarse quadrature grid */
    blaze::DiagonalMatrix<Matrix> _scatteringKernelFC;     /*! @brief Scattering kernel from fine to coarse mesh for FirstCollisionSource Calculation */

    VectorVector _solFC;            /*! @brief Stores the Flux of uncollided Part */
    VectorVector _solNewFC;         /*! @brief Helper Variable for flux of uncollided part */

    VectorVector _QFirstCollision;          /*! @brief First Collision Source (dim = [_nCells][_nqC] ) */

    // Fine Mesh Variables
    QuadratureBase * _quadratureFine;       /*! @brief Quadrature in fine mesh (used for uncollided part without local refinement) */
    unsigned _nqF;                          /*! @brief Nr of quadrature points in fine mesh */
    VectorVector _quadPointsF;              /*! @brief quadrature points of fine mesh */
    Vector _weightsF;                       /*! @brief weights of quadrature points of fine mesh */

    // coarse Mesh Variables
    unsigned _nqC;                          /*! @brief Nr of quadrature points in coarse mesh */
    VectorVector _quadPointsC;              /*! @brief Nr of quadrature points in coarse mesh */
    Vector _weightsC;                       /*! @brief weights of quadrature points of coarse mesh */
    VectorVectorU _neigh;                   /*! @brief Stores the neighbor dependency of quadrature points in coarse mesh */

    // Variables used for uncollided part --> new calculated for local refinement
    unsigned _nq;                           /*! @brief Nr of quadrature points used for uncollided part  (if NO local refinement: = _nqF)*/
    VectorVector _quadPoints;               /*! @brief quadrature points used for uncollided part (if NO local refinement: = _quadPointsF) */
    Vector _weights;                        /*! @brief weights of quadrature points used for uncollided part (if NO local refinement: = _weightsF) */

    // external Source
    std::vector<VectorVector> _QextC;       /*! @brief External Source calculated in coarse mesh */
    std::vector<VectorVector> _QextF;       /*! @brief External Source calculated in fine mesh */
    // _Q: external Source used in calculation (from solverbase)

    // local Refinement
    VectorU _refine;                        /*! @brief Vector that stores how often a quadPoint was refined, if refinement is needed*/
    VectorVectorU _quadC2F;                 /*! @brief VectorVector with fine quadrature points for each coarse quadrature point ( dim = [_nqC][nr of corresp. fine points]) */
    VectorU _quadIDs;                       /*! @brief VectorU with the IDs corresp to Coarse (ID<_nqC) or Fine (ID>=_nqC) quadarture mesh */
    blaze::CompressedMatrix<double> _solInter;       /*! @brief interpolation matrix for _solFC (dim = (_nq, _nqOld)) */

    unsigned _nqOld;                        /*! @brief Nr of Quadarture Points in the last iteration */
    VectorU _quadIDsOld;                    /*! @brief _quadIDs from last Itaertion */
};

#endif //FIRSTCOLLISIONCSDSNSOLVER_H
