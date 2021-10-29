//
// Created by chinsp on 29/10/21.
//

#ifndef FIRSTCOLLISIONMNSOLVER_H
#define FIRSTCOLLISIONMNSOLVER_H

#include "solverbase.h"

class EntropyBase;
class SphericalBase;
class OptimizerBase;

class FirstCollisionMNSolver : public SolverBase( settings )
{
private:

    virtual void IterPreprocessing( unsigned idx_pseudotime ) override;
    virtual void IterPostprocessing( ) override;

    void FluxUpdate() override; // necessary to define but not used
    virtual void FluxUpdateUncollided( );
    void FluxUpdateCollided( );

    void FVMUpdate( unsigned iter ) override;   // necessary to define but  not used
    virtual void FVMUpdateUncollided( unsigned iter );
    void FVMUpdateCollided( unsigned iter );

    void ComputeMoments();
    virtual void ComputeFirstCollisionSource( unsigned iter );
    void ComputeRealizableSolution( unsigned idx_cell );
    virtual void ComputeRadFlux() override;

    Vector ConstructFluxMN( unsigned idx_cell );
    Vector ConstructFluxSN( unsigned idx_cell );


public:
    // Constructor of FirstCollisonSolver
    FirstCollisionMNSolver( Config * settings );
    // Destructor
    virtual ~FirstCollisionMNSolver() {}


    void Solve() override; // overrides Basis Solver !!

private:
    void PrepareVolumeOutput() override;
    void WriteVolumeOutput( unsigned idx_pseudoTime ) override;


    // Variables

protected:

    VectorVector _solF;
    VectorVector _solNewF;

    // MN Variables
    unsigned _LMaxDegree;
    SphericalBase * _basis;
    unsigned _nTotalEntries;
    EntropyBase * _entropy;
    OptimizerBase * _optimizer;
    VectorVector _alpha;
    VectorVector _moments;
    Vector _scatterMatDiag;

    // coarse Variables
    VectorVector _quadPoints;
    Vector _weights;
    VectorVector _quadPointsSphere;
    Matrix _scatteringKernel;
    VectorVector _QFirstCollision;




};


#endif //FIRSTCOLLISIONMNSOLVER_H
