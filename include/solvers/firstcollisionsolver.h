//
// Created by chinsp on 29/10/21.
//

#ifndef FIRSTCOLLISIONSOLVER_H
#define FIRSTCOLLISIONSOLVER_H

class EntropyBase;
class SphericalBase;
class OptimizerBase;

class FirstCollisionSolver : public SolverBase( settings )
{
protected:

    virtual void IterPreprocessing( unsigned idx_pseudotime ) override;
    virtual void IterPostprocessing( ) override;

    void FluxUpdate() override;
    virtual void FluxUpdateUncollided( );
    void FluxUpdateCollided( );

    void FVMUpdate( unsigned iter ) override;
    virtual void FVMUpdateUncollided( unsigned iter );
    void FVMUpdateCollided( unsigned iter );

    void ComputeMoments();
    virtual void ComputeFirstCollisionSource( unsigned iter );
    void ComputeRealizableSolution( unsigned idx_cell );
    virtual void ComputeRadFlux() override;

    Vector ConstructFluxMN( unsigned idx_cell );
    Vector ConstructFluxSN( unsigned idx_cell, unsigned nq );


public:
    // Constructor of FirstCollisonSolver
    FirstCollisionSolver( Config * settings );
    // Destructor
    virtual ~FirstCollisionSolver() {}


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


    // SN Variables
    unsigned _nqF;
    VectorVector _quadPointsF;
    Vector _weightsF;
    Matrix _scatteringKernelFC;


    // coarse Variables
    VectorVector _quadPoints;
    Vector _weights;
    VectorVector _quadPointsSphere;
    Matrix _scatteringKernel;
    VectorVector _QFirstCollision;


};

#endif //FIRSTCOLLISIONSOLVER_H
