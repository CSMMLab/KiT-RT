#ifndef MNSOLVER_H
#define MNSOLVER_H

#include "solverbase.h"

class EntropyBase;
class SphericalBase;
class OptimizerBase;

class MNSolver : public Solver
{
  public:
    /**
     * @brief MNSolver constructor
     * @param settings stores all needed information
     */
    MNSolver( Config* settings );

    /*! @brief MNSolver destructor */
    virtual ~MNSolver();

  private:
    // --- Private member variables ---
    unsigned _nTotalEntries;    /*! @brief: Total number of equations in the system */
    unsigned short _LMaxDegree; /*! @brief: Max Order of Spherical Harmonics  */

    // Moment basis
    SphericalBase* _basis; /*! @brief: Class to compute and store current spherical harmonics basis */
    VectorVector _moments; /*! @brief: Moment Vector pre-computed at each quadrature point: dim= _nq x _nTotalEntries */

    // Scattering
    Vector _scatterMatDiag; /*! @brief: Diagonal of the scattering matrix (its a diagonal matrix by construction) */

    // Quadrature related members
    VectorVector _quadPoints;       /*!  @brief quadrature points, dim(_quadPoints) = (_nq,spatialDim) */
    Vector _weights;                /*!  @brief quadrature weights, dim(_weights) = (_nq) */
    VectorVector _quadPointsSphere; /*!  @brief (my,phi), dim(_quadPoints) = (_nq,2) */

    // Entropy Optimization related members
    EntropyBase* _entropy;     /*! @brief: Class to handle entropy functionals */
    VectorVector _alpha;       /*! @brief: Lagrange Multipliers for Minimal Entropy problem for each gridCell
                                           Layout: _nCells x _nTotalEntries*/
    OptimizerBase* _optimizer; /*! @brief: Class to solve minimal entropy problem */

    VectorVector _solDx; /*! @brief:  temporary storage of x-derivatives of alpha */
    VectorVector _solDy; /*! @brief:  temporary storage of y-derivatives of alpha */

    // ---- Private Member functions ---

    // IO
    void PrepareVolumeOutput() override;
    void WriteVolumeOutput( unsigned idx_pseudoTime ) override;

    // Solver
    void FVMUpdate( unsigned idx_energy ) override;
    void FluxUpdate() override;
    void IterPreprocessing( unsigned idx_pseudotime ) override;
    void IterPostprocessing( unsigned idx_energy );
    /*! @brief : Construct flux by computing the Moment of the  sum of FVM discretization at the interface of cell
     *  @param : idx_cell = current cell id
     *  @returns : sum over all neighbors of flux for all moments at interface of idx_cell, idx_neighbor */
    Vector ConstructFlux( unsigned idx_cell );
    /*! @brief Corrects the solution _sol[idx_cell] to be realizable w.r.t. the reconstructed entropy (eta'(alpha*m))
        @param idx_cell = cell where the correction happens*/
    void ComputeRealizableSolution( unsigned idx_cell );

    // Initialization of the solver
    /*! @brief : Pre-Compute Moments at all quadrature points. */
    void ComputeMoments();
    /*! @brief:  fucntion for computing and setting up EV matrix for scattering kernel */
    void ComputeScatterMatrix();

    // Helper
    /*! @brief:  Computes the radiative flux from the solution vector of the moment system */
    void ComputeRadFlux();
};
#endif    // MNSOLVER_H
