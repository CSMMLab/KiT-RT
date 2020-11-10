#ifndef MNSOLVER_H
#define MNSOLVER_H

#include "solverbase.h"

class EntropyBase;
class SphericalHarmonics;
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
    ~MNSolver();

    void Solve() override; /*! @brief Solve functions runs main time loop */

  private:
    unsigned _nTotalEntries;    /*! @brief: Total number of equations in the system */
    unsigned short _LMaxDegree; /*! @brief: Max Order of Moments */

    // Moment basis
    SphericalHarmonics* _basis; /*! @brief: Class to compute and store current spherical harmonics basis */
    VectorVector _moments;      /*! @brief: Moment Vector pre-computed at each quadrature point: dim= _nq x _nTotalEntries */

    // Right hand side members
    Vector _scatterMatDiag; /*! @brief: Diagonal of the scattering matrix (its a diagonal matrix by construction) */
    // TODO: Source

    // quadrature related members
    VectorVector _quadPoints;       /*!  @brief quadrature points, dim(_quadPoints) = (_nq,spatialDim) */
    Vector _weights;                /*!  @brief quadrature weights, dim(_weights) = (_nq) */
    VectorVector _quadPointsSphere; /*!  @brief (my,phi), dim(_quadPoints) = (_nq,2) */

    // Entropy Optimization related members
    EntropyBase* _entropy;     /*! @brief: Class to handle entropy functionals */
    VectorVector _alpha;       /*! @brief: Lagrange Multipliers for Minimal Entropy problem for each gridCell
                                           Layout: _nCells x _nTotalEntries*/
    OptimizerBase* _optimizer; /*! @brief: Class to solve minimal entropy problem */

    // ---- Member functions ---

    // IO
    /*! @brief Initializes the output groups and fields of this solver and names the fields */
    void PrepareOutputFields() override;

    /*! @brief Function that writes NN Training Data in a .csv file */
    void WriteNNTrainingData( unsigned idx_pseudoTime );

    /*! @brief Function that prepares VTK export and csv export of the current solver iteration
        @returns: Mass of current iteration
    */
    double WriteOutputFields( unsigned idx_pseudoTime ) override;

    /*! @brief : computes the global index of the moment corresponding to basis function (l,k)
     *  @param : degree l, it must hold: 0 <= l <=_nq
     *  @param : order k, it must hold: -l <=k <= l
     *  @returns : global index
     */
    int GlobalIndex( int l, int k ) const;

    /*! @brief : Construct flux by computing the Moment of the  sum of FVM discretization at the interface of cell
     *  @param : idx_cell = current cell id
     *  @returns : sum over all neighbors of flux for all moments at interface of idx_cell, idx_neighbor */
    Vector ConstructFlux( unsigned idx_cell );

    /*! @brief : Pre-Compute Moments at all quadrature points. */
    void ComputeMoments();

    /*! @brief:  fucntion for computing and setting up EV matrix for scattering kernel */
    void ComputeScatterMatrix();

    /*! @brief Corrects the solution _sol[idx_cell] to be realizable w.r.t. the reconstructed entropy (eta'(alpha*m))
        @param idx_cell = cell where the correction happens*/
    void ComputeRealizableSolution( unsigned idx_cell );

    // Solver
    void FVMUpdate( VectorVector& psiNew, unsigned idx_energy ) override;
    void FluxUpdate( VectorVector& psiNew ) override;
    void IterPreprocessing() override;
};
#endif    // MNSOLVER_H
