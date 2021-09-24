#ifndef SOLVER_H
#define SOLVER_H

// include Matrix, Vector definitions
#include "common/globalconstants.h"
#include "common/typedef.h"

// externals
#include "spdlog/spdlog.h"

// Forward Declarations
class NumericalFlux;
class Mesh;
class Config;
class ProblemBase;
class QuadratureBase;
class Reconstructor;

/*! @brief Base class for all solvers. */
class SolverBase
{
  protected:
    NumericalFlux* _g;     /*!< @brief class for numerical flux */
    Config* _settings;     /*!< @brief config class for global information */
    ProblemBase* _problem; /*!< @brief problem class for initial conditions */

    // --------- Often used variables of member classes for faster access ----

    // Time or Energystepping
    unsigned _maxIter;   /*!< @brief number of time steps, for non CSD, this equals _nEnergies, for _csd, _maxIter = _nEnergies-1*/
    unsigned _nEnergies; /*!< @brief number of energysteps, number of nodal energy values for CSD */
    double _dE;          /*!< @brief energy/time step size */
    Vector _energies;    /*!< @brief energy groups used in the simulation [keV] */

    // Mesh related members
    Mesh* _mesh;                               /*!< @brief mesh object for writing out information */
    unsigned _nCells;                          /*!< @brief number of spatial cells */
    std::vector<BOUNDARY_TYPE> _boundaryCells; /*!< @brief boundary type for all cells, dim(_boundary) = (_NCells) */
    std::vector<double> _areas;                /*!< @brief surface area of all spatial cells, dim(_areas) = _NCells */
    std::vector<std::vector<Vector>>
        _normals; /*!< @brief edge normals multiplied by edge length, dim(_normals) = (_NCells,nEdgesPerCell,spatialDim) */
    std::vector<std::vector<unsigned>> _neighbors; /*!< @brief edge neighbor cell ids, dim(_neighbors) = (_NCells,nEdgesPerCell) */

    // slope related params
    Reconstructor* _reconstructor;                        /*!< @brief reconstructor object for high-order scheme */
    unsigned _reconsOrder;                                /*!< @brief reconstruction order (current: 1 & 2) */
    VectorVector _cellMidPoints;                          /*!< @brief middle point locations of elements */
    std::vector<std::vector<Vector>> _interfaceMidPoints; /*!< @brief middle point locations of edges */

    std::vector<double> _density; /*!< @brief patient density, dim(_density) = _nCells (only for csdsolver) */
    Vector _s;                    /*!< @brief stopping power, dim(_s) = _maxIter (only for csdsolver) */

    std::vector<VectorVector> _Q; /*!< @brief external source term. Dim(_Q) = _maxIter x (_nCells x _nSystem) */
    VectorVector _sigmaS;         /*!< @brief scattering cross section for all energies. len: _nEnergies x _nCells */
    VectorVector _sigmaT;         /*!< @brief total cross section for all energies.  len: _nEnergies x _nCells*/

    // quadrature related numbers
    QuadratureBase* _quadrature; /*!< @brief pointer to quadrature class */
    unsigned _nq;                /*!< @brief number of quadrature points */

    // Solution related members
    VectorVector _sol;     /*!< @brief solution of the PDE, e.g. angular flux or moments */
    VectorVector _solNew;  /*!< @brief VectorVector to store the new flux and later the new solution per iteration */
    VectorVector _solDx;   /*!< @brief  solution gradient for each cell and each system var in x direction */
    VectorVector _solDy;   /*!< @brief  solution gradient for each cell and each system var in x direction */
    VectorVector _limiter; /*! < @brief slope limiter at cell and system pt */
    Vector _fluxNew;       /*!< @brief Vector to store the new Flux. Dim _nCells */
    Vector _flux;          /*!< @brief Vector to store the old Flux. Dim _nCells*/

    std::vector<double> _solverOutput; /*!< @brief LEGACY: Outputfield for solver ==> Will be replaced by _outputFields in the near future */

    // Output related members
    std::vector<std::vector<std::vector<double>>> _outputFields; /*!< @brief Solver Output: dimensions (GroupID,FieldID,CellID).*/
    std::vector<std::vector<std::string>> _outputFieldNames;     /*!< @brief Names of the outputFields: dimensions (GroupID,FieldID) */
    // we will have to add a further dimension for quadPoints and weights once we start with multilevel SN

    // Output related members
    std::vector<double> _screenOutputFields;          /*!< @brief Solver Output: dimensions (FieldID). */
    std::vector<std::string> _screenOutputFieldNames; /*!< @brief Names of the outputFields: dimensions (FieldID) */

    // Output related members
    std::vector<double> _historyOutputFields;          /*!< @brief Solver Output: dimensions (FieldID). */
    std::vector<std::string> _historyOutputFieldNames; /*!< @brief Names of the outputFields: dimensions (FieldID) */

    // ---- Member functions ----

    // Solver
    /*! @brief Performs preprocessing steps before the pseudo time iteration is started*/
    virtual void SolverPreprocessing();
    /*! @brief Performs preprocessing for the current solver iteration
        @param idx_iter current (peudo) time iteration */
    virtual void IterPreprocessing( unsigned idx_iter ) = 0;
    /*! @brief Performs postprocessing for the current solver iteration */
    virtual void IterPostprocessing( unsigned idx_pseudotime ) = 0;
    /*! @brief Constructs  the flux update for the current iteration and stores it in psiNew*/
    virtual void FluxUpdate() = 0;
    /*! @brief Computes the finite Volume update step for the current iteration
         @param idx_iter  current (peudo) time iteration */
    virtual void FVMUpdate( unsigned idx_iter ) = 0;

    // High order spatial reconstruction
    void ComputeGradients( unsigned nSys );               /*!<  @brief Computes the cell gradient using green gauss theorem */
    void ComputeLimiter();                                /*!< @brief Compute the slope limiter for all cells */
    void ComputeInterfaceValue( const unsigned idx_nbr ); /*!< @brief Compute the slope values at the interface to neighbor idx_nbr */

    // Helper
    /*! @brief ComputeTimeStep calculates the maximal stable time step using the cfl number
        @param cfl Courant-Friedrichs-Levy condition number */
    double ComputeTimeStep( double cfl ) const;
    /*! @brief Computes the flux of the solution to check conservation properties */
    virtual void ComputeRadFlux() = 0;

    // IO
    /*! @brief Initializes the output groups and fields of this solver and names the fields */
    virtual void PrepareVolumeOutput() = 0;
    /*! @brief Function that prepares VTK export and csv export of the current solver iteration
        @param idx_iter  current (pseudo) time iteration */
    virtual void WriteVolumeOutput( unsigned idx_iter ) = 0;
    /*! @brief Save Output solution at given energy (pseudo time) to VTK file. Write frequency is given by
               option VOLUME_OUTPUT_FREQUENCY. Always prints last iteration without iteration affix.
        @param idx_iter  current (pseudo) time iteration */
    void PrintVolumeOutput( int idx_iter ) const;
    /*! @brief Initialized the output fields and their Names for the screenoutput */
    void PrepareScreenOutput();

    /*! @brief Function that writes screen and history output fields
        @param idx_iter  current (pseudo) time iteration */
    void WriteScalarOutput( unsigned idx_iter );
    /*! @brief Prints ScreenOutputFields to Screen and to logger. Write frequency is given by
               option SCREEN_OUTPUT_FREQUENCY. Always prints last iteration.
        @param idx_iter  current (pseudo) time iteration */
    void PrintScreenOutput( unsigned idx_iter );
    /*! @brief Initialized the historyOutputFields and their Names for history output. Write frequency is given by
               option HISTORY_OUTPUT_FREQUENCY. Always prints last iteration. */
    void PrepareHistoryOutput();
    /*! @brief Prints HistoryOutputFields to logger
        @param idx_iter  current (pseudo) time iteration */
    void PrintHistoryOutput( unsigned idx_iter );
    /*! @brief Pre Solver Screen and Logger Output */
    void DrawPreSolverOutput();
    /*! @brief Post Solver Screen and Logger Output */
    void DrawPostSolverOutput();

  public:
    /*! @brief Solver constructor
     *  @param settings config class that stores all needed config information */
    SolverBase( Config* settings );

    virtual ~SolverBase();

    /*! @brief Create constructor
     *  @param settings config class that stores all needed config information
     *  @return pointer to SolverBase */
    static SolverBase* Create( Config* settings );

    /*! @brief Solve functions runs main iteration loop. Components of the solve loop are pure virtual and subclassed by the child solvers.  */
    virtual void Solve();

    /*! @brief Save Output solution to VTK file */
    void PrintVolumeOutput() const;    // Only for debugging purposes.
};
#endif    // SOLVER_H
