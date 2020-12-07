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

class Solver
{
  protected:
    Mesh* _mesh;           /*! @brief mesh object for writing out information */
    NumericalFlux* _g;     /*! @brief class for numerical flux */
    Config* _settings;     /*! @brief config class for global information */
    ProblemBase* _problem; /*! @brief problem class for initial conditions */

    // --------- Often used variables of member classes for faster access ----

    unsigned _nEnergies;             /*! @brief number of energy/time steps, number of nodal energy values for CSD */
    double _dE;                      /*! @brief energy/time step size */
    Vector _energies;                // energy groups used in the simulation [keV]
    std::vector<double> _density;    // patient density, dim(_density) = _nCells
    Vector _s;                       // stopping power, dim(_s) = _nTimeSteps
    std::vector<VectorVector> _Q;    /*!  @brief  external source term */

    VectorVector _sigmaS; /*!  @brief scattering cross section for all energies */
    VectorVector _sigmaT; /*!  @brief total cross section for all energies */

    // quadrature related numbers
    QuadratureBase* _quadrature; /*! @brief quadrature to create members below */
    unsigned _nq;                /*! @brief number of quadrature points */

    // VectorVector _quadPoints;    /*!  @brief quadrature points, dim(_quadPoints) = (_nSystem,spatialDim) */
    // Vector _weights;             /*!  @brief quadrature weights, dim(_weights) = (_NCells) */

    // Mesh related members
    unsigned _nCells;                          /*! @brief number of spatial cells */
    std::vector<BOUNDARY_TYPE> _boundaryCells; /*! boundary type for all cells, dim(_boundary) = (_NCells) */
    std::vector<double> _areas;                /*! @brief surface area of all spatial cells, dim(_areas) = _NCells */
    /*! @brief edge normals multiplied by edge length, dim(_normals) = (_NCells,nEdgesPerCell,spatialDim) */
    std::vector<std::vector<Vector>> _normals;
    /*! @brief edge neighbor cell ids, dim(_neighbors) = (_NCells,nEdgesPerCell) */
    std::vector<std::vector<unsigned>> _neighbors;

    // slope related params
    Reconstructor* _reconstructor; /*! @brief reconstructor class for high order scheme */
    unsigned _reconsOrder;
    VectorVector _psiDx;
    VectorVector _psiDy;
    VectorVector _cellMidPoints;
    std::vector<std::vector<Vector>> _interfaceMidPoints;

    // Solution related members
    VectorVector _sol;                 /*! @brief solution of the PDE, e.g. angular flux or moments */
    std::vector<double> _solverOutput; /*! @brief LEGACY: Outputfield for solver ==> Will be replaced by _outputFields in the near future */

    // Output related members
    std::vector<std::vector<std::vector<double>>> _outputFields; /*! @brief: Solver Output: dimensions (GroupID,FieldID,CellID).*/
    std::vector<std::vector<std::string>> _outputFieldNames;     /*! @brief: Names of the outputFields: dimensions (GroupID,FieldID) */
    // we will have to add a further dimension for quadPoints and weights once we start with multilevel SN

    // Output related members
    std::vector<double> _screenOutputFields;          /*! @brief: Solver Output: dimensions (FieldID). */
    std::vector<std::string> _screenOutputFieldNames; /*! @brief: Names of the outputFields: dimensions (FieldID) */

    // Output related members
    std::vector<double> _historyOutputFields;          /*! @brief: Solver Output: dimensions (FieldID). */
    std::vector<std::string> _historyOutputFieldNames; /*! @brief: Names of the outputFields: dimensions (FieldID) */

    // Internal Members
    VectorVector _solNew; /*! @brief: VectorVector to store the new flux and later the new solution per iteration */    // REPLACES psiNEW
    Vector _fluxNew; /*! @brief: Vector to store the new Flux. Dim _nCells */
    Vector _flux;    /*! @brief: Vector to store the old Flux. Dim _nCells*/

    // ---- Member functions ----

    // Solver
    /*! @brief Performs preprocessing for the current solver iteration */
    virtual void IterPreprocessing( unsigned idx_pseudotime ) = 0;
    /*! @brief Performs postprocessing for the current solver iteration */
    virtual void IterPostprocessing() = 0;
    /*! @brief Constructs  the flux update for the current iteration and stores it in psiNew*/
    virtual void FluxUpdate() = 0;
    /*! @brief Computes the finite Volume update step for the current iteration */
    virtual void FVMUpdate( unsigned idx_energy ) = 0;

    // Helper
    /*! @brief ComputeTimeStep calculates the maximal stable time step */
    double ComputeTimeStep( double cfl ) const;
    /*! @brief: Computes the flux of the solution to check conservation properties */
    virtual void ComputeRadFlux() = 0;

    // IO
    /*! @brief Initializes the output groups and fields of this solver and names the fields */
    virtual void PrepareVolumeOutput() = 0;
    /*! @brief Function that prepares VTK export and csv export of the current solver iteration  */
    virtual void WriteVolumeOutput( unsigned iteration ) = 0;
    /*! @brief Save Output solution at given energy (pseudo time) to VTK file. Write frequency is given by
               option VOLUME_OUTPUT_FREQUENCY. Always prints last iteration without iteration affix.*/
    void PrintVolumeOutput( int currEnergy ) const;
    /*! @brief: Initialized the output fields and their Names for the screenoutput */
    void PrepareScreenOutput();
    /*! @brief Function that writes screen and history output fields */
    void WriteScalarOutput( unsigned iteration );
    /*! @brief Prints ScreenOutputFields to Screen and to logger. Write frequency is given by
               option SCREEN_OUTPUT_FREQUENCY. Always prints last iteration. */
    void PrintScreenOutput( unsigned iteration );
    /*! @brief: Initialized the historyOutputFields and their Names for history output. Write frequency is given by
               option HISTORY_OUTPUT_FREQUENCY. Always prints last iteration. */
    void PrepareHistoryOutput();
    /*! @brief Prints HistoryOutputFields to logger */
    void PrintHistoryOutput( unsigned iteration );

  public:
    /*! @brief Solver constructor
     *  @param settings stores all needed information */
    Solver( Config* settings );

    ~Solver();

    /*! @brief Create constructor
     *  @param settings stores all needed information
     *  @return pointer to Solver */
    static Solver* Create( Config* settings );

    /*! @brief Solve functions runs main time loop */
    virtual void Solve();

    /*! @brief Save Output solution to VTK file */
    void PrintVolumeOutput() const;    // Only for debugging purposes.
};
#endif    // SOLVER_H
