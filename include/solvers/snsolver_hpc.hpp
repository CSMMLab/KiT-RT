#ifndef SNSOLVERHPC_H
#define SNSOLVERHPC_H

// include Matrix, Vector definitions
#include "common/globalconstants.hpp"
#include "common/typedef.hpp"

// externals
#include "spdlog/spdlog.h"

// Forward Declarations
class Config;

/*! @brief Base class for all solvers. */
class SNSolverHPC
{

  private:
    double _currTime;  /*!< @brief wall-time after current iteration */
    Config* _settings; /*!< @brief config class for global information */

    // Time
    unsigned _nIter; /*!< @brief number of time steps, for non CSD, this equals _nEnergies, for _csd, _maxIter =_nEnergies-1*/
    double _dT;      /*!< @brief energy/time step size */

    // Mesh related members, memory optimized
    unsigned _nCells;  /*!< @brief number of spatial cells */
    unsigned _nSystem; /*!< @brief number of equations in the transport system, i.e. num quad pts */
    unsigned _nq;      /*!< @brief number of quadrature points */
    unsigned _spatialDim;
    unsigned _nEdgesPerCell;

    std::vector<BOUNDARY_TYPE> _boundaryCells; /*!< @brief boundary type for all cells, dim(_boundary) =
                                                  (_NCells) */
    std::vector<double> _areas;                /*!< @brief surface area of all spatial cells,
                                                  dim(_areas) = _NCells */
    std::vector<double> _normals;              /*!< @brief edge normals multiplied by edge length,
                                                               dim(_normals) = (_NCells,nEdgesPerCell,spatialDim) */
    std::vector<unsigned> _neighbors;          /*!< @brief edge neighbor cell ids, dim(_neighbors) = (_NCells,nEdgesPerCell) */
    std::vector<double> _cellMidPoints;        /*!< @brief dim _nCells x _spatialDim */
    std::vector<double> _interfaceMidPoints;   /*!< @brief dim: _nCells x _nEdgesPerCell x _spatialDim */

    unsigned _temporalOrder;      /*!< @brief temporal order (current: 1 & 2) */
    unsigned _spatialOrder;       /*!< @brief spatial order (current: 1 & 2) */
    std::vector<double> _solDx;   /*!< @brief dim = _nCells x _nSys x _spatialDim*/
    std::vector<double> _limiter; /*!< @brief dim = _nCells x _nSys */

    // Scattering, absorption and source
    std::vector<double> _sigmaS;           /*!< @brief dim: _nCells */
    std::vector<double> _sigmaT;           /*!< @brief dim: _nCells */
    std::vector<double> _source;           /*!< @brief dim: _nCells x _nSys */
    std::vector<double> _scatteringKernel; /*!< @brief dim: _nSys x _nSys  */

    // quadrature related numbers
    std::vector<double> _quadPts;     /*!< @brief dim: _nSys x _spatialDim*/
    std::vector<double> _quadWeights; /*!< @brief dim: _nSys*/

    // Solution related members
    std::vector<double> _sol;    /*!< @brief dim = _nCells x _nSys */
    std::vector<double> _solNew; /*!< @brief dim = _nCells x _nSys */

    // Output related members
    std::vector<std::vector<std::vector<double>>> _outputFields; /*!< @brief Solver Output: dimensions
                   (GroupID,FieldID,CellID).*/
    std::vector<std::vector<std::string>> _outputFieldNames;     /*!< @brief Names of the outputFields: dimensions
                           (GroupID,FieldID) */

    std::vector<double> _screenOutputFields;          /*!< @brief Solver Output: dimensions (FieldID). */
    std::vector<std::string> _screenOutputFieldNames; /*!< @brief Names of the outputFields: dimensions (FieldID) */

    std::vector<double> _historyOutputFields;          /*!< @brief Solver Output: dimensions (FieldID). */
    std::vector<std::string> _historyOutputFieldNames; /*!< @brief Names of the outputFields: dimensions (FieldID) */

    // ---- Member functions ----

    // Solver
    /*! @brief Performs preprocessing steps before the pseudo time iteration is
     * started*/
    void SolverPreprocessing();
    /*! @brief Performs preprocessing for the current solver iteration
        @param idx_iter current (peudo) time iteration */
    void IterPreprocessing( unsigned idx_iter ) = 0;
    /*! @brief Performs postprocessing for the current solver iteration */
    void IterPostprocessing( unsigned idx_iter );
    /*! @brief Constructs  the flux update for the current iteration and stores it
     * in psiNew*/
    void FluxUpdate() = 0;
    /*! @brief Computes the finite Volume update step for the current iteration
         @param idx_iter  current (peudo) time iteration */
    void FVMUpdate( unsigned idx_iter ) = 0;
    /*! @brief Computes the finite Volume update step for the current iteration
         @param idx_iter  current (peudo) time iteration */
    void RKUpdate( std::vector<double>& sol0, std::vector<double> sol_rk );

    // Helper
    /*! @brief ComputeTimeStep calculates the maximal stable time step using the
       cfl number
        @param cfl Courant-Friedrichs-Levy condition number */
    double ComputeTimeStep( double cfl ) const;
    /*! @brief Computes the flux of the solution to check conservation properties
     */
    void ComputeScalarFlux() = 0;

    // IO
    /*! @brief Initializes the output groups and fields of this solver and names
     * the fields */
    void PrepareVolumeOutput() = 0;
    /*! @brief Function that prepares VTK export and csv export of the current
       solver iteration
        @param idx_iter  current (pseudo) time iteration */
    void WriteVolumeOutput( unsigned idx_iter ) = 0;
    /*! @brief Save Output solution at given energy (pseudo time) to VTK file.
       Write frequency is given by option VOLUME_OUTPUT_FREQUENCY. Always prints
       last iteration without iteration affix.
        @param idx_iter  current (pseudo) time iteration */
    void PrintVolumeOutput( int idx_iter ) const;
    /*! @brief Initialized the output fields and their Names for the screenoutput
     */
    void PrepareScreenOutput();

    /*! @brief Function that writes screen and history output fields
        @param idx_iter  current (pseudo) time iteration */
    void WriteScalarOutput( unsigned idx_iter );
    /*! @brief Prints ScreenOutputFields to Screen and to logger. Write frequency
       is given by option SCREEN_OUTPUT_FREQUENCY. Always prints last iteration.
        @param idx_iter  current (pseudo) time iteration */
    void PrintScreenOutput( unsigned idx_iter );
    /*! @brief Initialized the historyOutputFields and their Names for history
       output. Write frequency is given by option HISTORY_OUTPUT_FREQUENCY. Always
       prints last iteration. */
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
    SNSolverHPC( Config* settings );

    ~SNSolverHPC();

    /*! @brief Solve functions runs main iteration loop. Components of the solve
     * loop are pure  and subclassed by the child solvers.  */
    void Solve();

    /*! @brief Save Output solution to VTK file */
    void PrintVolumeOutput() const;    // Only for debugging purposes.
};
#endif    // SNSOLVERHPC_H
