#ifndef SNSOLVERHPC_H
#define SNSOLVERHPC_H

// include Matrix, Vector definitions
#include "common/globalconstants.hpp"
#include "common/typedef.hpp"

// externals
#include "spdlog/spdlog.h"

// Forward Declarations
class Config;
class Mesh;
class ProblemBase;

/*! @brief Base class for all solvers. */
class SNSolverHPC
{

  private:
    int _rank;
    int _numProcs;
    unsigned long _localNSys;
    unsigned long _startSysIdx;
    unsigned long _endSysIdx;

    double _currTime;  /*!< @brief wall-time after current iteration */
    Config* _settings; /*!< @brief config class for global information */
    Mesh* _mesh;
    ProblemBase* _problem;

    // Time
    unsigned long _nIter; /*!< @brief number of time steps, for non CSD, this equals _nEnergies, for _csd, _maxIter =_nEnergies-1*/
    double _dT;           /*!< @brief energy/time step size */
    int _idx_start_iter;  /*!< @brief index of first iteration */
    // Mesh related members, memory optimized
    unsigned long _nCells; /*!< @brief number of spatial cells */
    unsigned long _nSys;   /*!< @brief number of equations in the transport system, i.e. num quad pts */
    unsigned long _nq;     /*!< @brief number of quadrature points */
    unsigned long _nDim;
    unsigned long _nNbr;
    unsigned long _nNodes;

    std::vector<double> _areas;         /*!< @brief surface area of all spatial cells,
                                           dim(_areas) = _NCells */
    std::vector<double> _normals;       /*!< @brief edge normals multiplied by edge length,
                                                        dim(_normals) = (_NCells,nEdgesPerCell,spatialDim) */
    std::vector<unsigned> _neighbors;   /*!< @brief edge neighbor cell ids, dim(_neighbors) = (_NCells,nEdgesPerCell) */
    std::vector<double> _cellMidPoints; /*!< @brief dim _nCells x _dim */

    std::vector<double> _interfaceMidPoints;             /*!< @brief dim: _nCells x _nEdgesPerCell x _dim */
    std::vector<BOUNDARY_TYPE> _cellBoundaryTypes;       /*!< @brief dim: _nCells x _nEdgesPerCell x _dim */
    std::map<unsigned, std::vector<double>> _ghostCells; /*!< @brief Vector of ghost cells for boundary conditions. CAN BE MORE EFFICIENT */
    std::vector<double> _relativeInterfaceMidPt;         /*!< @brief dim _nCells * _nNbr * _nDim */
    std::vector<double> _relativeCellVertices;           /*!< @brief dim _nCells * _nNbr * _nDim */

    std::map<unsigned, bool> _ghostCellsReflectingY; /*!< map that indicates if a ghostcell has a fixed value or is a mirroring boundary */
    std::map<unsigned, bool> _ghostCellsReflectingX; /*!< map that indicates if a ghostcell has a fixed value or is a mirroring boundary */
    std::vector<unsigned> _quadratureYReflection;    /*!< map that gives a Reflection against the y axis for the velocity ordinates */
    std::vector<unsigned> _quadratureXReflection;    /*!< map that gives a Reflection against the y axis for the velocity ordinates */

    unsigned _temporalOrder;      /*!< @brief temporal order (current: 1 & 2) */
    unsigned _spatialOrder;       /*!< @brief spatial order (current: 1 & 2) */
    std::vector<double> _solDx;   /*!< @brief dim = _nCells x _nSys x _dim*/
    std::vector<double> _limiter; /*!< @brief dim = _nCells x _nSys */

    // Scattering, absorption and source
    std::vector<double> _sigmaS;           /*!< @brief dim: _nCells */
    std::vector<double> _sigmaT;           /*!< @brief dim: _nCells */
    std::vector<double> _source;           /*!< @brief dim: _nCells x _nSys */
    std::vector<double> _scatteringKernel; /*!< @brief dim: _nSys x _nSys  */

    // quadrature related numbers
    std::vector<double> _quadPts;     /*!< @brief dim: _nSys x _dim*/
    std::vector<double> _quadWeights; /*!< @brief dim: _nSys*/

    // Solution related members
    std::vector<double> _sol;  /*!< @brief dim = _nCells x _nSys */
    std::vector<double> _flux; /*!< @brief dim = _nCells x _nSys */

    // Output related members
    std::vector<double> _scalarFlux; /*!< @brief dim = _nCells  */

    // Lattice QOIS
    unsigned _nOutputMoments;

    double _mass;
    double _rmsFlux;
    double _curAbsorptionLattice;    /*!< @brief Absorption of particles at Lattice checkerboard regions at current time step */
    double _totalAbsorptionLattice;  /*!< @brief Absorption of particles at Lattice checkerboard regions integrated until current time step */
    double _curMaxAbsorptionLattice; /*!< @brief Maximum pointwise absorption of particles at Lattice checkerboard regions  until current time step */
    double _curScalarOutflow;        /*!< @brief Outflow over whole boundary at current time step */
    double _totalScalarOutflow;      /*!< @brief Outflow over whole boundary integrated until current time step */
    double _curMaxOrdinateOutflow;   /*!< @brief Maximum ordinate-wise ouftlow  over boundary over all time steps */
    std::vector<double> _localMaxOrdinateOutflow; /*!< @brief Maximum ordinate-wise ouftlow  over boundary over all time steps */
    double _curScalarOutflowPeri1;                /*!< @brief Outflow over whole boundary at current time step */
    double _totalScalarOutflowPeri1;              /*!< @brief Outflow over whole boundary integrated until current time step */
    double _curScalarOutflowPeri2;                /*!< @brief Outflow over whole boundary at current time step */
    double _totalScalarOutflowPeri2;              /*!< @brief Outflow over whole boundary integrated until current time step */
    // helper
    std::map<unsigned, std::vector<unsigned>> _cellsLatticePerimeter1;
    std::map<unsigned, std::vector<unsigned>> _cellsLatticePerimeter2;
    std::vector<bool> _isPerimeterLatticeCell1;
    std::vector<bool> _isPerimeterLatticeCell2;

    // Hohlraum QOIS
    double _totalAbsorptionHohlraumCenter;
    double _totalAbsorptionHohlraumVertical;
    double _totalAbsorptionHohlraumHorizontal;
    double _curAbsorptionHohlraumCenter;
    double _curAbsorptionHohlraumVertical;
    double _curAbsorptionHohlraumHorizontal;
    double _varAbsorptionHohlraumGreen;

    std::vector<std::vector<unsigned>> _probingCellsHohlraum; /*!< @brief Indices of cells that contain a probing sensor */
    std::vector<double> _probingMoments;                      /*!< @brief Solution Momnets at the probing cells that contain a probing sensor */
    unsigned _probingMomentsTimeIntervals;                    /*!< @brief Solution Momnets at the probing cells that contain a probing sensor */

    unsigned _nProbingCellsLineGreen;               /*!< @brief Number of sampling cells that contain a probing sensor for the sliding window */
    std::vector<unsigned> _probingCellsLineGreen;   /*!< @brief Indices of cells that contain a probing sensor for the sliding window */
    std::vector<double> _absorptionValsLineSegment; /*!< @brief Avg Absorption value at the sampleing points of lineGreen */

    unsigned _nProbingCellsBlocksGreen;
    std::vector<std::vector<unsigned>> _probingCellsBlocksGreen; /*!< @brief Indices of cells that contain a probing sensor blocks */
    std::vector<double> _absorptionValsBlocksGreen;              /*!< @brief Avg Absorption value at the sampleing blocks of lineGreen */

    // Design parameters
    std::vector<double> _cornerUpperLeftGreen; /*!< @brief Coord of corner of the green area (minus thickness/2 of it) relative to the green center */
    std::vector<double> _cornerLowerLeftGreen; /*!< @brief Coord of corner of the green area (minus thickness/2 of it) relative to the green center */
    std::vector<double>
        _cornerUpperRightGreen; /*!< @brief Coord of corner of the green area (minus thickness/2 of it) relative to the green center */
    std::vector<double>
        _cornerLowerRightGreen; /*!< @brief Coord of corner of the green area (minus thickness/2 of it) relative to the green center */

    double _thicknessGreen;           /*!< @brief thickness of the green area */
    std::vector<double> _centerGreen; /*!< @brief Center of the Hohlraum   */

    // Output
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
    void FluxOrder1();
    void FluxOrder2();

    void FVMUpdate();

    void IterPostprocessing();

    void SetGhostCells(); /*!< @brief Sets vector of ghost cells for

    // Helper
    /*! @brief ComputeTimeStep calculates the maximal stable time step using the
    cfl number
    @param cfl Courant-Friedrichs-Levy condition number */
    double ComputeTimeStep( double cfl ) const;

    // IO
    /*! @brief Initializes the output groups and fields of this solver and names
     * the fields */
    void PrepareVolumeOutput();
    /*! @brief Function that prepares VTK export and csv export of the current
       solver iteration
        @param idx_iter  current (pseudo) time iteration */
    void WriteVolumeOutput( unsigned idx_iter );
    /*! @brief Save Output solution at given energy (pseudo time) to VTK file.
       Write frequency is given by option VOLUME_OUTPUT_FREQUENCY. Always prints
       last iteration without iteration affix.
        @param idx_iter  current (pseudo) time iteration */
    void PrintVolumeOutput( int idx_iter );
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

    // Helper
    unsigned long Idx2D( unsigned long idx1, unsigned long idx2, unsigned long len2 );
    unsigned long Idx3D( unsigned long idx1, unsigned long idx2, unsigned long idx3, unsigned long len2, unsigned long len3 );
    bool IsAbsorptionLattice( double x, double y ) const;
    void ComputeCellsPerimeterLattice();

    void SetProbingCellsLineGreen();
    void ComputeQOIsGreenProbingLine();
    std::vector<unsigned> linspace2D( const std::vector<double>& start, const std::vector<double>& end, unsigned num_points );

  public:
    /*! @brief Solver constructor
     *  @param settings config class that stores all needed config information */
    SNSolverHPC( Config* settings );

    ~SNSolverHPC();

    /*! @brief Solve functions runs main iteration loop. Components of the solve
     * loop are pure  and subclassed by the child solvers.  */
    void Solve();

    /*! @brief Save Output solution to VTK file */
    void PrintVolumeOutput() const {};    // Only for debugging purposes.
};
#endif    // SNSOLVERHPC_H
