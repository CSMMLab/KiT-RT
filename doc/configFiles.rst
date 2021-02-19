Configuration files
-----------------------

Configuration files are written in `TOML-style <https://github.com/toml-lang/toml>`_
. A number of examples (link examples) can be found in the ``./code/input`` folder. Users can also handover their own files, which should follow the following specifications:

********************
Parameters
********************

Below you can see the possible parameters that can be specified in a config file. Which of these parameters are required depends on the chosen application and solver



.. code-block:: c++

    // --- Options ---
    // File Structure
    std::string _inputDir;     // Directory for input files
    std::string _outputDir;    // Directory for output files
    std::string _outputFile;   // Name of output file
    std::string _logDir;       // Directory of log file
    std::string _logFileName;  // Name of log file
    std::string _meshFile;     // Name of mesh file
    std::string _ctFile;       // Name of CT file

    // Quadrature
    QUAD_NAME _quadName;       // Quadrature Name
    unsigned short _quadOrder; // Quadrature Order
    unsigned _nQuadPoints;

    // Mesh
    unsigned _nCells;

    // Solver
    double _CFL;                     // CFL Number for Solver
    double _tEnd;                    // Final Time for Simulation 
    PROBLEM_NAME _problemName;       // Name of predefined Problem   
    SOLVER_NAME _solverName;         // Name of the used Solver 
    ENTROPY_NAME _entropyName;       // Name of the used Entropy Functional 
    unsigned short _maxMomentDegree; // Maximal Order of Moments for PN and MN Solver 
    unsigned short _reconsOrder;     // Spatial Order of Accuracy for Solver 

    // Linesource
    double _sigmaS; // Scattering coeffient for Linesource test case 

    /*If true, very low entries (10^-10 or smaller) of the flux matrices will be set to zero,
     * to improve floating point accuracy */
    bool _cleanFluxMat;

    /*If true, the SN Solver uses all Gauss pts in the quadrature */
    bool _allGaussPts; 
    
    /*If true, continuous slowing down approximation will be used */
    bool _csd;                 

    std::string _hydrogenFile; // Name of hydrogen cross section file 
    std::string _oxygenFile;   // Name of oxygen cross section file 

    // Boundary Conditions
    /*List of all Pairs (marker, BOUNDARY_TYPE), e.g. (farfield,DIRICHLET).
     *Each Boundary Conditions must have an entry in enum BOUNDARY_TYPE*/
    std::vector<std::pair<std::string, BOUNDARY_TYPE>> _boundaries;

    unsigned short _nMarkerDirichlet;          // Number of Dirichlet BC markers. Enum entry: DIRICHLET 
    unsigned short _nMarkerNeumann;            // Number of Neumann BC markers. Enum entry: Neumann 
    std::vector<std::string> _MarkerDirichlet; // Dirichlet BC markers. 
    std::vector<std::string> _MarkerNeumann;   // Neumann BC markers. 

    // Scattering Kernel
    KERNEL_NAME _kernelName; // Scattering Kernel Name

    // Optimizer
    OPTIMIZER_NAME _entropyOptimizerName; // Choice of optimizer 
    double _optimizerEpsilon;             // termination criterion epsilon for Newton Optmizer 
    unsigned short _newtonIter;           // Maximal Number of newton iterations 
    double _newtonStepSize;               // Stepsize factor for newton optimizer 
    unsigned short _newtonLineSearchIter; // Maximal Number of line search iterations for newton optimizer 
    bool _newtonFastMode;                 // If true, we skip the NewtonOptimizer for quadratic entropy and assign alpha = u 

    // Output Options
    unsigned short _nVolumeOutput;            // Number of volume outputs 
    std::vector<VOLUME_OUTPUT> _volumeOutput; // Output groups for volume output
    unsigned short _volumeOutputFrequency;    // Frequency of vtk write of volume output

    unsigned short _nScreenOutput;            // Number of screen outputs 
    std::vector<SCALAR_OUTPUT> _screenOutput; // Output groups for screen output
    unsigned short _screenOutputFrequency;    // Frequency of screen output

    unsigned short _nHistoryOutput;            // Number of screen outputs 
    std::vector<SCALAR_OUTPUT> _historyOutput; // Output groups for screen output
    unsigned short _historyOutputFrequency;    // Frequency of screen output




