/*!
 * @file config.h
 * @brief Class to handle all options and their pre and postprocessing.
 *         DO NOT CREATE SETTERS FOR THIS CLASS! ALL OPTIONS ARE CONSTANT (after SetPostprocessing).
 *
 * @author S. Schotthöfer
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <filesystem>
#include <map>
#include <vector>

#include "globalconstants.hpp"

// Forward declaration
class OptionBase;

/*!
 * @class Config
 * @brief Main class for defining the problem; basically this class reads the configuration file, and
 *        stores all the information.
 */

class Config
{
  private:
    std::string _fileName; /*!< @brief Name of the current file without extension */
    bool _baseConfig;

    // int _commRank, _commSize; /*!< @brief MPI rank and size.*/    // Not yet used!!

    // --- Options ---
    // File Structure
    std::string _inputDir;    /*!< @brief Directory for input files*/
    std::string _outputDir;   /*!< @brief Directory for output files*/
    std::string _outputFile;  /*!< @brief Name of output file*/
    std::string _logDir;      /*!< @brief Directory of log file*/
    std::string _logFileName; /*!< @brief Name of log file*/
    std::string _meshFile;    /*!< @brief Name of mesh file*/
    std::string _ctFile;      /*!< @brief Name of CT file*/

    // Quadrature
    QUAD_NAME _quadName;       /*!< @brief Quadrature Name*/
    unsigned short _quadOrder; /*!< @brief Quadrature Order*/
    unsigned _nQuadPoints;     /*!< @brief Number of quadrature points. (Deprecated)*/
    // std::vector<double> _1dIntegrationBounds; /*!< @brief Quadrature Order*/

    // Mesh
    unsigned _nCells;    /*!< @brief Number of cells in the mesh */
    unsigned short _dim; /*!< @brief spatial dimensionality of the mesh/test case */

    // Boundary Conditions
    /*!< @brief List of all Pairs (marker, BOUNDARY_TYPE), e.g. (farfield,DIRICHLET).
         Each Boundary Conditions must have an entry in enum BOUNDARY_TYPE*/
    std::vector<std::pair<std::string, BOUNDARY_TYPE>> _boundaries;
    unsigned short _nMarkerDirichlet;          /*!< @brief Number of Dirichlet BC markers. Enum entry: DIRICHLET */
    unsigned short _nMarkerNeumann;            /*!< @brief Number of Neumann BC markers. Enum entry: Neumann */
    std::vector<std::string> _MarkerDirichlet; /*!< @brief Dirichlet BC markers. */
    std::vector<std::string> _MarkerNeumann;   /*!< @brief Neumann BC markers. */

    // Solver
    double _CFL;                     /*!< @brief CFL Number for Solver*/
    double _tEnd;                    /*!< @brief Final Time for Simulation */
    PROBLEM_NAME _problemName;       /*!< @brief Name of predefined Problem   */
    SOLVER_NAME _solverName;         /*!< @brief Name of the used Solver */
    ENTROPY_NAME _entropyName;       /*!< @brief Name of the used Entropy Functional */
    unsigned short _maxMomentDegree; /*!< @brief Maximal Order of Moments for PN and MN Solver */
    unsigned short _reconsOrder;     /*!< @brief Spatial Order of Accuracy for Solver */
    bool _realizabilityRecons;       /*!< @brief Turns realizability reconstruction on/off for u sampling and MN solver */

    /*!< @brief If true, very low entries (10^-10 or smaller) of the flux matrices will be set to zero,
     * to improve floating point accuracy */
    bool _cleanFluxMat;
    bool _allGaussPts; /*!< @brief If true, the SN Solver uses all Gauss pts in the quadrature */
    bool _csd;         /*!< @brief If true, continuous slowing down approximation will be used */

    // --- Problems ---

    // Linesource
    double _sigmaS; /*!< @brief Scattering coeffient for Linesource test case */

    // Database ICRU
    std::string _dataDir; /*!< @brief material directory */
    // ElectronRT
    std::string _hydrogenFile;      /*!< @brief Name of hydrogen cross section file path*/
    std::string _oxygenFile;        /*!< @brief Name of oxygen cross section file path */
    std::string _stoppingPowerFile; /*!< @brief Name of stopping power file path */

    // CSD
    double _maxEnergyCSD; /*!< @brief Maximum energy for CSD simulation */

    // --- other variables ---
    // Scattering Kernel
    KERNEL_NAME _kernelName; /*!< @brief Scattering Kernel Name*/

    // Spherical Basis
    SPHERICAL_BASIS_NAME _sphericalBasisName; /*!< @brief Name of the basis on the unit sphere */

    // Optimizer
    OPTIMIZER_NAME _entropyOptimizerName; /*!< @brief Choice of optimizer */
    double _optimizerEpsilon;             /*!< @brief termination criterion epsilon for Newton Optmizer */
    unsigned long _newtonIter;            /*!< @brief Maximal Number of newton iterations */
    double _newtonStepSize;               /*!< @brief Stepsize factor for newton optimizer */
    unsigned long _newtonLineSearchIter;  /*!< @brief Maximal Number of line search iterations for newton optimizer */
    bool _newtonFastMode;                 /*!< @brief If true, we skip the NewtonOptimizer for quadratic entropy and assign alpha = u */
    double _regularizerGamma;             /*!< @brief Regularization parameter for the regularized closure */
    // NeuralModel
    unsigned short _neuralModel; /*!< @brief  Version number of the employed neural model */
    // Output Options
    unsigned short _nVolumeOutput;            /*!< @brief Number of volume outputs */
    std::vector<VOLUME_OUTPUT> _volumeOutput; /*!< @brief Output groups for volume output*/
    unsigned short _volumeOutputFrequency;    /*!< @brief Frequency of vtk write of volume output*/

    unsigned short _nScreenOutput;            /*!< @brief Number of screen outputs */
    std::vector<SCALAR_OUTPUT> _screenOutput; /*!< @brief Output groups for screen output*/
    unsigned short _screenOutputFrequency;    /*!< @brief Frequency of screen output*/

    unsigned short _nHistoryOutput;            /*!< @brief Number of screen outputs */
    std::vector<SCALAR_OUTPUT> _historyOutput; /*!< @brief Output groups for screen output*/
    unsigned short _historyOutputFrequency;    /*!< @brief Frequency of screen output*/

    // Data Generator Settings
    /*!< @brief Check, if data generator mode is active. If yes, no solver is called, but instead the data generator is executed */
    bool _dataGeneratorMode;
    SAMPLER_NAME _sampler;            /*!< @brief Sampling mode for regression or classification datasets */
    unsigned long _tainingSetSize;    /*!< @brief Size of training data set for data generator */
    bool _sizeByDimension;            /*!< @brief If true, the value of _trainingSetSize is the number of gridpoints in one dimension */
    unsigned long _maxValFirstMoment; /*!< @brief Size of training data set for data generator */
    double _RealizableSetEpsilonU0;   /*!< @brief Distance to 0 of the sampled moments to the boundary of the realizable set */
    double _RealizableSetEpsilonU1;   /*!< @brief norm(u_1)/u_0 !< _RealizableSetEpsilonU1 */
    bool _normalizedSampling;         /*!< @brief Flag for sampling of normalized moments, i.e. u_0 =1 */
    bool _alphaSampling;              /*!< @brief Flag for sampling alpha instead of u */
    double _alphaBound;               /*!< @brief The norm boundary for the sampling range of alpha*/
    double _minEVAlphaSampling;       /*!< @brief Rejection sampling criterion is a minimal eigenvalue threshold */
    bool _sampleUniform;              /*!< @brief If true, samples uniform, if false, sampleswith cutoff normal distribution */
    double _maxSamplingVelocity;      /*!< @brief The lower bound for the velocity space in the 1D classification sampler */
    // double _minSamplingVelocity;      /*!< @brief The upper bound for the velocity space in the 1D classification sampler */
    double _maxSamplingTemperature; /*!< @brief The lower bound for the interval to draw temperatures for the 1D classification sampler */
    double _minSamplingTemperature; /*!< @brief The upper bound for the interval to draw temperatures for the 1D classification sampler */
    unsigned short _nTemperatures;  /*!< @brief The number of sampling temperatures for the kinetic density sampler */

    // --- Parsing Functionality and Initializing of Options ---
    /*!
     * @brief Set default values for all options not yet set.
     */
    void SetDefault( void );

    /*!
     * @brief Set the config options.
     *        ==> Set new config options here.
     */
    void SetConfigOptions( void );

    /*!
     * @brief Set the config file parsing.
     */
    void SetConfigParsing( std::string case_filename );

    /*!
     * @brief Config file screen output.
     */
    void SetOutput( void );

    /*!
     * @brief Initializes pointers to null
     */
    void SetPointersNull( void );

    /*!
     * @brief Config file postprocessing.
     */
    void SetPostprocessing( void );

    /*!
     * @brief breaks an input line from the config file into a set of tokens
     * @param str the input line string
     * @param option_name the name of the option found at the beginning of the line
     * @param option_value the tokens found after the "=" sign on the line
     * @return false if the line is empty or a commment, true otherwise
     */
    bool TokenizeString( std::string& str, std::string& option_name, std::vector<std::string>& option_value );

    /*--- all_options is a map containing all of the options. This is used during config file parsing
     to track the options which have not been set (so the default values can be used). Without this map
     there would be no list of all the config file options. ---*/
    std::map<std::string, bool> _allOptions;

    /*--- brief param is a map from the option name (config file string) to its decoder (the specific child
     class of OptionBase that turns the string into a value) ---*/
    std::map<std::string, OptionBase*> _optionMap;

    // ---- Option Types ----

    // All of the addXxxOptions take in the name of the option, and a refernce to the field of that option
    // in the option structure. Depending on the specific type, it may take in a default value, and may
    // take in extra options. The addXxxOptions mostly follow the same pattern, so please see addDoubleOption
    // for detailed comments.
    //
    // List options are those that can be an unknown number of elements, and also take in a reference to
    // an integer. This integer will be populated with the number of elements of that type unmarshaled.
    //
    // Array options are those with a fixed number of elements.
    //
    // List and Array options should also be able to be specified with the string "NONE" indicating that there
    // are no elements. This allows the option to be present in a config file but left blank.

    /*!< @brief addDoubleOption creates a config file parser for an option with the given name whose
     value can be represented by a su2double.*/

    // Simple Options
    void AddBoolOption( const std::string name, bool& option_field, bool default_value );

    void AddDoubleOption( const std::string name, double& option_field, double default_value );

    void AddIntegerOption( const std::string name, int& option_field, int default_value );

    void AddLongOption( const std::string name, long& option_field, long default_value );

    void AddStringOption( const std::string name, std::string& option_field, std::string default_value );

    void AddUnsignedLongOption( const std::string name, unsigned long& option_field, unsigned long default_value );

    void AddUnsignedShortOption( const std::string name, unsigned short& option_field, unsigned short default_value );

    // enum types work differently than all of the others because there are a small number of valid
    // string entries for the type. One must also provide a list of all the valid strings of that type.
    template <class Tenum>
    void AddEnumOption( const std::string name, Tenum& option_field, const std::map<std::string, Tenum>& enum_map, Tenum default_value );

    // List Options
    void AddStringListOption( const std::string name, unsigned short& input_size, std::vector<std::string>& option_field );

    template <class Tenum>
    void AddEnumListOption( const std::string name,
                            unsigned short& num_marker,
                            std::vector<Tenum>& option_field,
                            const std::map<std::string, Tenum>& enum_map );

    // Initialize the cmdline and file logger
    void InitLogger();

  public:
    /*!
     * @brief Constructor of the class which reads the input file.
     */
    Config( std::string case_filename );

    /*!
     * @brief Destructor of the class.
     */
    ~Config( void );

    // ---- Getters for option values ----

    /*!
     * @brief Get Value of this option.
     *        Please keep alphabetical order within each subcategory
     */
    // File structure
    std::string inline GetCTFile() const { return std::filesystem::path( _ctFile ).lexically_normal(); }

    std::string inline GetLogDir() const { return std::filesystem::path( _logDir ).lexically_normal(); }
    std::string inline GetLogFile() const { return std::filesystem::path( _logFileName ).lexically_normal(); }
    std::string inline GetMeshFile() const { return std::filesystem::path( _meshFile ).lexically_normal(); }
    std::string inline GetOutputDir() const { return std::filesystem::path( _outputDir ).lexically_normal(); }
    std::string inline GetOutputFile() const { return std::filesystem::path( _outputFile ).lexically_normal(); }

    // Problem Files
    std::string inline GetHydrogenFile() const { return std::filesystem::path( _hydrogenFile ).lexically_normal(); }
    std::string inline GetOxygenFile() const { return std::filesystem::path( _oxygenFile ).lexically_normal(); }
    std::string inline GetStoppingPowerFile() const { return std::filesystem::path( _stoppingPowerFile ).lexically_normal(); }
    std::string inline GetDataDir() const { return std::filesystem::path( _dataDir ).lexically_normal(); }

    // Quadrature Structure
    unsigned GetNQuadPoints() { return _nQuadPoints; }
    QUAD_NAME inline GetQuadName() const { return _quadName; }
    unsigned short inline GetQuadOrder() const { return _quadOrder; }

    // Mesh Structure
    unsigned GetNCells() { return _nCells; }
    unsigned short GetDim() { return _dim; }

    // Solver Structure
    double inline GetCFL() const { return _CFL; }
    bool inline GetCleanFluxMat() const { return _cleanFluxMat; }
    ENTROPY_NAME inline GetEntropyName() const { return _entropyName; }
    unsigned short inline GetMaxMomentDegree() const { return _maxMomentDegree; }
    PROBLEM_NAME inline GetProblemName() const { return _problemName; }
    unsigned inline GetReconsOrder() { return _reconsOrder; }
    SOLVER_NAME inline GetSolverName() const { return _solverName; }
    double inline GetTEnd() const { return _tEnd; }
    bool inline GetSNAllGaussPts() const { return _allGaussPts; }
    bool inline GetIsCSD() const { return _csd; }
    bool inline GetRealizabilityReconstruction() { return _realizabilityRecons; }

    // Linesource
    double inline GetSigmaS() const { return _sigmaS; }

    // CSD
    double inline GetMaxEnergyCSD() const { return _maxEnergyCSD; }
    //  Optimizer
    double inline GetNewtonOptimizerEpsilon() const { return _optimizerEpsilon; }
    unsigned long inline GetNewtonIter() const { return _newtonIter; }
    double inline GetNewtonStepSize() const { return _newtonStepSize; }
    unsigned long inline GetNewtonMaxLineSearches() const { return _newtonLineSearchIter; }
    bool inline GetNewtonFastMode() const { return _newtonFastMode; }
    OPTIMIZER_NAME inline GetOptimizerName() const { return _entropyOptimizerName; }
    double inline GetRegularizerGamma() const { return _regularizerGamma; }
    // Neural Closure
    unsigned short inline GetNeuralModel() { return _neuralModel; }

    // Boundary Conditions
    BOUNDARY_TYPE GetBoundaryType( std::string nameMarker ) const; /*!< @brief Get Boundary Type of given marker */

    // Scattering Kernel
    KERNEL_NAME inline GetKernelName() const { return _kernelName; }

    // Basis name
    SPHERICAL_BASIS_NAME inline GetSphericalBasisName() const { return _sphericalBasisName; }
    // Output Structure
    std::vector<VOLUME_OUTPUT> inline GetVolumeOutput() { return _volumeOutput; }
    unsigned short inline GetNVolumeOutput() { return _nVolumeOutput; }
    unsigned short inline GetVolumeOutputFrequency() { return _volumeOutputFrequency; }

    std::vector<SCALAR_OUTPUT> inline GetScreenOutput() { return _screenOutput; }
    unsigned short inline GetNScreenOutput() { return _nScreenOutput; }
    unsigned short inline GetScreenOutputFrequency() { return _screenOutputFrequency; }

    std::vector<SCALAR_OUTPUT> inline GetHistoryOutput() { return _historyOutput; }
    unsigned short inline GetNHistoryOutput() { return _nHistoryOutput; }
    unsigned short inline GetHistoryOutputFrequency() { return _historyOutputFrequency; }

    // Data generator
    bool inline GetDataGeneratorMode() { return _dataGeneratorMode; }
    SAMPLER_NAME inline GetSamplerName() { return _sampler; }
    unsigned long inline GetTrainingDataSetSize() { return _tainingSetSize; }
    bool inline GetSizeByDimension() { return _sizeByDimension; }
    unsigned long inline GetMaxValFirstMoment() { return _maxValFirstMoment; }
    double GetRealizableSetEpsilonU0() { return _RealizableSetEpsilonU0; }
    double GetRealizableSetEpsilonU1() { return _RealizableSetEpsilonU1; }
    bool inline GetNormalizedSampling() { return _normalizedSampling; }
    bool inline GetAlphaSampling() { return _alphaSampling; }
    bool inline GetUniformSamlping() { return _sampleUniform; }
    double inline GetAlphaSamplingBound() { return _alphaBound; }
    double inline GetMinimalEVBound() { return _minEVAlphaSampling; }
    // double inline GetMinimalSamplingVelocity() { return _minSamplingVelocity; }
    double inline GetMaximalSamplingVelocity() { return _maxSamplingVelocity; }
    double inline GetMinimalSamplingTemperature() { return _minSamplingTemperature; }
    double inline GetMaximalSamplingTemperature() { return _maxSamplingTemperature; }
    unsigned short inline GetNSamplingTemperatures() { return _nTemperatures; }

    // ---- Setters for option structure
    // This section is dangerous
    // Quadrature Structure
    void SetNQuadPoints( unsigned nq ) { _nQuadPoints = nq; }        /*!< @brief Never change the nq! This is only for the test framework. */
    void SetQuadName( QUAD_NAME quadName ) { _quadName = quadName; } /*!< @brief Never change the quadName! This is only for the test framework. */
    void SetQuadOrder( unsigned quadOrder ) {
        _quadOrder = quadOrder;
    }                                                               /*!< @brief Never change the quadOrder! This is only for the test framework. */
    void SetSNAllGaussPts( bool useall ) { _allGaussPts = useall; } /*!< @brief Never change the this! This is only for the test framework. */
    // Mesh Structure
    void SetNCells( unsigned nCells ) { _nCells = nCells; }
};

#endif    // CONFIG_H
