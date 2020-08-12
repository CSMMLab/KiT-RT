/*!
 * @file CConfig.h
 * @brief Classes for different Options in rtsn
 * @author S. Schotthöfer
 *
 * Disclaimer: This class structure was copied and modifed with open source permission from SU2 v7.0.3 https://su2code.github.io/
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <filesystem>
#include <map>
#include <vector>

#include "globalconstants.h"

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
    std::string _fileName; /*!< @rief Name of the current file without extension */
    bool _baseConfig;

    // int _commRank, _commSize; /*!< @brief MPI rank and size.*/    // Not yet used!!

    // --- Options ---
    // File Structure
    std::string _inputDir;   /*!< @brief Directory for input files*/
    std::string _outputDir;  /*!< @brief Directory for output files*/
    std::string _outputFile; /*!< @brief Name of output file*/
    std::string _logDir;     /*!< @brief Directory of log file*/
    std::string _meshFile;   /*!< @brief Name of mesh file*/
    std::string _ctFile;     /*!< @brief Name of CT file*/

    // Quadrature
    QUAD_NAME _quadName;       /*!< @brief Quadrature Name*/
    unsigned short _quadOrder; /*!< @brief Quadrature Order*/
    unsigned _nQuadPoints;

    // Mesh
    unsigned _nCells;

    // Solver
    double _CFL;                     /*!< @brief CFL Number for Solver*/
    double _tEnd;                    /*!< @brief Final Time for Simulation */
    PROBLEM_NAME _problemName;       /*!< @brief Name of predefined Problem   */
    SOLVER_NAME _solverName;         /*!< @brief Name of the used Solver */
    ENTROPY_NAME _entropyName;       /*!< @brief Name of the used Entropy Functional */
    unsigned short _maxMomentDegree; /*!< @brief Maximal Order of Moments for PN and MN Solver */
    unsigned short _reconsOrder;     /*!< @brief Spatial Order of Accuracy for Solver */

    /*!< @brief If true, very low entries (10^-10 or smaller) of the flux matrices will be set to zero,
     * to improve floating point accuracy */
    bool _cleanFluxMat;

    bool _csd;                 /*!< @brief If true, continuous slowing down approximation will be used */
    std::string _hydrogenFile; /*!< @brief Name of hydrogen cross section file */
    std::string _oxygenFile;   /*!< @brief Name of oxygen cross section file */

    // Boundary Conditions
    /*!< @brief List of all Pairs (marker, BOUNDARY_TYPE), e.g. (farfield,DIRICHLET).
         Each Boundary Conditions must have an entry in enum BOUNDARY_TYPE*/
    std::vector<std::pair<std::string, BOUNDARY_TYPE>> _boundaries;
    unsigned short _nMarkerDirichlet;          /*!< @brief Number of Dirichlet BC markers. Enum entry: DIRICHLET */
    unsigned short _nMarkerNeumann;            /*!< @brief Number of Neumann BC markers. Enum entry: Neumann */
    std::vector<std::string> _MarkerDirichlet; /*!< @brief Dirichlet BC markers. */
    std::vector<std::string> _MarkerNeumann;   /*!< @brief Neumann BC markers. */

    // Scattering Kernel
    KERNEL_NAME _kernelName; /*!< @brief Scattering Kernel Name*/

    // Optimizer
    OPTIMIZER_NAME _entropyOptimizerName; /*!< @brief Choice of optimizer */
    double _optimizerEpsilon;             /*!< @brief termination criterion epsilon for Newton Optmizer */
    unsigned short _newtonIter;           /*!< @brief Maximal Number of newton iterations */
    double _newtonStepSize;               /*!< @brief Stepsize factor for newton optimizer */
    unsigned short _newtonLineSearchIter; /*!< @brief Maximal Number of line search iterations for newton optimizer */
    bool _newtonFastMode;                 /*!< @brief If true, we skip the NewtonOptimizer for quadratic entropy and assign alpha = u */

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
     * @param[in] str - the input line string
     * @param[out] option_name - the name of the option found at the beginning of the line
     * @param[out] option_value - the tokens found after the "=" sign on the line
     * @return false if the line is empty or a commment, true otherwise
     */
    bool TokenizeString( std::string& str, std::string& option_name, std::vector<std::string>& option_value );

    /*--- all_options is a map containing all of the options. This is used during config file parsing
     to track the options which have not been set (so the default values can be used). Without this map
     there would be no list of all the config file options. ---*/

    std::map<std::string, bool> _allOptions;

    /*--- brief param is a map from the option name (config file string) to its decoder (the specific child
     class of COptionBase that turns the string into a value) ---*/

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
    void AddStringListOption( const std::string name, unsigned short& num_marker, std::vector<std::string>& option_field );

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
    std::string inline GetMeshFile() const { return std::filesystem::path( _meshFile ).lexically_normal(); }
    std::string inline GetOutputDir() const { return std::filesystem::path( _outputDir ).lexically_normal(); }
    std::string inline GetOutputFile() const { return std::filesystem::path( _outputFile ).lexically_normal(); }
    std::string inline GetLogDir() const { return std::filesystem::path( _logDir ).lexically_normal(); }
    std::string inline GetCTFile() const { return std::filesystem::path( _ctFile ).lexically_normal(); }
    std::string inline GetHydrogenFile() const { return std::filesystem::path( _hydrogenFile ).lexically_normal(); }
    std::string inline GetOxygenFile() const { return std::filesystem::path( _oxygenFile ).lexically_normal(); }

    // Quadrature Structure
    QUAD_NAME inline GetQuadName() const { return _quadName; }
    unsigned short inline GetQuadOrder() const { return _quadOrder; }
    void SetNQuadPoints( unsigned nq ) { _nQuadPoints = nq; }
    unsigned GetNQuadPoints() { return _nQuadPoints; }

    // Mesh Structure
    void SetNCells( unsigned nCells ) { _nCells = nCells; }
    unsigned GetNCells() { return _nCells; }

    // Solver Structure
    unsigned short inline GetMaxMomentDegree() const { return _maxMomentDegree; }
    double inline GetCFL() const { return _CFL; }
    double inline GetTEnd() const { return _tEnd; }
    PROBLEM_NAME inline GetProblemName() const { return _problemName; }
    SOLVER_NAME inline GetSolverName() const { return _solverName; }
    ENTROPY_NAME inline GetEntropyName() const { return _entropyName; }
    bool inline GetCleanFluxMat() const { return _cleanFluxMat; }
    unsigned GetReconsOrder() { return _reconsOrder; }
    bool inline IsCSD() const { return _csd; }
    unsigned inline GetMaxMomentDegree() { return _maxMomentDegree; }

    //  Optimizer
    OPTIMIZER_NAME inline GetOptimizerName() const { return _entropyOptimizerName; }
    double inline GetNewtonOptimizerEpsilon() const { return _optimizerEpsilon; }
    unsigned inline GetNewtonIter() const { return _newtonIter; }
    double inline GetNewtonStepSize() const { return _newtonStepSize; }
    unsigned inline GetMaxLineSearches() const { return _newtonLineSearchIter; }
    bool inline GetNewtonFastMode() const { return _newtonFastMode; }

    // Boundary Conditions
    BOUNDARY_TYPE GetBoundaryType( std::string nameMarker ) const; /*! @brief Get Boundary Type of given marker */

    // Scattering Kernel
    KERNEL_NAME inline GetKernelName() const { return _kernelName; }

    // Output Structure
};

#endif    // CONFIG_H
