/*!
 * \file CConfig.h
 * \brief Classes for different Options in rtsn
 * \author S. Schotthoefer
 *
 * Disclaimer: This class structure was copied and modifed with open source permission from SU2 v7.0.3 https://su2code.github.io/
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <map>
#include <mpi.h>

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_sinks.h"
#include "spdlog/spdlog.h"

#include "globalconstants.h"
#include "optionstructure.h"

/*!
 * \class Config
 * \brief Main class for defining the problem; basically this class reads the configuration file, and
 *        stores all the information.
 */

class Config
{
  private:
    std::string _fileName; /*!< \brief Name of the current file without extension */
    bool _baseConfig;

    int _commRank, _commSize; /*!< \brief MPI rank and size.*/    // Not yet used!!

    // --- Options ---
    // File Structure
    std::string _inputDir;   /*!< \brief Directory for input files*/
    std::string _outputDir;  /*!< \brief Directory for output files*/
    std::string _outputFile; /*!< \brief Name of output file*/
    std::string _logDir;     /*!< \brief Directory of log file*/
    std::string _meshFile;   /*!< \brief Name of mesh file*/

    // Quadrature
    QUAD_NAME _quadName;       /*!< \brief Quadrature Name*/
    unsigned short _quadOrder; /*!< \brief Quadrature Order*/
    unsigned _nQuadPoints;

    // Mesh
    unsigned _nCells;

    // Solver
    double _CFL;  /*!< \brief CFL Number for Solver*/
    double _tEnd; /*!< \brief Final Time for Simulation */
    PROBLEM_NAME _problemName;

    // Boundary Conditions
    /*!< \brief List of all Pairs (marker, BOUNDARY_TYPE), e.g. (farfield,DIRICHLET).
         Each Boundary Conditions must have an entry in enum BOUNDARY_TYPE*/
    std::vector<std::pair<std::string, BOUNDARY_TYPE>> _boundaries;
    unsigned short _nMarkerDirichlet;          /*!< \brief Number of Dirichlet BC markers. Enum entry: DIRICHLET */
    std::vector<std::string> _MarkerDirichlet; /*!< \brief Dirichlet BC markers. */

    // --- Parsing Functionality and Initializing of Options ---
    /*!
     * \brief Set default values for all options not yet set.
     */
    void SetDefault( void );

    /*!
     * \brief Set the config options.
     *        ==> Set new config options here.
     */
    void SetConfigOptions( void );

    /*!
     * \brief Set the config file parsing.
     */
    void SetConfigParsing( char case_filename[MAX_STRING_SIZE] );

    /*!
     * \brief Config file screen output.
     */
    void SetOutput( void );

    /*!
     * \brief Initializes pointers to null
     */
    void SetPointersNull( void );

    /*!
     * \brief Config file postprocessing.
     */
    void SetPostprocessing( void );

    /*!
     * \brief breaks an input line from the config file into a set of tokens
     * \param[in] str - the input line string
     * \param[out] option_name - the name of the option found at the beginning of the line
     * \param[out] option_value - the tokens found after the "=" sign on the line
     * \return false if the line is empty or a commment, true otherwise
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

    /*!< \brief addDoubleOption creates a config file parser for an option with the given name whose
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
    void InitLogger( spdlog::level::level_enum terminalLogLvl, spdlog::level::level_enum fileLogLvl );

  public:
    /*!
     * \brief Constructor of the class which reads the input file.
     */
    Config( char case_filename[MAX_STRING_SIZE] );

    /*!
     * \brief Destructor of the class.
     */
    ~Config( void );

    // ---- Getters for option values ----

    /*!
     * \brief Get Value of this option.
     *        Please keep alphabetical order within each subcategory
     */
    // File structure
    std::string inline GetMeshFile() const { return _meshFile; }
    std::string inline GetOutputDir() const { return _outputDir; }
    std::string inline GetOutputFile() const { return _outputFile; }
    std::string inline GetLogDir() const { return _logDir; }

    // Quadrature Structure
    QUAD_NAME inline GetQuadName() const { return _quadName; }
    unsigned short inline GetQuadOrder() const { return _quadOrder; }
    void SetNQuadPoints( unsigned nq ) { _nQuadPoints = nq; }
    unsigned GetNQuadPoints() { return _nQuadPoints; }

    // Mesh Structure
    void SetNCells( unsigned nCells ) { _nCells = nCells; }
    unsigned GetNCells() { return _nCells; }

    // Solver Structure
    double inline GetCFL() const { return _CFL; }
    double inline GetTEnd() const { return _tEnd; }
    PROBLEM_NAME inline GetProblemName() const { return _problemName; }

    // Boundary Conditions
    BOUNDARY_TYPE GetBoundaryType( std::string nameMarker ) const; /*! \brief Get Boundary Type of given marker */

    // Output Structure
};

#endif    // CONFIG_H
