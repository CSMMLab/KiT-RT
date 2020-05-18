/*!
 * \file CConfig.h
 * \brief Classes for different Options in rtsn
 * \author S. Schotthoefer
 *
 * Disclaimer: This class structure was copied and modifed with open source permission from SU2 v7.0.3 https://su2code.github.io/
 */

#ifndef CONFIG_H
#define CONFIG_H

#include "GlobalConstants.h"
#include "OptionStructure.h"
#include <filesystem>
#include <map>

/*!
 * \class CConfig
 * \brief Main class for defining the problem; basically this class reads the configuration file, and
 *        stores all the information.
 */

class CConfig
{
  private:
    std::string _fileName; /*!< \brief Name of the current file without extension */
    bool _base_config;

    int _comm_rank, _comm_size; /*!< \brief MPI rank and size.*/

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
    // Solver
    double _CFL;  /*!< \brief CFL Number for Solver*/
    double _tEnd; /*!< \brief Final Time for Simulation */
    // Boundary Conditions
    /*!< \brief List of all Pairs (marker, BOUNDARY_TYPE), e.g. (farfield,DIRICHLET).
         Each Boundary Conditions must have an entry in enum BOUNDARY_TYPE*/
    std::vector<std::pair<std::string, BOUNDARY_TYPE>> _boundaries;
    unsigned short _nMarkerDirichlet; /*!< \brief Number of Dirichlet BC markers. Enum entry: DIRICHLET */
    std::string* _MarkerDirichlet;    /*!< \brief Dirichlet BC markers. */

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

    std::map<std::string, bool> all_options;

    /*--- brief param is a map from the option name (config file string) to its decoder (the specific child
     class of COptionBase that turns the string into a value) ---*/

    std::map<std::string, COptionBase*> option_map;

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
    void addBoolOption( const std::string name, bool& option_field, bool default_value );

    void addDoubleOption( const std::string name, double& option_field, double default_value );

    void addIntegerOption( const std::string name, int& option_field, int default_value );

    void addLongOption( const std::string name, long& option_field, long default_value );

    void addStringOption( const std::string name, std::string& option_field, std::string default_value );

    void addUnsignedLongOption( const std::string name, unsigned long& option_field, unsigned long default_value );

    void addUnsignedShortOption( const std::string name, unsigned short& option_field, unsigned short default_value );

    // enum types work differently than all of the others because there are a small number of valid
    // string entries for the type. One must also provide a list of all the valid strings of that type.
    template <class Tenum>
    void addEnumOption( const std::string name, Tenum& option_field, const std::map<std::string, Tenum>& enum_map, Tenum default_value );

    // List Options
    void addStringListOption( const std::string name, unsigned short& num_marker, std::string*& option_field );

  public:
    /*!
     * \brief Constructor of the class which reads the input file.
     */
    CConfig( char case_filename[MAX_STRING_SIZE] );

    /*!
     * \brief Destructor of the class.
     */
    ~CConfig( void );

    // ---- Getters for option values ----

    /*!
     * \brief Get Value of this option.
     *        Please keep alphabetical order within each subcategory
     */
    // File structure
    std::string inline GetMeshFile() const { return _inputDir + _meshFile; }
    std::string inline GetOutputDir() const {
        if( _outputDir.at( _outputDir.size() - 1 ) != '/' )
            return _outputDir + "/";
        else
            return _outputDir;
    }
    std::string inline GetOutputFile() const { return _outputFile; }
    std::string inline GetLogDir() const { return _logDir; }

    // Quadrature Structure
    QUAD_NAME inline GetQuadName() const { return _quadName; }
    unsigned short inline GetQuadOrder() const { return _quadOrder; }

    // Solver Structure
    double inline GetCFL() const { return _CFL; }
    double inline GetTEnd() const { return _tEnd; }

    // Boundary Conditions
    BOUNDARY_TYPE GetBoundaryType( std::string nameMarker ) const; /*! \brief Get Boundary Type of given marker */

    // Output Structure
};

#endif    // CONFIG_H
