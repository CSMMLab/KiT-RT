/*!
 * \file CConfig.cpp
 * \brief Classes for different Options in rtsn
 * \author S. Schotthoefer
 *
 * Disclaimer: This class structure was copied and modifed with open source permission from SU2 v7.0.3 https://su2code.github.io/
 */

#include "settings/CConfig.h"
#include "settings/GlobalConstants.h"
#include "toolboxes/CRTSNError.h"
#include "toolboxes/NIoToolBox.h"
#include "toolboxes/TextProcessingToolbox.h"
#include <cassert>
#include <fstream>

using namespace std;

CConfig::CConfig( char case_filename[MAX_STRING_SIZE] ) {

    /*--- Set the case name to the base config file name without extension ---*/

    _fileName = NIoToolBox::split( string( case_filename ), '.' )[0];
    _inputDir = std::filesystem::path( case_filename ).parent_path().string();
    if( _inputDir.substr( _inputDir.size() - 1 ) != "/" ) _inputDir.append( "/" );

    _base_config = true;

    /*--- Store MPI rank and size ---*/
    // TODO with MPI implementation

    /*--- Initialize pointers to Null---*/

    SetPointersNull();

    /*--- Reading config options  ---*/

    SetConfigOptions();

    /*--- Parsing the config file  ---*/

    SetConfigParsing( case_filename );

    /*--- Set the default values for all of the options that weren't set ---*/

    SetDefault();

    /*--- Get the Mesh Valueso--- */

    // val_nDim = GetnDim(Mesh_FileName, Mesh_FileFormat);
    // TODO

    /*--- Configuration file postprocessing ---*/

    SetPostprocessing();

    /*--- Configuration file boundaries/markers setting ---*/

    // SetBoundary(); IDK how boundaries are implemented, but i think this should be treated here

    /*--- Configuration file output ---*/

    // if ((rank == MASTER_NODE))
    SetOutput();
}

CConfig::~CConfig( void ) {
    // Delete all introduced arrays!

    // Boundary Conditions
    if( _MarkerDirichlet != NULL ) delete[] _MarkerDirichlet;
}

// ---- Add Options ----

// Simple Options
void CConfig::addBoolOption( const string name, bool& option_field, bool default_value ) {
    // Check if the key is already in the map. If this fails, it is coder error
    // and not user error, so throw.
    assert( option_map.find( name ) == option_map.end() );

    // Add this option to the list of all the options
    all_options.insert( pair<string, bool>( name, true ) );

    // Create the parser for a bool option with a reference to the option_field and the desired
    // default value. This will take the string in the config file, convert it to a bool, and
    // place that bool in the memory location specified by the reference.
    COptionBase* val = new COptionBool( name, option_field, default_value );

    // Create an association between the option name ("CFL") and the parser generated above.
    // During configuration, the parsing script will get the option name, and use this map
    // to find how to parse that option.
    option_map.insert( pair<string, COptionBase*>( name, val ) );
}

void CConfig::addDoubleOption( const string name, double& option_field, double default_value ) {
    assert( option_map.find( name ) == option_map.end() );
    all_options.insert( pair<string, bool>( name, true ) );
    COptionBase* val = new COptionDouble( name, option_field, default_value );
    option_map.insert( pair<string, COptionBase*>( name, val ) );
}

void CConfig::addIntegerOption( const string name, int& option_field, int default_value ) {
    assert( option_map.find( name ) == option_map.end() );
    all_options.insert( pair<string, bool>( name, true ) );
    COptionBase* val = new COptionInt( name, option_field, default_value );
    option_map.insert( pair<string, COptionBase*>( name, val ) );
}

void CConfig::addLongOption( const string name, long& option_field, long default_value ) {
    assert( option_map.find( name ) == option_map.end() );
    all_options.insert( pair<string, bool>( name, true ) );
    COptionBase* val = new COptionLong( name, option_field, default_value );
    option_map.insert( pair<string, COptionBase*>( name, val ) );
}

void CConfig::addStringOption( const string name, string& option_field, string default_value ) {
    assert( option_map.find( name ) == option_map.end() );
    all_options.insert( pair<string, bool>( name, true ) );
    COptionBase* val = new COptionString( name, option_field, default_value );
    option_map.insert( pair<string, COptionBase*>( name, val ) );
}

void CConfig::addUnsignedLongOption( const string name, unsigned long& option_field, unsigned long default_value ) {
    assert( option_map.find( name ) == option_map.end() );
    all_options.insert( pair<string, bool>( name, true ) );
    COptionBase* val = new COptionULong( name, option_field, default_value );
    option_map.insert( pair<string, COptionBase*>( name, val ) );
}

void CConfig::addUnsignedShortOption( const string name, unsigned short& option_field, unsigned short default_value ) {
    assert( option_map.find( name ) == option_map.end() );
    all_options.insert( pair<string, bool>( name, true ) );
    COptionBase* val = new COptionUShort( name, option_field, default_value );
    option_map.insert( pair<string, COptionBase*>( name, val ) );
}

// enum types work differently than all of the others because there are a small number of valid
// string entries for the type. One must also provide a list of all the valid strings of that type.
template <class Tenum>
void CConfig::addEnumOption( const string name, Tenum& option_field, const map<string, Tenum>& enum_map, Tenum default_value ) {
    assert( option_map.find( name ) == option_map.end() );
    all_options.insert( pair<string, bool>( name, true ) );
    COptionBase* val = new COptionEnum<Tenum>( name, enum_map, option_field, default_value );
    option_map.insert( pair<string, COptionBase*>( name, val ) );
    return;
}

// List Options
void CConfig::addStringListOption( const string name, unsigned short& num_marker, std::string*& option_field ) {
    assert( option_map.find( name ) == option_map.end() );
    all_options.insert( pair<string, bool>( name, true ) );
    COptionBase* val = new COptionStringList( name, num_marker, option_field );
    option_map.insert( pair<string, COptionBase*>( name, val ) );
}

// ---- Getter Functions ----

BOUNDARY_TYPE CConfig::GetBoundaryType( std::string name ) const {
    for( unsigned i = 0; i < _boundaries.size(); ++i ) {
        if( name == _boundaries[i].first ) return _boundaries[i].second;
    }
    return BOUNDARY_TYPE::INVALID;
}

// ---- Setter Functions ----

void CConfig::SetDefault() {
    /*--- Set the default values for all of the options that weren't set ---*/

    for( map<string, bool>::iterator iter = all_options.begin(); iter != all_options.end(); ++iter ) {
        if( option_map[iter->first]->GetValue().size() == 0 ) option_map[iter->first]->SetDefault();
    }
}

void CConfig::SetConfigOptions() {

    /* BEGIN_CONFIG_OPTIONS */

    /*!\par CONFIG_CATEGORY: Problem Definition \ingroup Config */
    /*--- Options related to problem definition and partitioning ---*/

    // File Structure related options
    /*!\brief OUTPUT_DIR \n DESCRIPTION: Relative Directory of output files \n DEFAULT "/out" \ingroup Config.*/
    addStringOption( "OUTPUT_DIR", _outputDir, string( "/out" ) );
    /*!\brief OUTPUT_FILE \n DESCRIPTION: Name of output file \n DEFAULT "output" \ingroup Config.*/
    addStringOption( "OUTPUT_FILE", _outputFile, string( "output" ) );
    /*!\brief LOG_DIR \n DESCRIPTION: Relative Directory of log files \n DEFAULT "/out" \ingroup Config.*/
    addStringOption( "LOG_DIR", _logDir, string( "/out" ) );
    /*!\brief MESH_FILE \n DESCRIPTION: Name of mesh file \n DEFAULT "" \ingroup Config.*/
    addStringOption( "MESH_FILE", _meshFile, string( "mesh.su2" ) );

    // Quadrature relatated options
    /*!\brief QUAD_TYPE \n DESCRIPTION: Type of Quadrature rule \n Options: see \link QUAD_NAME \endlink \n DEFAULT: QUAD_MonteCarlo \ingroup Config*/
    addEnumOption( "QUAD_TYPE", _quadName, Quadrature_Map, QUAD_MonteCarlo );
    /*!\brief QUAD_ORDER \n DESCRIPTION: Order of Quadrature rule \n DEFAULT 2 \ingroup Config.*/
    addUnsignedShortOption( "QUAD_ORDER", _quadOrder, 2 );

    // Solver related options
    /*!\brief CFL \n DESCRIPTION: CFL number \n DEFAULT 1.0 \ingroup Config.*/
    addDoubleOption( "CFL_NUMBER", _CFL, 1.0 );
    /*!\brief TIME_FINAL \n DESCRIPTION: Final time for simulation \n DEFAULT 1.0 \ingroup Config.*/
    addDoubleOption( "TIME_FINAL", _tEnd, 1.0 );

    // Mesh related options
    // Boundary Markers
    /*!\brief BC_DIRICHLET\n DESCRIPTION: Dirichlet wall boundary marker(s) \ingroup Config*/
    addStringListOption( "BC_DIRICHLET", _nMarkerDirichlet, _MarkerDirichlet );
}

void CConfig::SetConfigParsing( char case_filename[MAX_STRING_SIZE] ) {
    string text_line, option_name;
    ifstream case_file;
    vector<string> option_value;

    /*--- Read the configuration file ---*/

    case_file.open( case_filename, ios::in );

    if( case_file.fail() ) {
        CRTSNError::Error( "The configuration file (.cfg) is missing!!", CURRENT_FUNCTION );
    }

    string errorString;

    int err_count     = 0;     // How many errors have we found in the config file
    int max_err_count = 30;    // Maximum number of errors to print before stopping

    map<string, bool> included_options;

    /*--- Parse the configuration file and set the options ---*/

    while( getline( case_file, text_line ) ) {

        if( err_count >= max_err_count ) {
            errorString.append( "too many errors. Stopping parse" );

            cout << errorString << endl;
            throw( 1 );
        }

        if( TokenizeString( text_line, option_name, option_value ) ) {

            /*--- See if it's a python option ---*/

            if( option_map.find( option_name ) == option_map.end() ) {
                string newString;
                newString.append( option_name );
                newString.append( ": invalid option name" );
                newString.append( ". Check current RTSN options in config_template.cfg." );
                newString.append( "\n" );
                errorString.append( newString );
                err_count++;
                continue;
            }

            /*--- Option exists, check if the option has already been in the config file ---*/

            if( included_options.find( option_name ) != included_options.end() ) {
                string newString;
                newString.append( option_name );
                newString.append( ": option appears twice" );
                newString.append( "\n" );
                errorString.append( newString );
                err_count++;
                continue;
            }

            /*--- New found option. Add it to the map, and delete from all options ---*/

            included_options.insert( pair<string, bool>( option_name, true ) );
            all_options.erase( option_name );

            /*--- Set the value and check error ---*/

            string out = option_map[option_name]->SetValue( option_value );
            if( out.compare( "" ) != 0 ) {
                errorString.append( out );
                errorString.append( "\n" );
                err_count++;
            }
        }
    }

    /*--- See if there were any errors parsing the config file ---*/

    if( errorString.size() != 0 ) {
        CRTSNError::Error( errorString, CURRENT_FUNCTION );
    }

    case_file.close();
}

void CConfig::SetPointersNull( void ) {
    // All pointer valued options should be set to NULL here
    _MarkerDirichlet = NULL;
}

void CConfig::SetPostprocessing() {
    // Check if there are contradictive or false options set.

    // Regroup Boundary Conditions to  std::vector<std::pair<std::string, BOUNDARY_TYPE>> _boundaries;
    {
        for( int i = 0; i < _nMarkerDirichlet; i++ ) {
            _boundaries.push_back( std::pair<std::string, BOUNDARY_TYPE>( _MarkerDirichlet[i], DIRICHLET ) );
        }
    }
}

void CConfig::SetOutput() {
    // Set Output for settings, i.e. feedback on what has been chosen
}

bool CConfig::TokenizeString( string& str, string& option_name, vector<string>& option_value ) {
    const string delimiters( " (){}:,\t\n\v\f\r" );
    // check for comments or empty string
    string::size_type pos, last_pos;
    pos = str.find_first_of( "%" );
    if( ( str.length() == 0 ) || ( pos == 0 ) ) {
        // str is empty or a comment line, so no option here
        return false;
    }
    if( pos != string::npos ) {
        // remove comment at end if necessary
        str.erase( pos );
    }

    // look for line composed on only delimiters (usually whitespace)
    pos = str.find_first_not_of( delimiters );
    if( pos == string::npos ) {
        return false;
    }

    // find the equals sign and split string
    string name_part, value_part;
    pos = str.find( "=" );
    if( pos == string::npos ) {
        cerr << "Error in TokenizeString(): "
             << "line in the configuration file with no \"=\" sign." << endl;
        cout << "Look for: " << str << endl;
        cout << "str.length() = " << str.length() << endl;
        throw( -1 );
    }
    name_part  = str.substr( 0, pos );
    value_part = str.substr( pos + 1, string::npos );
    // cout << "name_part  = |" << name_part  << "|" << endl;
    // cout << "value_part = |" << value_part << "|" << endl;

    // the first_part should consist of one string with no interior delimiters
    last_pos = name_part.find_first_not_of( delimiters, 0 );
    pos      = name_part.find_first_of( delimiters, last_pos );
    if( ( name_part.length() == 0 ) || ( last_pos == string::npos ) ) {
        cerr << "Error in CConfig::TokenizeString(): "
             << "line in the configuration file with no name before the \"=\" sign." << endl;
        throw( -1 );
    }
    if( pos == string::npos ) pos = name_part.length();
    option_name = name_part.substr( last_pos, pos - last_pos );
    last_pos    = name_part.find_first_not_of( delimiters, pos );
    if( last_pos != string::npos ) {
        cerr << "Error in TokenizeString(): "
             << "two or more options before an \"=\" sign in the configuration file." << endl;
        throw( -1 );
    }
    StringToUpperCase( option_name );

    // cout << "option_name = |" << option_name << "|" << endl;
    // cout << "pos = " << pos << ": last_pos = " << last_pos << endl;

    // now fill the option value vector
    option_value.clear();
    last_pos = value_part.find_first_not_of( delimiters, 0 );
    pos      = value_part.find_first_of( delimiters, last_pos );
    while( string::npos != pos || string::npos != last_pos ) {
        // add token to the vector<string>
        option_value.push_back( value_part.substr( last_pos, pos - last_pos ) );
        // skip delimiters
        last_pos = value_part.find_first_not_of( delimiters, pos );
        // find next "non-delimiter"
        pos = value_part.find_first_of( delimiters, last_pos );
    }
    if( option_value.size() == 0 ) {
        cerr << "Error in TokenizeString(): "
             << "option " << option_name << " in configuration file with no value assigned." << endl;
        throw( -1 );
    }

#if 0
  cout << "option value(s) = ";
  for (unsigned int i = 0; i < option_value.size(); i++)
    cout << option_value[i] << " ";
  cout << endl;
#endif

    // look for ';' DV delimiters attached to values
    vector<string>::iterator it;
    it = option_value.begin();
    while( it != option_value.end() ) {
        if( it->compare( ";" ) == 0 ) {
            it++;
            continue;
        }

        pos = it->find( ';' );
        if( pos != string::npos ) {
            string before_semi = it->substr( 0, pos );
            string after_semi  = it->substr( pos + 1, string::npos );
            if( before_semi.empty() ) {
                *it = ";";
                it++;
                option_value.insert( it, after_semi );
            }
            else {
                *it = before_semi;
                it++;
                vector<string> to_insert;
                to_insert.push_back( ";" );
                if( !after_semi.empty() ) to_insert.push_back( after_semi );
                option_value.insert( it, to_insert.begin(), to_insert.end() );
            }
            it = option_value.begin();    // go back to beginning; not efficient
            continue;
        }
        else {
            it++;
        }
    }
#if 0
  cout << "option value(s) = ";
  for (unsigned int i = 0; i < option_value.size(); i++)
    cout << option_value[i] << " ";
  cout << endl;
#endif
    // remove any consecutive ";"
    it                = option_value.begin();
    bool semi_at_prev = false;
    while( it != option_value.end() ) {
        if( semi_at_prev ) {
            if( it->compare( ";" ) == 0 ) {
                option_value.erase( it );
                it           = option_value.begin();
                semi_at_prev = false;
                continue;
            }
        }
        if( it->compare( ";" ) == 0 ) {
            semi_at_prev = true;
        }
        else {
            semi_at_prev = false;
        }
        it++;
    }

#if 0
  cout << "option value(s) = ";
  for (unsigned int i = 0; i < option_value.size(); i++)
    cout << option_value[i] << " ";
  cout << endl;
#endif
    return true;
}
