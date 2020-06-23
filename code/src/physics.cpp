#include "physics.h"
#include "toolboxes/errormessages.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

Physics::Physics( Config* settings ) {
    // @TODO set up physical model from settings?
    // @TODO What else belongs in physics class?
    // load H2O cross sections
    //  LoadXSH2O( "test1", "test2" );

    // Load physics from ENDL_H file
}

void Physics::LoadXSH2O( std::string fileName1, std::string fileName2 ) {
    // @TODO
}

VectorVector Physics::GetScatteringXS( std::vector<double> energies, std::vector<double> density, std::vector<double> Omegas ) {
    // @TODO
    VectorVector scattering_XS;
    return scattering_XS;
}

VectorVector Physics::GetTotalXS( std::vector<double> energies, std::vector<double> density ) {
    // @TODO
    VectorVector total_XS;
    return total_XS;
}

VectorVector Physics::GetStoppingPower( std::vector<double> energies, std::vector<double> sH2O ) {
    // @TODO
    VectorVector stopping_power;
    return stopping_power;
}

Physics* Physics::Create( Config* settings ) { return new Physics( settings ); }

void Physics::ReadENDL_H( string filename ) {
    // open txt file
    string text_line, option_name;
    ifstream case_file;
    vector<string> option_value;

    /*--- Read the configuration file ---*/

    case_file.open( filename, ios::in );

    if( case_file.fail() ) {
        ErrorMessages::Error( "The ENDL_H file is missing!!", CURRENT_FUNCTION );
    }

    map<string, bool> included_options;

    /*--- Parse the physics file and save the values ---*/

    // list of entries (dataLine) of a datapackage, later becomes one entry in the list _data;
    list<vector<double>> dataPackList;
    vector<vector<double>> headerPack( 2 );

    // indices
    unsigned pack_idx        = 0;
    unsigned curr_HeaderLine = 0;

    while( getline( case_file, text_line ) ) {

        std::stringstream ss( text_line );    // give line to stringstream
        list<double> linelist;
        for( double double_in; ss >> double_in; ) {    // parse line
            linelist.push_back( double_in );
        }

        // Save content of line
        vector<double> dataLine( linelist.size() );
        bool header         = false;
        unsigned lineLenght = linelist.size();
        // check, if line is a dataline
        if( lineLenght == 2 || lineLenght == 3 ) {
            dataLine[0] = linelist.front();
            linelist.pop_front();
            dataLine[1] = linelist.front();
            linelist.pop_front();
        }
        if( lineLenght == 3 ) {
            dataLine[2] = linelist.front();
            linelist.pop_front();
        }
        // check, if line is  a header
        if( lineLenght == 8 || lineLenght == 9 ) {
            for( unsigned data_idx = 0; data_idx < lineLenght; data_idx++ ) {
                dataLine[data_idx] = linelist.front();
                linelist.pop_front();
                header = true;
            }
        }

        // save dataLine in vector;
        if( header ) {
            headerPack[curr_HeaderLine] = dataLine;
            curr_HeaderLine++;
        }
        else {
            dataPackList.push_back( dataLine );
        }

        // check if line is an endline, then reset indices and save packs
        if( lineLenght == 1 ) {

            _headers.push_back( headerPack );

            vector<vector<double>> dataPack( dataPackList.size() );
            for( unsigned idx = 0; idx < dataPack.size(); idx++ ) {
                dataPack[idx] = dataPackList.front();
                dataPackList.pop_front();
            }
            _data.push_back( dataPack );

            pack_idx++;
            curr_HeaderLine = 0;
        }
    }
    cout << pack_idx;
}
