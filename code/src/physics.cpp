#include "physics.h"

#include "spline.h"
#include "toolboxes/errormessages.h"

#include <fstream>
#include <iostream>
#include <map>
#include <string>

Physics::Physics() {
    // LoadDatabase
}

void Physics::LoadDatabase( std::string fileName_H, std::string fileName_O, std::string fileName_stppower ) {
    VectorVector transport_XS_H;
    VectorVector transport_XS_O;

    VectorVector scattering_XS_H;
    VectorVector scattering_XS_O;

    VectorVector total_scat_XS_H;
    VectorVector total_scat_XS_O;

    std::vector<VectorVector> headers_H;
    std::vector<VectorVector> data_H;
    std::vector<VectorVector> headers_O;
    std::vector<VectorVector> data_O;

    // Read data for H
    std::tie( headers_H, data_H ) = ReadENDL( fileName_H );

    for( unsigned i = 0; i < headers_H.size(); ++i ) {
        auto header = headers_H[i];
        auto data   = data_H[i];
        if( header[1][0] == 7 && header[1][1] == 0 ) {
            transport_XS_H = data;
        }
        // Integrated elastic scattering XS
        else if( header[1][0] == 10 && header[1][1] == 0 ) {
            total_scat_XS_H = data;
        }
        // Angular distribution large angle elastic scattering XS
        else if( header[1][0] == 8 && header[1][1] == 22 ) {
            scattering_XS_H = data;
        }
    }

    // Read data for O
    std::tie( headers_O, data_O ) = ReadENDL( fileName_O );

    // Find required quantities and transfer to matrix
    for( unsigned i = 0; i < headers_O.size(); ++i ) {
        auto header = headers_O[i];
        auto data   = data_O[i];
        if( header[1][0] == 7 && header[1][1] == 0 ) {
            transport_XS_O = data;
        }
        // Integrated elastic scattering XS
        else if( header[1][0] == 10 && header[1][1] == 0 ) {
            total_scat_XS_O = data;
        }
        // Angular distribution large angle elastic scattering XS
        else if( header[1][0] == 8 && header[1][1] == 22 ) {
            scattering_XS_O = data;
        }
    }

    if( transport_XS_H.size() != transport_XS_O.size() )
        ErrorMessages::Error( "transport_XS of hydrogen and oxygen have different sizes", CURRENT_FUNCTION );
    if( scattering_XS_H.size() != scattering_XS_O.size() )
        ErrorMessages::Error( "scattering_XS of hydrogen and oxygen have different sizes", CURRENT_FUNCTION );
    if( total_scat_XS_H.size() != total_scat_XS_O.size() )
        ErrorMessages::Error( "total_XS of hydrogen and oxygen have different sizes", CURRENT_FUNCTION );

    // Combine values for H and O according to mixture ratio in water
    for( unsigned i = 0; i < transport_XS_H.size(); ++i )
        _xsTransportH2O.push_back( 0.11189400 * transport_XS_H[i] + 0.88810600 * transport_XS_O[i] );
    for( unsigned i = 0; i < scattering_XS_H.size(); ++i ) _xsH2O.push_back( 0.11189400 * scattering_XS_H[i] + 0.88810600 * scattering_XS_O[i] );
    for( unsigned i = 0; i < total_scat_XS_H.size(); ++i ) _xsTotalH2O.push_back( 0.11189400 * total_scat_XS_H[i] + 0.88810600 * total_scat_XS_O[i] );

    _stpowH2O = ReadStoppingPowers( fileName_stppower );
}

VectorVector Physics::GetScatteringXS( Vector energies, Vector density, Vector angle ) {

    VectorVector scattering_XS;
    return scattering_XS;
}

VectorVector Physics::GetTotalXS( Vector energies, Vector density ) {
    VectorVector total_XS;
    double xsH2O_i;

    Spline interp;
    interp.set_points( _xsTotalH2O[0], _xsTotalH2O[1], false );    // false == linear interpolation

    for( unsigned i = 0; i < energies.size(); i++ ) {

        if( energies[i] < _xsTotalH2O[1][0] ) {
            xsH2O_i = _xsTotalH2O[1][0];
        }
        else if( energies[i] > _xsTotalH2O[1][energies.size() - 1] ) {
            xsH2O_i = _xsTotalH2O[1][energies.size() - 1];
        }
        else {
            // Linear interpolation
            xsH2O_i = interp( energies[i] );
        }

        total_XS[i] = xsH2O_i * density;
    }

    return total_XS;
}

VectorVector Physics::GetStoppingPower( Vector energies, Vector density ) {
    VectorVector stopping_power;
    double stpw_H2O_i;

    Spline interp;
    interp.set_points( _stpowH2O[0], _stpowH2O[1], false );    // false == linear interpolation

    for( unsigned i = 0; i < energies.size(); i++ ) {

        if( energies[i] < _stpowH2O[1][0] ) {
            stpw_H2O_i = _stpowH2O[1][0];
        }
        else if( energies[i] > _stpowH2O[1][energies.size() - 1] ) {
            stpw_H2O_i = _stpowH2O[1][energies.size() - 1];
        }
        else {
            // Linear interpolation
            stpw_H2O_i = interp( energies[i] );
        }

        stopping_power[i] = stpw_H2O_i * density;
    }
    return stopping_power;
}

VectorVector Physics::GetTransportXS( Vector energies, Vector density ) {
    VectorVector transport_XS;
    double xsH2O_i;

    Spline interp;
    interp.set_points( _xsTransportH2O[0], _xsTransportH2O[1], false );    // false == linear interpolation

    for( unsigned i = 0; i < energies.size(); i++ ) {

        if( energies[i] < _xsTransportH2O[1][0] ) {
            xsH2O_i = _xsTransportH2O[1][0];
        }
        else if( energies[i] > _xsTransportH2O[1][energies.size() - 1] ) {
            xsH2O_i = _xsTransportH2O[1][energies.size() - 1];
        }
        else {
            // Linear interpolation
            xsH2O_i = interp( energies[i] );
        }

        transport_XS[i] = xsH2O_i * density;
    }
    return transport_XS;
}

Physics* Physics::Create() { return new Physics(); }

std::tuple<std::vector<VectorVector>, std::vector<VectorVector>> Physics::ReadENDL( std::string filename ) {
    std::vector<VectorVector> _headers;
    std::vector<VectorVector> _data;

    // open txt file
    std::string text_line, option_name;
    std::ifstream case_file;
    std::vector<std::string> option_value;

    /*--- Read the configuration file ---*/

    case_file.open( filename, std::ios::in );

    if( case_file.fail() ) {
        ErrorMessages::Error( "The ENDL file is missing!!", CURRENT_FUNCTION );
    }

    std::map<std::string, bool> included_options;

    /*--- Parse the physics file and save the values ---*/

    // list of entries (dataLine) of a datapackage, later becomes one entry in the list _data;
    std::list<Vector> dataPackList;
    VectorVector headerPack( 2 );

    // indices
    unsigned pack_idx        = 0;
    unsigned curr_HeaderLine = 0;

    while( getline( case_file, text_line ) ) {

        std::stringstream ss( text_line );    // give line to stringstream
        std::list<double> linelist;
        for( double double_in; ss >> double_in; ) {    // parse line
            linelist.push_back( double_in );
        }

        // Save content of line
        Vector dataLine( linelist.size() );
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

            VectorVector dataPack( dataPackList.size() );
            for( unsigned idx = 0; idx < dataPack.size(); idx++ ) {
                dataPack[idx] = dataPackList.front();
                dataPackList.pop_front();
            }
            _data.push_back( dataPack );

            pack_idx++;
            curr_HeaderLine = 0;
        }
    }
    case_file.close();
    return { _headers, _data };
}

VectorVector Physics::ReadStoppingPowers( std::string fileName ) {
    VectorVector stp_powers;
    std::string text_line;

    std::vector<double> e, stp_power;

    std::ifstream file_stream;
    file_stream.open( fileName, std::ios::in );

    if( file_stream.fail() ) {
        ErrorMessages::Error( "The Stopping Power file is missing!!", CURRENT_FUNCTION );
    }

    while( getline( file_stream, text_line ) ) {

        std::stringstream ss( text_line );    // give line to stringstream
        std::list<double> linelist;
        for( double double_in; ss >> double_in; ) {    // parse line
            linelist.push_back( double_in );
        }

        e.push_back( linelist.front() );    // append elements from line to column vectors
        linelist.pop_front();
        stp_power.push_back( linelist.front() );
        linelist.pop_front();
    }
    file_stream.close();

    stp_powers.push_back( Vector( e.size(), e.data() ) );
    stp_powers.push_back( Vector( stp_power.size(), stp_power.data() ) );

    return stp_powers;
}
