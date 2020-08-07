#include "physics.h"

Physics::Physics( std::string fileName_H, std::string fileName_O, std::string fileName_stppower ) {
    _xsScatteringH2O.resize( 2 );
    _xsTotalH2O.resize( 2 );
    _xsTransportH2O.resize( 2 );
    LoadDatabase( fileName_H, fileName_O, fileName_stppower );
}

Physics::~Physics() {}

void Physics::LoadDatabase( std::string fileName_H, std::string fileName_O, std::string fileName_stppower ) {
    VectorVector transport_XS_H;
    VectorVector transport_XS_O;

    VectorVector scattering_XS_H;
    VectorVector scattering_XS_O;

    VectorVector total_XS_H;
    VectorVector total_XS_O;

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
            total_XS_H = data;
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
            total_XS_O = data;
        }
        // Angular distribution large angle elastic scattering XS
        else if( header[1][0] == 8 && header[1][1] == 22 ) {
            scattering_XS_O = data;
        }
    }

    _xsScatteringH2O[H] = scattering_XS_H;
    _xsScatteringH2O[O] = scattering_XS_O;

    _xsTransportH2O[H] = transport_XS_H;
    _xsTransportH2O[O] = transport_XS_O;

    _xsTotalH2O[H] = total_XS_H;
    _xsTotalH2O[O] = total_XS_O;

    _stpowH2O = ReadStoppingPowers( fileName_stppower );
}

VectorVector Physics::GetScatteringXS( Vector energies, Vector angle ) {
    std::vector<std::vector<double>> tmpH, tmpO;          // vectorvector which stores data at fixed energies
    std::vector<std::vector<double>> xsHGrid, xsOGrid;    // matrix which stores tensorized data for given angular grid, original energy grid
    VectorVector xsH2OGridGrid;                           // matrix which stores tensorized data for given angular and energy grid
    std::vector<double> tmpAngleGridH( angle.size() );
    std::vector<double> tmpAngleGridO( angle.size() );
    Vector tmpEnergyGridH( energies.size() );
    Vector tmpEnergyGridO( energies.size() );
    std::vector<double> tmp1H, tmp1O;
    std::vector<double> energiesOrig;    // original energy grid
    double energyCurrent = _xsScatteringH2O[H][0][0];

    // build grid at original energies and given angular grid
    for( unsigned i = 0; i < _xsScatteringH2O[H].size(); ++i ) {
        // split vector into subvectors with identical energy
        if( abs( _xsScatteringH2O[H][i][0] - energyCurrent ) < 1e-12 ) {
            if( tmpH.empty() ) tmpH.resize( 2 );
            tmpH[0].push_back( _xsScatteringH2O[H][i][1] );
            tmpH[1].push_back( _xsScatteringH2O[H][i][2] );
        }
        else {
            // new energy section starts in _xsH2O. Now we interpolate the values in tmp at given angular grid
            // interpolate vector at different angles for fixed current energies
            Interpolation xsH( tmpH[0], tmpH[1], Interpolation::linear );
            for( unsigned k = 0; k < angle.size(); k++ ) {
                tmpAngleGridH[k] = xsH( angle[k] );
            }
            xsHGrid.push_back( tmpAngleGridH );

            // reset current energy
            energiesOrig.push_back( energyCurrent );
            energyCurrent = _xsScatteringH2O[H][i][0];
            tmpH.clear();
        }
    }

    energyCurrent = _xsScatteringH2O[O][0][0];
    for( unsigned i = 0; i < _xsScatteringH2O[O].size(); ++i ) {
        // split vector into subvectors with identical energy
        if( abs( _xsScatteringH2O[O][i][0] - energyCurrent ) < 1e-12 ) {
            if( tmpO.empty() ) tmpO.resize( 2 );
            tmpO[0].push_back( _xsScatteringH2O[O][i][1] );
            tmpO[1].push_back( _xsScatteringH2O[O][i][2] );
        }
        else {
            // new energy section starts in _xsH2O. Now we interpolate the values in tmp at given angular grid
            // interpolate vector at different angles for fixed current energies
            Interpolation xsO( tmpO[0], tmpO[1], Interpolation::linear );
            for( unsigned k = 0; k < angle.size(); k++ ) {
                tmpAngleGridO[k] = xsO( angle[k] );
            }
            xsOGrid.push_back( tmpAngleGridO );

            // reset current energy
            energyCurrent = _xsScatteringH2O[O][i][0];
            tmpO.clear();
        }
    }

    // perform interpolation at fixed original energy for new energy grid

    for( unsigned j = 0; j < angle.size(); ++j ) {    // loop over all angles
        tmp1H.resize( 0 );
        tmp1O.resize( 0 );
        // store all values at given angle j for all original energies
        for( unsigned i = 0; i < energiesOrig.size(); ++i ) {
            tmp1H.push_back( xsHGrid[i][j] );
            tmp1O.push_back( xsOGrid[i][j] );
        }
        Interpolation xsH( energiesOrig, tmp1H, Interpolation::linear );
        Interpolation xsO( energiesOrig, tmp1O, Interpolation::linear );
        for( unsigned k = 0; k < energies.size(); k++ ) {
            // Linear interpolation
            tmpEnergyGridH[k] = xsH( energies[k] );
            tmpEnergyGridO[k] = xsO( energies[k] );
        }
        xsH2OGridGrid.push_back( _H20MassFractions[H] * tmpEnergyGridH + _H20MassFractions[O] * tmpEnergyGridO );
    }

    VectorVector xs( energies.size(), Vector( angle.size() ) );
    for( unsigned i = 0; i < energies.size(); ++i ) {
        for( unsigned j = 0; j < angle.size(); ++j ) {
            xs[i][j] = xsH2OGridGrid[j][i];
        }
    }

    return xs;
}

VectorVector Physics::GetTotalXS( Vector energies, Vector density ) {
    VectorVector total_XS( energies.size() );

    Interpolation xsH( _xsTotalH2O[H][0], _xsTotalH2O[H][1], Interpolation::linear );
    Interpolation xsO( _xsTotalH2O[O][0], _xsTotalH2O[O][1], Interpolation::linear );

    for( unsigned i = 0; i < energies.size(); i++ ) {
        total_XS[i] = ( _H20MassFractions[H] * xsH( energies[i] ) + _H20MassFractions[O] * xsO( energies[i] ) ) * density;
    }

    return total_XS;
}

Vector Physics::GetTotalXSE( Vector energies ) {
    Vector total_XS( energies.size() );

    Interpolation xsH( _xsTotalH2O[H][0], _xsTotalH2O[H][1], Interpolation::linear );
    Interpolation xsO( _xsTotalH2O[O][0], _xsTotalH2O[O][1], Interpolation::linear );

    for( unsigned i = 0; i < energies.size(); i++ ) {
        total_XS[i] = ( _H20MassFractions[H] * xsH( energies[i] ) + _H20MassFractions[O] * xsO( energies[i] ) );
    }

    return total_XS;
}

Vector Physics::GetStoppingPower( Vector energies ) {
    Vector stopping_power( energies.size() );

    Interpolation pwr( _stpowH2O[0], _stpowH2O[1], Interpolation::linear );

    for( unsigned i = 0; i < energies.size(); i++ ) {
        stopping_power[i] = pwr( energies[i] );
    }

    return stopping_power;
}

VectorVector Physics::GetTransportXS( Vector energies, Vector density ) {
    VectorVector transport_XS( energies.size() );

    Interpolation xsH( _xsTransportH2O[H][0], _xsTransportH2O[H][1], Interpolation::linear );
    Interpolation xsO( _xsTransportH2O[O][0], _xsTransportH2O[O][1], Interpolation::linear );

    for( unsigned i = 0; i < energies.size(); i++ ) {
        for( unsigned i = 0; i < energies.size(); i++ ) {
            transport_XS[i] = ( _H20MassFractions[H] * xsH( energies[i] ) + _H20MassFractions[O] * xsO( energies[i] ) ) * density;
        }
    }

    return transport_XS;
}

Vector Physics::GetTransportXSE( Vector energies ) {
    Vector transport_XS( energies.size() );

    Interpolation xsH( _xsTransportH2O[H][0], _xsTransportH2O[H][1], Interpolation::linear );
    Interpolation xsO( _xsTransportH2O[O][0], _xsTransportH2O[O][1], Interpolation::linear );

    for( unsigned i = 0; i < energies.size(); i++ ) {
        for( unsigned i = 0; i < energies.size(); i++ ) {
            transport_XS[i] = ( _H20MassFractions[H] * xsH( energies[i] ) + _H20MassFractions[O] * xsO( energies[i] ) );
        }
    }

    return transport_XS;
}

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
        else {
            // save dataLine in vector;
            if( header ) {
                headerPack[curr_HeaderLine] = dataLine;
                curr_HeaderLine++;
            }
            else {
                dataPackList.push_back( dataLine );
            }
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
