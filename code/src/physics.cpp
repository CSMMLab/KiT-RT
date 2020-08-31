#include "physics.h"

Physics::Physics( std::string fileName_H, std::string fileName_O, std::string fileName_stppower ) {
    _xsScatteringH2O.resize( 2 );
    _xsTotalH2O.resize( 2 );
    LoadDatabase( fileName_H, fileName_O, fileName_stppower );
}

Physics::~Physics() {}

void Physics::LoadDatabase( std::string fileName_H, std::string fileName_O, std::string fileName_stppower ) {
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

        // Integrated elastic scattering XS
        if( header[1][0] == 10 && header[1][1] == 0 ) {
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
        // Integrated elastic scattering XS
        if( header[1][0] == 10 && header[1][1] == 0 ) {
            total_XS_O = data;
        }
        // Angular distribution large angle elastic scattering XS
        else if( header[1][0] == 8 && header[1][1] == 22 ) {
            scattering_XS_O = data;
        }
    }

    for( auto& sXS : scattering_XS_H ) sXS[1] = 1.0 - sXS[1];    // data is given as x=1-cos(theta) in [0,2] -> shift to [-1,1] interval
    for( auto& sXS : scattering_XS_O ) sXS[1] = 1.0 - sXS[1];

    _xsScatteringH2O[H] = scattering_XS_H;
    _xsScatteringH2O[O] = scattering_XS_O;

    _xsTotalH2O[H] = total_XS_H;
    _xsTotalH2O[O] = total_XS_O;

    if( !( fileName_stppower.empty() ) ) _stpowH2O = ReadStoppingPowers( fileName_stppower );
}

std::vector<Matrix> Physics::GetScatteringXS( const Vector& energies, const Matrix& angle ) {
    std::vector<Matrix> out( energies.size(), Matrix( angle.rows(), angle.columns() ) );
    Vector angleVec( angle.columns() * angle.rows() );
    // store Matrix with mu values in vector format to call GetScatteringXS
    for( unsigned i = 0; i < angle.rows(); ++i ) {
        for( unsigned j = 0; j < angle.columns(); ++j ) {
            angleVec[i * angle.columns() + j] = angle( i, j );
        }
    }

    VectorVector outVec = GetScatteringXS( energies, angleVec );

    // rearrange output to matrix format
    for( unsigned n = 0; n < energies.size(); ++n ) {
        for( unsigned i = 0; i < angle.rows(); ++i ) {
            for( unsigned j = 0; j < angle.columns(); ++j ) {
                out[n]( i, j ) = outVec[n][i * angle.columns() + j];
            }
        }
    }
    return out;
}

VectorVector Physics::GetScatteringXS( Vector energies, Vector angle ) {
    auto scatter_XS_H = _xsScatteringH2O[H];
    auto scatter_XS_O = _xsScatteringH2O[O];

    Vector integratedScatteringXS = GetTotalXSE( energies );

    std::vector<double> dataEnergy;
    std::vector<Vector> intermediateGrid;
    unsigned idx = 0;
    for( unsigned i = 0; i < scatter_XS_H.size(); i += idx ) {
        idx = 1;
        std::vector<double> dataAngle{ scatter_XS_H[i][1] }, dataXS{ scatter_XS_H[i][2] };
        dataEnergy.push_back( scatter_XS_H[i][0] );
        while( scatter_XS_H[i + idx][0] == scatter_XS_H[i + idx - 1][0] ) {
            dataAngle.push_back( scatter_XS_H[i + idx][1] );
            dataXS.push_back( scatter_XS_H[i + idx][2] );
            idx++;
            if( i + idx >= scatter_XS_H.size() ) break;
        }
        std::reverse( dataAngle.begin(), dataAngle.end() );
        std::reverse( dataXS.begin(), dataXS.end() );

        Interpolation interp( dataAngle, dataXS, Interpolation::linear );
        intermediateGrid.push_back( interp( angle ) );
        double integral_sXS = 0.0;
        auto sXS            = intermediateGrid.back();
        for( unsigned a = 1; a < angle.size(); ++a ) {
            integral_sXS += 0.5 * ( angle[a] - angle[a - 1] ) * ( sXS[a] + sXS[a - 1] );
        }
        // intermediateGrid.back() /= integral_sXS;    // re-normalize to yield a valid distribution in a statistical sense; behaves poorly due to
        // great
        // differences in order of magnitude
    }
    VectorVector xs( energies.size(), Vector( angle.size() ) );
    for( unsigned j = 0; j < angle.size(); ++j ) {
        std::vector<double> xsPerE( dataEnergy.size() );
        for( unsigned i = 0; i < dataEnergy.size(); ++i ) {
            xsPerE[i] = intermediateGrid[i][j];
        }
        Interpolation interp( dataEnergy, xsPerE, Interpolation::linear );
        Vector eGrid = interp( energies );
        for( unsigned i = 0; i < energies.size(); ++i ) {
            // xs[i][j] = eGrid[i] * integratedScatteringXS[i];    // multiply with the integrated scattering cross section to force the
            // integration over the angle to yield integratedScatteringXS
        }
    }

    return xs;
}

VectorVector Physics::GetTotalXS( Vector energies, Vector density ) {
    VectorVector total_XS( energies.size() );

    Interpolation xsH( _xsTotalH2O[H][0], _xsTotalH2O[H][1], Interpolation::linear );
    Interpolation xsO( _xsTotalH2O[O][0], _xsTotalH2O[O][1], Interpolation::linear );

    for( unsigned i = 0; i < energies.size(); i++ ) {
        total_XS[i] = ( H2OMassFractions[H] * xsH( energies[i] ) * 1e-24 + H2OMassFractions[O] * xsO( energies[i] ) * 1e-24 ) * density;
    }

    return total_XS;
}

VectorVector Physics::ReorderScatteringXS( const VectorVector& data ) const {
    VectorVector ret( data[0].size(), Vector( data.size() ) );
    for( unsigned i = 0; i < data.size(); i++ ) {
        for( unsigned j = 0; j < data[0].size(); j++ ) {
            ret[j][i] = data[i][j];
        }
    }
    return ret;
}

VectorVector Physics::ReorderTotalXS( const VectorVector& data ) const {
    VectorVector ret( data[0].size(), Vector( data.size() ) );
    for( unsigned i = 0; i < data.size(); i++ ) {
        for( unsigned j = 0; j < data[0].size(); j++ ) {
            ret[j][i] = data[i][j];
        }
    }
    return ret;
}

Vector Physics::GetTotalXSE( Vector energies ) {
    Vector total_XS( energies.size() );

    auto total_XS_H = ReorderTotalXS( _xsTotalH2O[H] );
    auto total_XS_O = ReorderTotalXS( _xsTotalH2O[O] );

    Interpolation xsH( total_XS_H[0], total_XS_H[1], Interpolation::linear );
    Interpolation xsO( total_XS_O[0], total_XS_O[1], Interpolation::linear );

    for( unsigned i = 0; i < energies.size(); i++ ) {
        total_XS[i] = ( H2OMassFractions[H] * xsH( energies[i] ) + H2OMassFractions[O] * xsO( energies[i] ) ) * 1e-24;
    }

    return total_XS;
}

Vector Physics::GetStoppingPower( Vector energies ) {
    if( _stpowH2O.empty() ) {
        return ComputeStoppingPower( energies );
    }
    else {
        Vector stopping_power( energies.size() );
        Interpolation pwr( _stpowH2O[0], _stpowH2O[1], Interpolation::linear );
        for( unsigned i = 0; i < energies.size(); i++ ) {
            stopping_power[i] = pwr( energies[i] );
        }
        return stopping_power;
    }
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

Vector Physics::ComputeStoppingPower( const Vector& energies ) const {

    Vector stoppingPower( energies.size() );
    for( unsigned i = 0; i < energies.size(); ++i ) {
        long double e = energies[i] * 1e6 * ELECTRON_VOLT;    // MeV -> J
        long double J = 0, S = 0;
        for( int i = 0; i < 2; i++ ) {
            if( H2OAtomicNumbers[i] <= 6 )
                J = ELECTRON_VOLT * 11.5 * H2OAtomicNumbers[i];
            else
                J = ELECTRON_VOLT * ( 9.76 * H2OAtomicNumbers[i] + 58.8 * std::pow( H2OAtomicNumbers[i], -0.19 ) );
            S += H2OMassFractions[i] * H2OAtomicNumbers[i] / H2OAtomicWeights[i] * ( AVOGADRO_CONSTANT * 1e3 ) * std::log( SQRT05E / J * e );
        }
        S = C1 * S / e;

        stoppingPower[i] = static_cast<double>( S / ( ELECTRON_VOLT * H2OMassDensity * 1e5 ) );
    }
    /* OpenMC formula
    Vector stoppingPower( energies.size() );
    for( unsigned i = 0; i < energies.size(); ++i ) {
        long double S_rad = 0, S_col = 0;
        for( int i = 0; i < 2; i++ ) {
            S_rad += H20MassFractions[i] * ( 1.0 );
            long double beta  = 0;
            long double T     = 0;
            long double I     = 0;
            long double tau   = T / ELECTRON_MASS;
            long double F     = ( 1 - beta * beta ) * ( 1 + tau * tau / 8 - ( 2 * tau + 1 * M_LN2 ) );
            long double delta = 0;
            S_col += H20MassFractions[i] * ( 2.0 * PI * ELECTRON_RADIUS * ELECTRON_ENERGY / beta * AVOGADRO_CONSTANT * H20AtomicNumbers[i] /
                                             H2OAtomicWeights[i] * ( std::log( ( T * T ) / ( I * I ) ) + std::log( 1.0 + tau / 2.0 ) + F - delta )
    );
        }
        stoppingPower[i] = static_cast<double>( S_rad + S_col );
    }
    */
    return stoppingPower;
}
