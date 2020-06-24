#include "physics.h"
using namespace std;
using blaze::CompressedVector;

#include "toolboxes/errormessages.h"
#include <fstream>
#include <iostream>
#include <string>

Physics::Physics (){;
//LoadDatabase
}

void Physics::LoadDatabase( std::string fileName_H, std::string fileName_O, std::string fileName_stppower){
    VectorVector transport_XS_H;
    VectorVector transport_XS_O;
	
	VectorVector scattering_XS_H;
    VectorVector scattering_XS_O;
	
	VectorVector total_scat_XS_H;
    VectorVector total_scat_XS_O;
	
	VectorVector stopp_pow_H;
    VectorVector stopp_pow_O;
	
	VectorVector header;
	VectorVector data;
    
	//Read data for H
	[headers_H,data_H] = ReadENDL(filename_H);
	
	//Find required quantities and transfer to matrix
	BOOST_FOREACH(boost::tie(header,data), boost::combine(headers_H, data_H)) {
		
		//Integrated elastic transport XS
		if (header[1][0]== 7 && header[1][1]==0) {
			transport_XS_H = data;
		}
		//Integrated elastic scattering XS
		else if (header[1][0]== 10 && header[1][1]==0) {
			total_scat_XS_H = data;
		}
		//Angular distribution large angle elastic scattering XS
		else if (header[1][0]== 8 && header[1][1]== 22) {
			scattering_XS_H = data;
		}	
	}
	
	//Read data for O
	[headers_O,data_O] = ReadENDL(filename_O);
	
	//Find required quantities and transfer to matrix
	BOOST_FOREACH(boost::tie(header,data), boost::combine(headers_O, data_O)) {
		
		//Integrated elastic transport XS
		if (header[1][0]== 7 && header[1][1]==0) {
			transport_XS_O = data;
		}
		//Integrated elastic scattering XS
		else if (header[1][0]== 10 && header[1][1]==0) {
			total_scat_XS_O = data;
		}
		//Angular distribution large angle elastic scattering XS
		else if (header[1][0]== 8 && header[1][1]== 22) {
			scattering_XS_O = data;
		}	
	}
	
	//Combine values for H and O according to mixture ratio in water 
    Physics::_xsTransportH2O = 0.11189400*transport_XS_H + 0.88810600*transport_XS_O;
	Physics::_xsH2O = 0.11189400*scattering_XS_H + 0.88810600*scattering_XS_O;
	Physics::_xsTotalH2O = 0.11189400*total_scat_XS_H + 0.88810600*total_scat_XS_O;
	Physics::_stpowH2O = 0.11189400*stopp_pow_XS_H + 0.88810600*stopp_pow_XS_O;
}

VectorVector Physics::GetScatteringXS (Vector energies,Vector density, Vector angle){

    VectorVector scattering_XS;
    return scattering_XS;
}


VectorVector Physics::GetTotalXS (Vector energies,Vector density) {
    VectorVector total_XS;
	
	if (energies < 1.000000000E-05) {
	
	}
	else if (energies > 1.000000000E+05) {
		
	}
	else {
	//Find upper and lower bound of energy value in table	
	CompressedVector<int>::Iterator upper( Physics::_xsTotalH2O(0).upperBound(energies));
	CompressedVector<int>::Iterator lower( Physics::_xsTotalH2O(0).lowerBound(energies));
	
	//Linear interpolation between upper and lower bound
	
	}
	
    return total_XS;
}

VectorVector Physics::GetStoppingPower( Vector energies ) {
    // @TODO
    VectorVector stopping_power;
    return stopping_power;
}


VectorVector Physics::GetTransportXS (Vector energies,Vector density){
    // @TODO
    VectorVector transport_XS;
    return transport_XS;
}

Physics* Physics::Create() { return new Physics(); }

std::tuple<std::list<VectorVector>,std::list<VectorVector>> Physics::ReadENDL( string filename ) {
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
    list<Vector> dataPackList;
    VectorVector headerPack( 2 );

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
	return {_headers,_data};
}

