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
	Physics::Physics::_stpowH2O = ReadStoppingPowers(fileName_stppower);
}

VectorVector Physics::GetScatteringXS (Vector energies,Vector density, Vector angle){

    VectorVector scattering_XS;
    return scattering_XS;
}


VectorVector Physics::GetTotalXS (Vector energies,Vector density) {
    VectorVector total_XS;
	double xsH2O_i;
	
	for (int i=0; i< energies.size(); i++){
		
		if (energies(i) < Physics::_xsTotalH2O(1)(0)) {
		xs_H2O_i = Physics::_xsTotalH2O(1)(0);
		}
		else if (energies(i) > Physics::_xsTotalH2O(1)(energies.size()-1)) {
		xs_H2O_i = 	_xsTotalH2O(1)(energies.size()-1);
		}
		else {
		//Find upper and lower bound of energy value in table	
		CompressedVector<int>::Iterator upper( Physics::_xsTotalH2O(0).upperBound(energies(i)));
		CompressedVector<int>::Iterator lower( Physics::_xsTotalH2O(0).lowerBound(energies(i)));
		
		//Linear interpolation between upper and lower bound
		xsH2O_i = Physics::_xsTotalH2O(1)((lower - Physics::_xsTotalH2O(0).begin())) + (energies(i) - Physics::_xsTotalH2O(0)((lower - Physics::_xsTotalH2O(0).begin()))) * (Physics::_xsTotalH2O(1)((upper - Physics::_xsTotalH2O(0).begin())) - Physics::_xsTotalH2O(1)((lower - Physics::_xsTotalH2O(0).begin()))) / (Physics::_xsTotalH2O(0)((upper - Physics::_xsTotalH2O(0).begin())) - Physics::_xsTotalH2O(0)((lower - Physics::_xsTotalH2O(0).begin())));
		}
		
		totalXS(i) = xsH2O_i*density;
	
	}
	
    return total_XS;
}

VectorVector Physics::GetStoppingPower(Vector energies,Vector density) {
    VectorVector stopping_power;
	double stpw_H2O_i;
	
	for (int i=0; i< energies.size(); i++){
		
		if (energies(i) < Physics::_stpowH2O(1)(0)) {
		stpw_H2O_i = Physics::_stpowH2O(1)(0);
		}
		else if (energies(i) > Physics::_stpowH2O(1)(energies.size()-1)) {
		stpw_H2O_i = Physics::_stpowH2O(1)(energies.size()-1);
		}
		else {
		//Find upper and lower bound of energy value in table	
		CompressedVector<int>::Iterator upper( Physics::_stpowH2O(0).upperBound(energies(i)));
		CompressedVector<int>::Iterator lower( Physics::_stpowH2O(0).lowerBound(energies(i)));
		
		//Linear interpolation between upper and lower bound
		stpw_H2O_i = Physics::_stpowH2O(1)((lower - Physics::_stpowH2O(0).begin())) + (energies(i) - Physics::_stpowH2O(0)((lower - Physics::_stpowH2O(0).begin()))) * (Physics::_stpowH2O(1)((upper - Physics::_stpowH2O(0).begin())) - Physics::_stpowH2O(1)((lower - Physics::_stpowH2O(0).begin()))) / (Physics::_stpowH2O(0)((upper - Physics::_stpowH2O(0).begin())) - Physics::_stpowH2O(0)((lower - Physics::_stpowH2O(0).begin())));
		}
		
		stopping_power(i) = stpw_H2O_i*density;
	
	}
    return stopping_power;
}


VectorVector Physics::GetTransportXS (Vector energies,Vector density){
    VectorVector transport_XS;
	double xsH2O_i;
	
	for (int i=0; i< energies.size(); i++){
		
		if (energies(i) < Physics::_xsTransportH2O(1)(0)) {
		xs_H2O_i = Physics::_xsTransportH2O(1)(0);
		}
		else if (energies(i) > Physics::_xsTransportH2O(1)(energies.size()-1)) {
		xs_H2O_i = 	Physics::_xsTransportH2O(1)(energies.size()-1);
		}
		else {
		//Find upper and lower bound of energy value in table	
		CompressedVector<int>::Iterator upper( Physics::_xsTransportH2O(0).upperBound(energies(i)));
		CompressedVector<int>::Iterator lower( Physics::_xsTransportH2O(0).lowerBound(energies(i)));
		
		//Linear interpolation between upper and lower bound
		xsH2O_i = Physics::_xsTransportH2O(1)((lower - Physics::_xsTransportH2O(0).begin())) + (energies(i) - Physics::_xsTransportH2O(0)((lower - Physics::_xsTransportH2O(0).begin()))) * (Physics::_xsTransportH2O(1)((upper - Physics::_xsTransportH2O(0).begin())) - Physics::_xsTransportH2O(1)((lower - Physics::_xsTransportH2O(0).begin()))) / (Physics::_xsTransportH2O(0)((upper - _Physics::_xsTransportH2O(0).begin())) - Physics::_xsTransportH2O(0)((lower - Physics::_xsTransportH2O(0).begin())));
		}
		
		transport_XS(i) = xsH2O_i*density;
	
	}
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
        ErrorMessages::Error( "The ENDL file is missing!!", CURRENT_FUNCTION );
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

VectorVector Physics::ReadStoppingPowers(std::string fileName) {
	VectorVector stp_powers;
	string text_line;
	
	ifstream file_stream;
	file_stream.open( fileName, ios::in );
	
	if( file_stream.fail() ) {
		ErrorMessages::Error( "The Stopping Power file is missing!!", CURRENT_FUNCTION );
    }
	
	while( getline(file_stream, text_line ) ) { 
	
	    std::stringstream ss( text_line );    // give line to stringstream
        list<double> linelist;
        for( double double_in; ss >> double_in; ) {    // parse line
            linelist.push_back( double_in );
        }
	
		stp_powers[0].append(linelist.front());   //append elements from line to column vectors 
		linelist.pop_front();
		stp_powers[1].append(linelist.front());
		linelist.pop_front();	
        
    }
    case_file.close();
}

