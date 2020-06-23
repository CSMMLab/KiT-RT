#include "physics.h"
using namespace std;
using blaze::CompressedVector;


Physics::Physics( Config* settings ) {
    // @TODO set up physical model from settings?
    // @TODO What else belongs in physics class?
    // load H2O cross sections
    LoadXSH2O( "test1", "test2" );
}

Physics::Physics (){
    //LoadXSH2O("test1","test2");
    //LoadTransportXSH2O("test3","test4");

}

void Physics::LoadXSH2O( std::string fileName1, std::string fileName2 ) {
    // @TODO
}

void Physics::LoadDatabase( std::string fileName){
    VectorVector transport_XS_H;
    VectorVector transport_XS_O;
	
	VectorVector scattering_XS_H;
    VectorVector scattering_XS_O;
	
	VectorVector total_scat_XS_H;
    VectorVector total_scat_XS_O;
	
	VectorVector stopp_pow_H;
    VectorVector stopp_pow_O;
	
    //Read data for H



    //Read data for O
	
	//Combine values for H and O according to mixture ratio in water 
    Physics::_xsTransportH2O = 0.11189400*transport_XS_H + 0.88810600*transport_XS_O;
	Physics::_xsH2O = 0.11189400*scattering_XS_H + 0.88810600*scattering_XS_O;
	Physics::_xsTotalH2O = 0.11189400*total_scat_XS_H + 0.88810600*total_scat_XS_O;
	Physics::_stpowH2O = 0.11189400*stopp_pow_XS_H + 0.88810600*stopp_pow_XS_O;
}

VectorVector Physics::GetScatteringXS (Vector energies,Vector density, Vector Theta){

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

VectorVector Physics::GetStoppingPower( std::vector<double> energies, std::vector<double> sH2O ) {
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
