#include "physics.h"

Physics::Physics (Settings* settings){
    // @TODO set up physical model from settings?
    // @TODO What else belongs in physics class?
    // load H2O cross sections
    LoadXSH2O("test1","test2");
}

void Physics::LoadXSH2O( std::string fileName1, std::string fileName2){
    // @TODO
}

VectorVector Physics::GetScatteringXS (std::vector<double> energies,std::vector<double> density, std::vector<double> Omegas){
    // @TODO
}

VectorVector Physics::GetTotalXS (std::vector<double> energies,std::vector<double> density){
    // @TODO 
}

VectorVector Physics::GetStoppingPower (std::vector<double> energies,std::vector<double> sH2O){
    // @TODO
}

Physics* Physics::Create( Settings* settings ) { return new Physics( settings ); }