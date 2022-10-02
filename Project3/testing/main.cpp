#include "../include/PenningTrap.hpp"
#include "../include/Particle.hpp"
#include "../../Project2/include/utils.hpp"
#include "include/tests.hpp"

int main(int argc, char* argv[]){
    std::string testString  = argv[1];
    
    if (testString == "PenningSetup"){
        testPenningSetup();
    }
    else if (testString == "OneParticleFE"){
        testOneParticleFE();
    }
    else if (testString == "OneParticleRK4"){
        testOneParticleRK4();
    }
    else if (testString == "DoubleSetupWithInternal"){
        testDoubleSetup(true);
    }
    else if (testString == "DoubleSetupWithoutInternal"){
        testDoubleSetup(false);
    };
    return 0;
}