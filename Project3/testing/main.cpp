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
    }
    else if (testString == "CompareStepsize"){
        double stepSize = atof(argv[2]);
        compareStepsize(stepSize);
    }
    else if (testString == "fractionWithin"){
        for (double f : {0.1, 0.4, 0.7}){
            arma::vec fraction = fractionWithin(f);
            std::cout << scientificFormat(fraction) << std::endl;
        }
    }
    else if (testString == "fractionWithinZoom"){
        arma::vec fractionWithoutCouloumb = fractionWithinZoom(false);
        std::cout << scientificFormat(fractionWithoutCouloumb) << std::endl;
        arma::vec fractionWithCouloumb = fractionWithinZoom(true);
        std::cout << scientificFormat(fractionWithCouloumb) << std::endl;
    }
    return 0;
}