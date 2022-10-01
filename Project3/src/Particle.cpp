#include "../include/Particle.hpp"


Particle::Particle(double chargeIn, double massIn, arma::vec positionIn, arma::vec velocityIn){

    charge = chargeIn;
    mass = massIn;
    position = positionIn;
    velocity = velocityIn;
}

void Particle::printCurrentPos(){

    std::cout << scientificFormat(position);
}

void Particle::printCurrentVel(){

    std::cout << scientificFormat(velocity);
}
