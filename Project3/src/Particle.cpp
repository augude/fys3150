#include "../include/Particle.hpp"


Particle::Particle(double chargeIn, double massIn, arma::vec postionIn, arma::vec velocityIn){

    charge_ = chargeIn;
    mass_ = massIn;
    postion_ = postionIn;
    velocity_ = velocityIn;
}

void Particle::printCurrentPos(){

    std::cout << scientificFormat(postion_);
}

void Particle::printCurrentVel(){

    std::cout << scientificFormat(velocity_);
}
