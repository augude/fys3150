#pragma once 
#include <armadillo>
#include "../../Project2/include/utils.hpp"

class Particle{
    
    public:
        double charge;
        double mass;
        arma::vec position;
        arma::vec velocity;

        Particle(double chargeIn, double massIn, arma::vec positionIn, arma::vec velocityIn);

        //print current positon of particle to screen
        void printCurrentPos();

        //print current velocity of particle to screen
        void printCurrentVel();
};