#pragma once 
#include <armadillo>
#include "../../Project2/include/utils.hpp"

class Particle{
    
    public:
        double charge_;
        double mass_;
        arma::vec postion_;
        arma::vec velocity_;

        Particle(double chargeIn, double massIn, arma::vec postionIn, arma::vec velocityIn);

        //print current positon of particle to screen
        void printCurrentPos();

        //print current velocity of particle to screen
        void printCurrentVel();
};