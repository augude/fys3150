#pragma once 
#include <iostream>
#include <armadillo>

class PenningTrap{
    public:
        double B0_;
        double V0_;
        double d_;

        PenningTrap(double B0In, double V0In, double dIn);

        arma::vec electricField(arma::vec & position);

        arma::vec magneticField(arma::vec & position);
};