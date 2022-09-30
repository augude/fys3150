#include "../include/PenningTrap.hpp"

PenningTrap::PenningTrap(double B0In, double V0In, double dIn){
    B0_ = B0In;
    V0_ = V0In;
    d_ = dIn;
}

arma::vec PenningTrap::electricField(arma::vec & position){
    arma::vec E = arma::vec(3);
    E(0) = position(0);
    E(1) = position(1);
    E(2) = -2*position(2);
    E = E*V0_/(2*d_*d_);
    return E;
}

arma::vec PenningTrap::magneticField(arma::vec & postion){
    arma::vec B = arma::vec(3);
    B(0) = 0;
    B(1) = 0;
    B(2) = B0_;
    return B;
}