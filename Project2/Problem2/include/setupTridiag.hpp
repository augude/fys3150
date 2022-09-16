#pragma once
#include <iostream>
#include <armadillo>
#include <map>
#include <math.h>

//function returns a maps of pairs of eigenvalues (doubles) and eienvectors (arma::vec)
std::map <double, arma::vec> evalsArma(const arma::mat & A);
std::map <double, arma::vec> evalsExact(const arma::mat & A);
void testSetup(std::map <double, arma::vec> exact, std::map <double, arma::vec> approx);