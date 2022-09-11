#pragma once
#include <iostream>
#include <armadillo>
#include <map>
#include <math.h>

//function returns a maps of pairs of eigenvalues (doubles) and eienvectors (arma::vec)
std::map <double, arma::vec> evalsArma(int n);
std::map <double, arma::vec> evalsExact(int n);
void testSetup(int n);