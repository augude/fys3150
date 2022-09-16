#pragma once 
#include <iostream>
#include <armadillo>
#include <math.h>

void jacobiRotate(arma::mat& A, arma::mat& R, int k, int l, double eps);

void jacobiEigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged);