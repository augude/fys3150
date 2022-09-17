#pragma once
#include <iostream>
#include <armadillo>
#include <math.h>
#include <map>
#include <assert.h>

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, double eps);
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigval,
                        arma::mat& eigvec, const int maxiter, int& iterations, bool& converged);
