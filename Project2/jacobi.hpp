#ifndef jacobi_hpp
#define jacobi_hpp

#include "../include/includefiles.hpp"

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigval,
                        arma::mat& eigvec, const int maxiter, int& iterations, bool& converged);

#endif /* jacobi_hpp */
