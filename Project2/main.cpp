#include "../include/includefiles.hpp"
#include "jacobi.hpp"

int main()
{
    //Initialise a random example
    int N = 6;
    
    const arma::mat A(N, N, arma::fill::randu);
    
    arma::vec eigval(N);
    arma::mat eigvec(N, N, arma::fill::zeros);
    
    const int maxiter = 100;
    
    jacobi_eigensolver(A, 1e-08, eigval, eigvec, maxiter, 0, false);
    
}
