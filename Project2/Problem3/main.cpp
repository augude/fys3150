#include "include/maxSymMatrix.hpp"
#include <iostream>
#include <armadillo>

int main(){
    // Initialize the matrix
    arma::mat A(4,4, arma::fill::eye);
    A(3, 0) = 0.5;
    A(0, 3) = 0.5;
    A(2, 1) = -0.7;
    A(1, 2) = -0.7;

    // Find the maximum entry in the symmetric matrix and its location
    int k, l;
    double max_off_diag = max_offdiag_symmetric(A, k, l);

    // Print the results
    std::printf("The maximum off-diagonal element is %.2f at (%d, %d)\n", max_off_diag, k, l);
    
    return 0;
}