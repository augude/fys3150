#include "../include/includefiles.hpp"
#include "jacobi.hpp"

int main()
{
    //Initialise example
    int N = 6;
    double h = 1.0/(N + 1);
    
    arma::mat A(N, N, arma::fill::zeros);
    A(0, 0) = 2.0/(h*h);
    for (int i = 1; i < N; i ++){
                A(i, i) = 2.0/(h*h);
                A(i - 1, i) = -1.0/(h*h);
                A(i, i - 1) = -1.0/(h*h);
    }
    
    //Initialise eigenvalues and eigenvectors for analytical soltuion
    arma::vec eigval(N);
    arma::mat eigvec(N, N, arma::fill::zeros);
    
    //For rotation solution
    arma::vec eigval1(N);
    arma::mat eigvec1(N, N, arma::fill::zeros);
    
    double eps = 1e-10;
    const int maxiter = 1000;
    int iterations = 0;
    bool converged = false;
    
    //Analytical solution
    eig_sym(eigval, eigvec, A);
    
    //Jacobi solution
    jacobi_eigensolver(A, eps, eigval1, eigvec1, maxiter, iterations, converged);
    
    //Comparing with analytical values
    std::cout << "Analytical eigenvalues:" << std::endl;
    std::cout << eigval << std::endl;
    std::cout << "Jacobi eigenvalues:" << std::endl;
    std::cout << eigval1 << std::endl;
    std::cout << "Analytical eigenvectors:" << std::endl;
    std::cout << eigvec << std::endl;
    std::cout << "Jacobi eigenvectors:" << std::endl;
    std::cout << eigvec1 << std::endl;
    
    return 0;
}
