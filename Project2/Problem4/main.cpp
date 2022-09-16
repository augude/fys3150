#include "../include/jacobi.hpp"

int main()
{
    //Initialise example
    int N = 6;
    double h = 1.0/(N + 1);
    
    arma::mat A(N, N, arma::fill::zeros);
    for (int i = 1; i < N; i ++){
                A(i, i) = 2.0/(h*h);
                A(i - 1, i) = -1.0/(h*h);
                A(i, i - 1) = -1.0/(h*h);
    }
    
    arma::vec eigval(N);
    arma::mat eigvec(N, N, arma::fill::zeros);
    
    double eps = 1e-08;
    const int maxiter = 1000;
    int iterations = 0;
    bool converged = false;
    
    jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iterations, converged);
    
    return 0;
}
