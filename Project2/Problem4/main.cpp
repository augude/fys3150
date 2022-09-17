#include "include/jacobi.hpp"
#include "../Problem3/include/maxSymMatrix.hpp"
#include "../include/utils.hpp"
#include "../Problem2/include/setupTridiag.hpp"

int main(){
    int N = 10;
    
    arma::mat A = arma::mat(N, N);
    setupTridiag(A);

    double eps = 1e-8;
    arma::vec eigenvalues = arma::vec(N);
    arma::mat eigenvectors = arma::mat(N, N);
    int maxiter = 1e5;
    int iterations = 0;
    bool converged = false;

    jacobiEigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
    
    std::map <double, arma::vec> exact = evalsExact(A);
    std::map <double, arma::vec> approx;
    
    for (int i = 0; i < N; i++){
        approx[eigenvalues(i)] = -1.0*eigenvectors.col(i) * (eigenvectors(i, 0) < 0) + 1.0*eigenvectors.col(i) * (eigenvectors(i, 0) >= 0);
    }

    testSetup(exact, approx);
    
    arma::vec eigval;
    arma::mat eigvec;

    eig_sym(eigval, eigvec, A);

    return 0;
}