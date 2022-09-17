#include "../include/jacobi.hpp"
#include "../../Problem3/include/maxSymMatrix.hpp"
//#include "../../include/utils.hpp"

void jacobiRotate(arma::mat& A, arma::mat& R, int k, int l, double eps){
    int N = A.n_cols;

    //extract elements of A
    double A_ll = A(l, l);
    double A_kk = A(k, k);
    double A_kl = A(k, l);

    double tau = (A_ll - A_kk)/(2.0*A_kl);

    //compute tridigonometric values
    double t = (1.0/(tau + sqrt(1 + tau*tau))) * (tau > 0) + (-1.0/(-tau + sqrt(1 + tau*tau))) * (tau <= 0);
    double c = 1.0/sqrt(1 + t*t);
    double s = c*t;

    //check for nan-values
    assert(isfinite(t));
    assert(isfinite(s));
    assert(isfinite(c));
    
    A(k, k) = A_kk*c*c - 2*A_kl*c*s + A_ll*s*s;
    A(l, l) = A_ll*c*c + 2*A_kl*c*s + A_kk*s*s;
    A(k, l) = 0;
    A(l, k) = 0;

    for (int i = 0; i < N; i++){
        if (i != k && i != l){
            double A_ki = A(k, i);
            double A_li = A(l, i);

            A(i, k) = A_ki*c - A_li*s;
            A(k, i) = A(i, k);
            A(i, l) = A_li*c + A_ki*s;
            A(l, i) = A(i, l);
        }
        double R_ik = R(i, k);
        double R_il = R(i, l);
        R(i, k) = R_ik*c - R_il*s;
        R(i, l) = R_il*c + R_ik*s;
    }
}

void jacobiEigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged){

    int N = A.n_cols;
    
    //make a copy of a non-constant A
    arma::mat Ac = A;

    //make identity matrix R
    arma::mat R = arma::mat(N, N).fill(0.0);
    for (int i = 0; i < N; i++){
        R(i, i) = 1.0;
    }

    double offDiagMax;
    int k;
    int l;

    offDiagMax = max_offdiag_symmetric(Ac, k, l);

    while (offDiagMax > eps && iterations < maxiter){
        jacobiRotate(Ac, R, k, l, eps);
        offDiagMax = max_offdiag_symmetric(Ac, k, l);
        iterations ++;
    }

    if (iterations < maxiter){
        converged = true;
        std::cout << "The algorithm converged after " << iterations << " iterations." << std::endl;
    }
    else{
        std::cout << "The algorithm did not converge after " << iterations << " iterations." << std::endl;
    }

    //extract the eigenvectors and eigenvalues
    eigenvectors = normalise(R);
    eigenvalues = arma::diagvec(Ac);
    
    //sort eigenvectors
    arma::uvec sorting_indices = arma::sort_index(eigenvalues);
    arma::mat temp = eigenvectors;
    for (int i=0; i<N; i++)
    {
        int old_index = sorting_indices(i);
        eigenvectors.col(i) = temp.col(old_index);

    }
    //sort eigenvalues
    eigenvalues = arma::sort(eigenvalues);
}

