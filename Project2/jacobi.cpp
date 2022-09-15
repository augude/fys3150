#include "jacobi.hpp"
#include "max_offdiag.hpp"

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l)
{
    double c;
    double s;
    
    //Concerns all non-zero offdiagonal elements
    if (A(k, l) != 0.0){
        double t;
        double the;

        //Choosing theta such that all offdiagonal elements become zero
        the = (A(l, l) - A(k, k))/(2*A(k, l));

        //Must choose t to be the smaller of the two roots
        //Smaller root corresponds to a rotation angle smaller than pi/4
        //Using the form of the quadratic formular with the discriminant in the denominator, the smaller root can be written as
        if (the > 0){
            t = 1.0/(the + sqrt(the*the + 1));
        } else {
            t = -1.0/(the + sqrt(the*the + 1));
        }
            

        c = 1.0/sqrt(t*t + 1);
        s = c*t;

    } else {
        c = 1.0;
        s = 0.0;
    }

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k, k);
    a_ll = A(l, l);

    A(k, k) = c*c*a_kk - 2.0*c*s*A(k, l) + s*s*a_ll;
    A(l, l) = s*s*a_kk + 2.0*c*s*A(k, l) + c*c*a_ll;
    
    //Hard-coding zeros
    A(k, l) = 0.0;
    A(l, k) = 0.0;

    for (int i = 0; i < A.n_rows; i++) {
        if ( i != k && i != l ) {
            a_ik = A(i, k);
            a_il = A(i, l);
            A(i, k) = c*a_ik - s*a_il;
            A(k, i) = A(i, k);
            A(i, l) = c*a_il + s*a_ik;
            A(l, i) = A(i, l);
        }
    
        //Compute the new eigenvectors
        r_ik = R(i, k);
        r_il = R(i, l);
        R(i, k) = c*r_ik - s*r_il;
        R(i, l) = c*r_il + s*r_ik;
        
    }
    return;
}

void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigval,
                        arma::mat& eigvec, const int maxiter, int& iterations, bool& converged)
{
    //Initialise S as A and R as identity matrix
    int N = A.n_rows;
    arma::mat S = A;
    arma::mat R(N, N, arma::fill::eye);
    
    //Define k, l values for max off-diagonal
    int k, l;
    double maxoffdiag = max_offdiag(S, k, l);
    
    //Rotate until:
    //All max off-diagonal values become zero
    //Max number of iterations is reached
    //Convergence is reached
    while (fabs(maxoffdiag) > eps && iterations < maxiter && converged == false){
        maxoffdiag = max_offdiag(S, k, l);
        jacobi_rotate(S, R, k, l);
        
        //Consider the sum of squares of the off-diagonal elements to check for convergence of the Jacobi method
        double sum = 0.0;
        for (int a = 0; a < N; a++){
            sum += pow(arma::accu(S.col(a)) - S(a, a), 2);
            
            if (abs(sum - 2*pow(maxoffdiag, 2)) < eps){
                converged = true;
            }
        
        iterations++;
        }
    }

    //Write out iterations
    std::cout << "The number of iterations performed was " << iterations << "." << std::endl;
    
    //Write eigenvalues and eigenvectors as entries
    for (int j = 0; j < N; j++){
        eigval(j) = S(j, j);
        eigvec.col(j) = S.col(j);
    }
    
    std::cout << eigval << std::endl;
    std::cout << eigvec << std::endl;
    
    
    return;
}