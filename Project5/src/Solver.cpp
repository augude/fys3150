#include "../include/Solver.hpp"

Solver::Solver(int M_in, int N_in, double DT){
    M = M_in;
    N = N_in;
    dt = DT;
    h = 1./(M-1);
    r = std::complex<double>(0, dt/(2*h*h));

    /* u.zeros(m, m);
    v.zeros(m, m);
    int z = 1;

    // Initalise input matrix V with some constant
    for (int i = 1; i < m-1; i++){
        for (int j = 1; j < m-1; j++){
            v(i, j) = z;
        }
    } */

    // Only internal points
    u.zeros(M-2, M-2);
    v.ones(M-2, M-2); //Initalise input matrix V with some constant
}

int Solver::pair_to_single(int i, int j){
    return j*(M-2) + i;
}

void Solver::fill_matrices(std::complex<double> r, arma::cx_vec a, arma::cx_vec b){
    //matrices
    A.zeros((M-2)*(M-2), (M-2)*(M-2));
    B.zeros((M-2)*(M-2), (M-2)*(M-2));

    //submatrices
    arma::sp_cx_mat subDiag(M-2, M-2); subDiag.zeros();
    
    //diagonals
    a.zeros((M-2)*(M-2));
    b.zeros((M-2)*(M-2));

    //The r diagonals
    arma::cx_vec offDiag((M-2)*(M-3)); offDiag.fill(r);
    
    //Creating subDiag:
    arma::cx_vec r_vec(M-3); r_vec.fill(r);
    subDiag.diag(1)=r_vec;
    subDiag.diag(-1)=r_vec;

    arma::cx_double c_b;
    int k;

    //Creating the diagonals. 
    for (int i = 0; i < M-2; i++){
        for (int j = 0; j < M-2; j++){
            k = pair_to_single(i, j);
            std::complex<double> c_a(1.0, 4.0*r.imag() + 0.5*dt*v.at(i,j));
            c_b = std::conj(c_a);

            a.at(k) = c_a;
            b.at(k) = c_b;
        }
    }

    std::cout << a << std::endl;
    std::cout << b << std::endl;


    //Inserting into A and B
    A.diag(M-2)=-offDiag; A.diag(2-M)=-offDiag;
    B.diag(M-2)=offDiag; B.diag(2-M)=offDiag;

    for (int i=0;i<M-2;i++){
        int j=i*(M-2);
        A.submat(j,j,j+M-3,j+M-3)=-subDiag;
        B.submat(j,j,j+M-3,j+M-3)=subDiag;
    }

    A.diag() = a;
    B.diag() = b;
    std::cout << A << std::endl;
    std::cout << B << std::endl;
}


arma::cx_mat Solver::forward(arma::cx_mat A, arma::cx_mat B, arma::vec u)
        {
            arma::cx_vec b = B*u;
            arma::cx_vec x = arma::solve(A, b);
            return x;
        }
