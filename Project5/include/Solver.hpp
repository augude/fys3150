#include <iostream>
#include <armadillo>
#include <complex>
#include <cmath>

class Solver{

    public:
        int M, N; // Iter points for x,y axes and time
        double dt;
        double h;
        arma::mat u;
        arma::mat v;
        arma::cx_mat A;
        arma::cx_mat B;
        std::complex<double> r;
        arma::cx_vec a;
        arma::cx_vec b;

        //Constructor
        Solver(int M_in, int N_in, double T);

        //Converting double index to single index
        int pair_to_single(int i, int j);

        //Filling matrices for Crank-Nicholson
        void fill_matrices(std::complex<double> r, arma::cx_vec a, arma::cx_vec b);

        //Forwards the simulation one time step
        arma::cx_mat forward(arma::cx_mat A, arma::cx_mat B, arma::vec u);




    



};