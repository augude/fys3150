#pragma once
#include <iostream>
#include <armadillo>
#include <complex>
#include <cmath>
#include <fstream>
#include <assert.h>

class Solver{

    public:
        int M, N;
        double T; 
        double dt;
        double h;

        // State matrix and vector equivalent
        arma::cx_mat U;
        arma::cx_vec u;
        arma::cx_mat tmp;

        arma::mat V; // Input matrix (assumed real?)

        // Crank-Nicholson matrices
        arma::sp_cx_mat A;
        arma::sp_cx_mat B;
        std::complex<double> r;

        arma::cx_cube states; // State storage cube


        // Constructor
        Solver(double h_in, double dt_in, double T_in);

        // Converting matrix indices to vector index
        int pair_to_single(int i, int j);

        // Filling matrices for Crank-Nicholson
        void fill_matrices();

        // Sets the initial state matrix U
        void set_initial_state(double xc, double sigx, double px, double yc, double sigy, double py);

        // Sets the potential matrix V
        void set_potential(std::string filename, double v0);

        // Moves forward one time step
        void forward();

        // Solves for all time steps
        void solve();

        // Writes all state matrices 
        void write_to_file();

};