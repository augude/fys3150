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

        // initial states values
        double xc;
        double sigx;
        double px;
        double yc;
        double sigy;
        double py;

        // Constructor
        Solver(double h_in, double dt_in, double T_in, double xc_in, double sigx_in, double px_in, double yc_in, double sigy_in, double py_in);

        // Converting matrix indices to vector index
        int pair_to_single(int i, int j);

        // Filling matrices for Crank-Nicholson
        void fill_matrices();


        // Sets the initial state matrix U
        void set_initial_state();

        // Sets the potential matrix V
        void set_potential(std::string filename);

        // Moves forward one time step
        void forward();

        // Solves for all time steps with the option of a measurement at measure_time with the result (x, y)
        void solve(bool measure = false, double measure_time = 0, double x = 0, double y = 0);

        //measures the position of the particle to be (x, y)
        void measure_particle(double x, double y);

};