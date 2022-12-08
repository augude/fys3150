#include "../include/Solver.hpp"
#include "../include/progressbar.hpp"
using namespace std;
using namespace arma;

Solver::Solver(double h_in, double dt_in, double T_in, double xc_in, double sigx_in, double px_in, double yc_in, double sigy_in, double py_in)
{
    h = h_in;
    dt = dt_in;
    T = T_in;

    xc = xc_in;
    sigx = sigx_in;
    px = px_in;

    yc = yc_in;
    sigy = sigy_in;
    py = py_in;

    M = 1./h + 1;
    N = T/dt + 1;
    r = std::complex<double>(0, dt/(2*h*h));

    // To avoid issues with reshape, used advanced constructor for u vector
    // Any changes to matrix/vector form affects both
    U.zeros(M-2, M-2);
    u = arma::cx_vec(U.memptr(), U.n_elem, false, false);

    tmp.zeros(M-2, M-2);
    V.zeros(M-2, M-2);
    states.zeros(M-2, M-2, N);
}


int Solver::pair_to_single(int i, int j)
{
    return j*(M-2) + i;
}


void Solver::fill_matrices()
{
    cout << "Filling matrices ..." << endl;
    A.zeros((M-2)*(M-2), (M-2)*(M-2));
    B.zeros((M-2)*(M-2), (M-2)*(M-2));
    arma::cx_vec a((M-2)*(M-2), arma::fill::zeros);
    arma::cx_vec b((M-2)*(M-2), arma::fill::zeros);
    
    // Initalise and set size for submatrices
    arma::sp_cx_mat subDiag(M-2, M-2); 
    subDiag.zeros();
    
    // Creating r/-r diagonals
    arma::cx_vec rDiag((M-2)*(M-3)); 
    rDiag.fill(r);
    
    // Filling super/subdiagonals for submatrices with r/-r
    arma::cx_vec r_vec(M-3); 
    r_vec.fill(r);
    subDiag.diag(1) = r_vec;
    subDiag.diag(-1) = r_vec;

    // Filling diagonal vectors a and b
    arma::cx_double c_b;
    int k;

    for (int i = 0; i < M-2; i++){
        for (int j = 0; j < M-2; j++){
            k = pair_to_single(i, j);
            std::complex<double> c_a(1.0, 4.0*r.imag() + 0.5*dt*V.at(i,j));
            c_b = std::conj(c_a);

            a.at(k) = c_a;
            b.at(k) = c_b;
        }
    }

    // Inserting r/-r diagonals into A and B
    A.diag(M-2) = -rDiag;
    A.diag(2-M) = -rDiag;
    B.diag(M-2) = rDiag; 
    B.diag(2-M) = rDiag;

    // Inserting submatrices into A and B
    for (int i = 0; i < M-2; i++){
        int j = i*(M-2);
        A.submat(j, j, j + M-3, j + M-3) = -subDiag;
        B.submat(j, j, j + M-3, j + M-3) = subDiag;
    }

    // Inserting a and b vectors as diagonals for A and B
    A.diag() = a;
    B.diag() = b;
    cout << "Done" << endl;
}


void Solver::set_initial_state()
{
    double real;
    double imag;
    std::complex<double> e;

    for (int i = 0; i < M-2; i++){
        for (int j = 0; j < M-2; j++){
            real = -((i*h - xc)*(i*h - xc))/(2*sigx*sigx) - ((j*h - yc)*(j*h - yc))/(2*sigy*sigy);
            imag = px*(i*h - xc) + py*(j*h - yc);
            e = std::complex<double>(real, imag);
            U.at(i,j) = exp(e);
        }
    }

    // Elementwise multiplication and normalization
    U = U/sqrt(arma::accu(arma::conj(U) % U)); 
}


void Solver::forward()
{
    tmp = u;
    // Tridiag solver ideal
    arma::spsolve(u, A, B*tmp);
}


void Solver::measure_particle(double x, double y){
    sigx = sigx/4;
    sigy = sigy/4;
    xc = x; 
    yc = y;
    set_initial_state();
}


void Solver::solve(bool measure, double measure_time, double x, double y)
{
    int measure_index = N*measure_time/T;
    cout << "Solving eq ..." << endl;
    progressbar bar(N);
    for (int i = 0; i < N; i++){
        states.slice(i) = U;
        forward();        
        bar.update();
        if (i == measure_index && measure){
            measure_particle(x, y);
            cout << " A wild measurement occurred! Your particle is at (" << xc << ", " << yc << ")." << endl;         
        }
    }
}


void Solver::set_potential(std::string filename)
{
    cout << "Loading potential ..." << endl;
   V.load(filename, arma::raw_ascii);
   bool ok = V.load(filename, arma::raw_ascii);
   if (ok == false){
    cout << "Problem with loading file in line: " << __LINE__  << "in file: " << __FILE__ << endl;
   }
    assert(V.n_rows == (M - 2) && V.n_cols == (M - 2));
    cout << "Done" << endl;
}
