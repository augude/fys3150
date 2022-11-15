#include <iostream>
#include <armadillo>
#include <complex>
#include <cmath>

class Solver 
{
    public:
        int m, n; // Iter points for x,y axes and time
        double dt;
        double h;

        arma::mat u;
        arma::mat v;

        arma::cx_mat A;
        arma::cx_mat B;

        std::complex<double> r;

        /* Constructor */
        Solver(int M, int N, double DT)
        {
            m = M;
            n = N;
            dt = DT;
            h = 1./(N-1);
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
            u.zeros(m-2, m-2);
            v.ones(m-2, m-2); //Initalise input matrix V with some constant

        }

        /* Converting double index to single index */
        int pair_to_single(int i, int j)
        {
            return j*(m-2) + i;
        }

        /* Filling matrices for Crank-Nicholson */
        void fill_matrices()
        {
            A.zeros((m-2)*(m-2), (m-2)*(m-2));
            B.zeros((m-2)*(m-2), (m-2)*(m-2));
            arma::cx_vec a((m-2)*(m-2));
            arma::cx_vec b((m-2)*(m-2));

            std::complex<double> c_b;
            int k;

            for (int i = 0; i < m-2; i++){
                for (int j = 0; j < m-2; j++){
                    k = pair_to_single(i, j);
                    std::complex<double> c_a(1.0, 4.0*r.imag() + 0.5*dt*v.at(i,j));
                    c_b = std::conj(c_a);

                    a.at(k) = c_a;
                    b.at(k) = c_b;
                }
            }

            std::cout << a << std::endl;
            std::cout << b << std::endl;

            A.diag() = a;
            B.diag() = b;


            for (int i = (m-2); i < (m-2)*(m-2); i++){
                A.at(i - (m-2), i) = -r;
                A.at(i, i - (m-2)) = -r;

                B.at(i - (m-2), i) = r;
                B.at(i, i - (m-2)) = r;
            }

            for (int i = 1; i < (m-2)*(m-2); i++){
                if (i % (m-2) != 0){
                    A.at(i-1, i) = -r;
                    A.at(i, i-1) = -r;

                    B.at(i-1, i) = r;
                    B.at(i, i-1) = r;
                }
            }
        }

        arma::cx_mat forward(arma::cx_mat A, arma::cx_mat B, arma::vec u)
        {
            arma::cx_vec b = B*u;
            arma::cx_vec x = arma::solve(A, b);
            return x;
        }

};

int main()
{
    Solver obj(5, 5, 0.1);
    obj.fill_matrices();

    return 0;
}