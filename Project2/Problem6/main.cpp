#include "../include/utils.hpp"
#include "../Problem4/include/jacobi.hpp"
#include "../Problem3/include/maxSymMatrix.hpp"
#include "../Problem2/include/setupTridiag.hpp"
#include <vector>

int main(){
    std::vector <int> Ns = {9, 99};
    for (int N : Ns){
        arma::mat A = arma::mat(N, N);
        setupTridiag(A);
        double eps = 1e-3;
        arma::mat eigenvectors = arma::mat(N, N);
        arma::vec eigenvalues = arma::vec(N);
        int maxiter = 1e5;
        int iterations = 0;
        bool converged = false;
        jacobiEigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

        std::map <double, arma::vec> exact = evalsExact(A);
        std::string filename1 = std::string("eigenvectors") + std::to_string(N) + std::string(".txt");
        std::ofstream ofile1;
        ofile1.open(filename1);
        for (int i = 0; i < N; i++){
            ofile1 << scientificFormat(eigenvectors.col(i)) << std::endl;
        }
        
        for (auto it = exact.cbegin(); it != exact.cend(); ++it){
            ofile1 << scientificFormat(it -> second) << std::endl;
        }
        ofile1.close();

        std::string filename2 = std::string("eigenvalues") + std::to_string(N) + std::string(".txt");
        std::ofstream ofile2;
        ofile2.open(filename2);
        for (int i = 0; i < N; i++){
            ofile2 << scientificFormat(eigenvalues(i)) << std::endl;
        }

        for (auto it = exact.cbegin(); it != exact.cend(); ++it){
            ofile2 << scientificFormat(it -> first) << std::endl;
        }
        ofile2.close();
    }

    return 0;
}