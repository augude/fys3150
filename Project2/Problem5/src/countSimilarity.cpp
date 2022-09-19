#include "../include/countSimilarity.hpp"
#include "../../Problem4/include/jacobi.hpp"
#include "../../include/utils.hpp"

std::vector<int> countSimilarities(const std::vector<int> listOfN, double eps){
    /*
        Takes in a list of values of N
        Does the Jacobi rotation algo for each N and records the number of iterations in a vector.
    */
   bool converged = false;
   std::vector<int> listOfIter;
   for (int N : listOfN){
        arma::mat A = arma::mat(N, N); 
        arma::vec eigenvalues = arma::vec(N);
        arma::mat eigenvectors = arma::mat(N, N).fill(0.0);
        int maxiter = 1e5; //expects convergence in O(N^4)
        int iterations = 0;
        setupTridiag(A);
        // Generate random N*N matrix
        //arma::mat A = arma::mat(N, N).randn();  
        // Symmetrize the matrix by reflecting the upper triangle to lower triangle
        //A = arma::symmatu(A); 
        jacobiEigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
        listOfIter.push_back(iterations);
   }
   
    return listOfIter;
}