#include "include/setupTridiag.hpp"
#include "../include/utils.hpp"

int main(){
    int N = 6;
    arma::mat A = arma::mat(N, N).fill(0.0);
    setupTridiag(A);
    std::map <double, arma::vec> exact = evalsExact(A);
    std::map <double, arma::vec> arma = evalsArma(A);
    testSetup(exact, arma);
    
    return 0;
}