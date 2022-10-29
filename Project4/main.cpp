#include "include/lattice.hpp"
#include "include/mcmc.hpp"
#include "include/utils.hpp"

int main(){
    srand(time(NULL));
    int N = 2000; //number of MC cycles
    int L = 2; //size of lattice
    double T = 1.0;
    string filename = "validation2x2.csv";
    mcmc(T, L, N, filename);
    return 0;
}