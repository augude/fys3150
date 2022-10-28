#include "include/lattice.hpp"
#include "include/mcmc.hpp"
#include "include/utils.hpp"

int main(){
    srand(time(NULL));
    int N = 1000; //number of steps in Markov chain
    int Nl = 2; //size of lattice
    double T = 1.0;
    vec energies = vec(N + 1);
    vec magnetizations = vec(N + 1);
    
    Lattice s = Lattice(Nl, T, false);
    vec EM = s.energyMagnetization();
    energies(0) = EM(0);
    magnetizations(0) = EM(1);
    
    for (int i = 0; i < N; i++){
        Lattice s1 = mcmc(s);
        vec EM = s1.energyMagnetization();
        energies(i) = EM(0);
        magnetizations(i) = EM(1);
        s = s1; //update the state
    }
    double Z = 12 + 4*cosh(8);
    double teoMean = -32*sinh(8)/Z;
    cout << "Sample average of energies after " << N << " steps in Markov chain: " << endl;
    std::cout << scientificFormat(arma::mean(energies)) << std::endl;
    cout << "Theoretical expected value: " << endl;
    std::cout << scientificFormat(teoMean) << std::endl;
    return 0;
}