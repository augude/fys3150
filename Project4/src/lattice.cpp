#include "../include/lattice.hpp"

Lattice::Lattice(int N_, double T_, bool ordered){
    T = T_;
    N = N_;
    
    if (ordered){
        mat lattice = mat(N, N);
        lattice.fill(1);
        spins = lattice;
    }
    else{
        arma_rng::set_seed_random();
        mat lattice = randi<mat>(N, N, arma::distr_param(0, 1));
        for (int i = 0; i < N ; i++){
            for (int j = 0; j < N; j++){
                if (lattice(i, j) == 0){
                    lattice(i, j) = -1;
                }
            }
        }
        spins = lattice;
    }
    
    //saving the values for relative probabilities 
    energyDiff[-8] = exp(8/T);
    energyDiff[-4] = exp(4/T);
    energyDiff[0] = 1.0;
    energyDiff[4] = exp(4/T);
    energyDiff[8] = exp(-8/T);
    
}

vec Lattice::energyMagnetization(){
    vec EM = vec(2).fill(0); //vector to store energy and magnetization
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            int spin = spins(i, j); //spin at position (i, j)
            int spinRight = spins(i, (j + 1) % N); //spin to the right of (i, j) with periodic boundary conditions
            int spinUp = spins((i + 1) % N, j); //spin above (i, j) with periodic boundary conditions
            EM(0) += -spin*spinRight - spin*spinUp; 
            EM(1) += abs(spin);
        }
    }
    return EM;
}

int Lattice::energyij(int i, int j){
    int E = 0;
    E -= spins(i, (j + 1) % N); //bound to neighbour to the right
    E -= spins(i, (j - 1 + N) % N); //bound to neighbour to the ledt
    E -= spins((i + 1) % N, j); //bound to neighbour above 
    E -= spins((i - 1 + N) % N, j); //bound to neighbour below
    return E*spins(i, j);
}
