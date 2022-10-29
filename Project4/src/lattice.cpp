#include "../include/lattice.hpp"

Lattice::Lattice(int L_, double T_, bool ordered){
    T = T_;
    L = L_;
    
    if (ordered){
        mat lattice = mat(L, L);
        lattice.fill(1);
        spins = lattice;
    }
    else{
        arma_rng::set_seed_random();
        mat lattice = randi<mat>(L, L, arma::distr_param(0, 1));
        for (int i = 0; i < L ; i++){
            for (int j = 0; j < L; j++){
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
    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            int spin = spins(i, j); //spin at position (i, j)
            int spinRight = spins(i, (j + 1) % L); //spin to the right of (i, j) with periodic boundary conditions
            int spinUp = spins((i + 1) % L, j); //spin above (i, j) with periodic boundary conditions
            EM(0) += -spin*spinRight - spin*spinUp; 
            EM(1) += spin;
        }
    }
    return EM;
}

int Lattice::energyij(int i, int j){
    int E = 0;
    E -= spins(i, (j + 1) % L); //bound to neighbour to the right
    E -= spins(i, (j - 1 + L) % L); //bound to neighbour to the ledt
    E -= spins((i + 1) % L, j); //bound to neighbour above 
    E -= spins((i - 1 + L) % L, j); //bound to neighbour below
    return E*spins(i, j);
}
