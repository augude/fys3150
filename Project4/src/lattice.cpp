#include "../include/lattice.hpp"

Lattice::Lattice(int L_, double T_, bool ordered){
    //throw an assertion if L is not even
    assert(L % 2 == 0);
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
    
    //saving the values for the acceptance probabilities 
    energyDiff[-8] = 1.0;
    energyDiff[-4] = 1.0;
    energyDiff[0] = 1.0;
    energyDiff[4] = exp(-4/T);
    energyDiff[8] = exp(-8/T);
    
}

vec Lattice::energyMagnetization(){
    vec EM = vec(2).fill(0); //vector to store energy and magnetization
    //run through every row
    for (int j = 0; j < L; j++){
        //storing values that don't depend on the inner loop
        vec colj = spins.col(j);
        vec colright = spins.col((j + 1) % L);
        vec colleft = spins.col((j - 1 + L) % L);
        //run through every second coloumn starting at rownumber % 2
        for (int i = j % 2; i < L; i += 2){
            int pos = colj(i);
            int right = colright(i);
            int left = colleft(i);
            int up = colj((i + 1) % L);
            int down = colj((i - 1 + L) % L);
            EM(0) -= spins(i, j)*(right + left + up + down);
            EM(1) += spins(i, j)  + right;
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
