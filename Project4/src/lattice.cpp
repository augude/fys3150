#include "../include/lattice.hpp"
using namespace std;
using namespace arma; 

Lattice::Lattice(const int L_, double T_, bool ordered){
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
        for (int j = 0; j < L ; j++){
            for (int i = 0; i < L; i++){
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

void Lattice::energyMagnetization(){
    //run through every row
    for (int j = 0; j < L; j++){
        //storing values that don't depend on the inner loop
        vec colj = spins.col(j);
        vec colright = spins.col((j + 1) % L);
        //run through every second coloumn starting at rownumber % 2
        for (int i = 0; i < L; i ++){
            int pos = colj(i);
            int right = colright(i);
            int up = colj((i + 1) % L);
            E -= spins(i, j)*(right + up);
            M += spins(i, j);
        }
    }
}

int Lattice::energyij(int i, int j){
    int s = spins(i, j);
    int left = spins(i, (j + 1) % L); //neighbour to the right
    int right = spins(i, (j - 1 + L) % L); //neighbour to the left
    int up = spins((i + 1) % L, j); //neighbour above 
    int down = spins((i - 1 + L) % L, j); //neighbour below
    int E = -s*(left + right + up + down);
    return E;
}
