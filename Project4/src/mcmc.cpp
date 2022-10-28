#include "../include/mcmc.hpp"

Lattice mcmc(Lattice s){
    int N = s.N;
    arma::mat candidate = s.spins;
    int randomi = (int) rand() % N; //random integer indidies 
    int randomj = (int) rand() % N;
    candidate(randomi, randomj) = -1*candidate(randomi, randomj); //generate the candidate state

    //TODO: CAN PROBABLY SKIP THIS STEP SOMEHOW BASED ON THE RESULT OF PROBLEM 2
    int ed = -2*s.energyij(randomi, randomj); //calculate the energy difference between the two states
    double relativeProb = s.energyDiff[ed]; //get the relative probability from the energy difference using prerecorded values
    double A = min(1.0, relativeProb); //acceptance probability 
    double r = ((double) rand() / (double) RAND_MAX);
    if (r <= A){
        s.spins = candidate; //accept the candidate
    }
    return s;
}