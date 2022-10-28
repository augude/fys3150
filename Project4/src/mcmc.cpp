#include "../include/mcmc.hpp"

Lattice mcmc(Lattice s){
    //TODO: WHERE TO BEST DECLARE THESE?
    map<int, double> energyDiff;
    energyDiff[-8] = exp(8);
    energyDiff[-4] = exp(4);
    energyDiff[0] = 1.0;
    energyDiff[4] = exp(4);
    energyDiff[8] = exp(-8);
    srand(time(NULL)); //initialize random seed
    int N = s.N;
    Lattice candidate = s;
    int randomi = (int) rand() % N; //random indidies 
    int randomj = (int) rand() % N;
    candidate.spins(randomi, randomj) = -1*candidate.spins(randomi, randomj); //generate the candidate state
    //TODO: CAN PROBABLY SKIP THIS STEP SOMEHOW BASED ON THE RESULT OF PROBLEM 2
    int ed = -2*s.energyij(randomi, randomj); //calculate the energy difference between the two states
    double relativeProb = energyDiff[ed]; //get the relative probability from the energydifference using prerecorded values
    double A = min(1.0, relativeProb); //acceptance probability 
    double r = (double) rand();
}