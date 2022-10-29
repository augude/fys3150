#include "../include/mcmc.hpp"
#include "../include/utils.hpp"

Lattice mcCycle(Lattice s){
    int L = s.L;
    for (int i = 0; i < L*L; i++){
        arma::mat candidate = s.spins;
        int randomi = (int) rand() % L; //random integer indidies 
        int randomj = (int) rand() % L;
        candidate(randomi, randomj) = -1*candidate(randomi, randomj); //generate the candidate state

        int ed = -2*s.energyij(randomi, randomj); //calculate the energy difference between the two states
        double relativeProb = s.energyDiff[ed]; //get the relative probability from the energy difference using prerecorded values
        double A = min(1.0, relativeProb); //acceptance probability 
        double r = ((double) rand() / (double) RAND_MAX);
        if (r <= A){
            s.spins = candidate; //accept the candidate
        }
    }
    return s;
}

void mcmc(double T, int L, int numberCycles, string filename){
    //vectors for storing quantities
    vec energy = vec(numberCycles); 
    vec energy1mom = vec(numberCycles);
    vec energy2mom = vec(numberCycles);
    
    vec magnetization  = vec(numberCycles);
    vec mag1mom = vec(numberCycles);
    vec mag2mom = vec(numberCycles);

    vec heatCap = vec(numberCycles);
    vec suscep = vec(numberCycles);

    std::ofstream ofile;
    ofile.open(filename);

    //initialize an onordered 2x2 lattice at temperature T 
    Lattice s = Lattice(L, T, false);
    
    ofile << "energy,energy1mom,energy2mom,magnetization,magnetization1mom,magnetization2mom,heatCapacity,susceptibility,temperature,gridsize" << endl;
    for (int i = 0; i < numberCycles; i++){
        Lattice nextState = mcCycle(s); //do one whole mc cycle 
        energy(i) = nextState.energyMagnetization()(0);
        energy1mom(i) = arma::sum(energy)/(i + 1);
        energy2mom(i) = arma::sum(arma::square(energy)) /(i + 1);
        magnetization(i) = nextState.energyMagnetization()(1);
        mag1mom(i) = arma::sum(arma::abs(magnetization))/(i + 1);
        mag2mom(i) = arma::sum(arma::square(magnetization))/(i + 1);
        heatCap = 1/(L*L*T*T)*(energy2mom - arma::square(energy1mom));
        suscep = 1/(L*L*T)*(mag2mom - arma::square(mag1mom));

        //write values to file
        ofile << scientificFormat(energy(i)) << "," 
              << scientificFormat(energy1mom(i)) << ","
              << scientificFormat(energy2mom(i)) << ","
              << scientificFormat(magnetization(i)) << ","
              << scientificFormat(mag1mom(i)) << ","
              << scientificFormat(mag2mom(i)) << ","
              << scientificFormat(heatCap(i)) << ","
              << scientificFormat(suscep(i)) << ","
              << scientificFormat(T) << ","
              << scientificFormat(L) << endl;
        
        s.spins = nextState.spins; //update state and rerun 
    }
    ofile.close();
    
}
    
    