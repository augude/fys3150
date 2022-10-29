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
        double A = s.energyDiff[ed]; //acceptance probability from the energy difference using prerecorded values
        double r = ((double) rand() / (double) RAND_MAX);
        s.spins = s.spins * (r > A) + candidate * (r <= A); //accept the candidate if r <= A
    }
    return s;
}

void mcmc(double T, int L, int numberCycles, bool ordered, string filename){
    //vectors for storing quantities
    //vec energy = vec(numberCycles); 
    //vec energy1mom = vec(numberCycles);
    //vec energy2mom = vec(numberCycles);
    
   // vec magnetization  = vec(numberCycles);
   // vec mag1mom = vec(numberCycles);
   // vec mag2mom = vec(numberCycles);

    // vec heatCap = vec(numberCycles);
    // vec suscep = vec(numberCycles);

    std::ofstream ofile;
    ofile.open(filename);

    //initialize a lattice of size LxL at temperature T 
    Lattice s = Lattice(L, T, ordered);
    vec energyMag = s.energyMagnetization();
    double energy = energyMag(0);
    double energy1mom = energy;
    double energy2mom = energy*energy;

    double mag = energyMag(1);
    double mag1mom = abs(mag);
    double mag2mom = mag*mag;
    double heatCap = 1/(L*L*T*T)*(energy2mom - energy1mom*energy1mom);
    double suscep = 1/(L*L*T)*(mag2mom - mag1mom*mag1mom);
    
    ofile << "energy,energy1mom,energy2mom,magnetization,magnetization1mom,magnetization2mom,heatCapacity,susceptibility,temperature,gridsize" << endl;
    for (int i = 0; i < numberCycles; i++){
        Lattice nextState = mcCycle(s); //do one whole mc cycle 
        energyMag = nextState.energyMagnetization();
        energy = energyMag(0);
        energy1mom = (i*energy1mom + energy)/(i + 1);
        energy2mom = (i*energy2mom + energy*energy)/(i + 1);

        mag = energyMag(1);
        mag1mom = (i*mag1mom + abs(mag))/(i + 1);
        mag2mom = (i*mag2mom + mag*mag)/(i + 1);
        
        heatCap = 1/(L*L*T*T)*(energy2mom - energy1mom*energy1mom);
        suscep = 1/(L*L*T)*(mag2mom - mag1mom*mag1mom);

        //write values to file
        ofile << scientificFormat(energy) << "," 
              << scientificFormat(energy1mom) << ","
              << scientificFormat(energy2mom) << ","
              << scientificFormat(mag) << ","
              << scientificFormat(mag1mom) << ","
              << scientificFormat(mag2mom) << ","
              << scientificFormat(heatCap) << ","
              << scientificFormat(suscep) << ","
              << scientificFormat(T) << ","
              << scientificFormat(L) << endl;
        
        s.spins = nextState.spins; //update state and rerun 
    }
    ofile.close();
    
}
    
    