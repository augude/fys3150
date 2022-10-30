#include "../include/mcmc.hpp"
#include "../include/utils.hpp"

void mcCycle(Lattice& s, double& E, double& M){
    int L = s.L;
    for (int i = 0; i < L*L; i++){
        int randomi = (int) rand() % L; //random integer indidies 
        int randomj = (int) rand() % L;

        int ed = -2*s.energyij(randomi, randomj); //energy difference between the two states
        double A = s.energyDiff[ed]; //acceptance probability from the energy difference using prerecorded values
        double r = ((double) rand() / (double) RAND_MAX);
        
        if (r <= A){
            //update the state, energy and magnetization
            s.spins(randomi, randomj) *= -1; 
            int md = -2*s.spins(randomi, randomj); //magnatization difference between the two states
            E += ed;
            M += md;
        }
    }
}

void mcmc(double T, int L, int numberCycles, bool ordered, string filename){
   
    std::ofstream ofile;
    ofile.open(filename);

    //initialize a lattice of size LxL at temperature T 
    Lattice s = Lattice(L, T, ordered);

    double denomCv= L*L*T*T; //denominator in expression for Cv
    double denomChi = L*L*T; //denominator in expression for chi
    vec energyMag = s.energyMagnetization(); //initial energy and magnetization

    //values to be stored
    double energy = energyMag(0);
    double energy1mom = energy;
    double energy2mom = energy*energy;

    double mag = energyMag(1);
    double mag1mom = abs(mag);
    double mag2mom = mag*mag;
    double heatCap = 1/denomCv*(energy2mom - energy1mom*energy1mom);
    double suscep = 1/denomChi*(mag2mom - mag1mom*mag1mom);
    
    ofile << "energy,energy1mom,energy2mom,magnetization,magnetization1mom,magnetization2mom,heatCapacity,susceptibility,temperature,gridsize" << endl;
    for (int i = 0; i < numberCycles; i++){
        mcCycle(s, energy, mag); //do one whole mc cycle and update the state, energy and magnetization 
        energy1mom = (i*energy1mom + energy)/(i + 1);
        energy2mom = (i*energy2mom + energy*energy)/(i + 1);

        mag1mom = (i*mag1mom + abs(mag))/(i + 1);
        mag2mom = (i*mag2mom + mag*mag)/(i + 1);
        
        heatCap = 1/denomCv*(energy2mom - energy1mom*energy1mom);
        suscep = 1/denomChi*(mag2mom - mag1mom*mag1mom);

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
    }
    ofile.close();
    
}
    
    