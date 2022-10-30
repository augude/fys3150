#include "../include/mcmc.hpp"

void mcCycle(Lattice& s, int& E, int& M){
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

    vec energyMag = s.energyMagnetization(); //initial energy and magnetization

    //values to be stored
    int energy = energyMag(0);
    int mag = energyMag(1);
    
    ofile << "energy,magnetization,temperature,gridsize" << endl;
    ofile << energy << "," 
              << mag << ","
              << T << ","
              << L << endl;
    for (int i = 0; i < numberCycles; i++){
        mcCycle(s, energy, mag); //do one whole mc cycle and update the state, energy and magnetization 
        //write values to file
        ofile << energy << "," 
              << mag << endl;
    }
    ofile.close();
    
}
    
    