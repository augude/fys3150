#include "../include/mcmc.hpp"
using namespace std;
using namespace arma; 


void mcCycle(Lattice& s, int& E, int& M, vec randomi, vec randomj, vec r){
    int L = s.L;
    for (int i = 0; i < L*L; i++){
        
        int ed = -2*s.energyij(randomi[i], randomj[i]); //energy difference between the two states
        double A = s.energyDiff[ed]; //acceptance probability from the energy difference using prerecorded values
        
        if (r[i] <= A){
            //update the state, energy and magnetization
            s.spins(randomi[i], randomj[i]) *= -1; 
            int md = -2*s.spins(randomi[i], randomj[i]); //magnatization difference between the two states
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
        vec randomi = randi<vec>(L*L, arma::distr_param(0,L-1));
        vec randomj = randi<vec>(L*L, arma::distr_param(0,L-1));
        vec r(L*L, fill::randu);

        mcCycle(s, energy, mag, randomi, randomj, r); //do one whole mc cycle and update the state, energy and magnetization 
        //write values to file
        ofile << energy << "," 
              << mag << endl;
    }
    ofile.close();
    
}
    
    