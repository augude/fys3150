#pragma once 
#include <armadillo>
#include <iostream>
#include <map>
#include <cassert>

class Lattice{
    
    public:
        int L; //size of lattice
        double T; //temperature
        int E; //energy of spin ensamble
        int M; //magnetization of spin ensamble
        arma::mat spins; //ensamble of spins
        std::map<int, double> energyDiff; //map from energy differences to acceptance probabilities
  
        
        //initialize the spin ensamble of size N x N with either ordered spins(all +1) or unordered(random +1/-1)
        Lattice(const int L_, double T_, bool ordered);

        //calculate initial energy and magnetization of lattice
        void energyMagnetization();
        
        //calculate the energy of spin (i, j)
        int energyij(int i, int j);
};