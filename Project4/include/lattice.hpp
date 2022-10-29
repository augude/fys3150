#pragma once 
#include <armadillo>
#include <iostream>
#include <map>
#include <cassert>
using namespace std;
using namespace arma; 

class Lattice{
    
    public:
        int L; //size of lattice
        double T; //temperature
        mat spins; //ensamble of spins
        map<int, double> energyDiff; //map from energy differences to acceptance probabilities
  
        
        //initialize the spin ensamble of size N x N with either ordered spins(all +1) or unordered(random +1/-1)
        Lattice(const int L_, double T_, bool ordered);

        //calculate energy and magnetization of lattice
        vec energyMagnetization();
        
        //calculate the energy of spin (i, j)
        int energyij(int i, int j);
};