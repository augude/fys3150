#pragma once 
#include <armadillo>
#include <iostream>
using namespace std;
using namespace arma; 

class Lattice{
    
    public:
        int N; //size of lattice
        double T; //temperature
        mat spins; //ensamble of spins    

        //initialize the spin ensamble of size N x N with either ordered spins(all +1) or unordered(random +1/-1)
        Lattice(int N_, double T_, bool ordered);

        //calculate energy and magnetization of lattice
        vec energyMagnetization();
        
};