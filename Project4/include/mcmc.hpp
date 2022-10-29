#pragma once 
#include <string>
#include <fstream>
#include "lattice.hpp"

//performe one cycle (corresponds to N attemped spin flips) and update the state, energy and magnetization
void mcCycle(Lattice& s, double& E, double& M);

//performe many mccycles at temperature T of LxL-grid and print results to csv file 
void mcmc(double T, int L, int numberCycles, bool ordered, string filename);