#pragma once 
#include <string>
#include <fstream>
#include <vector>
#include <chrono>
#include "lattice.hpp"

//performe one cycle (corresponds to N attemped spin flips) and update the state, energy and magnetization
void mcCycle(Lattice& s, std::vector<int> randomi, std::vector<int> randomj, std::vector<double> r);

//performe many mccycles at temperature T of LxL-grid and write results to csv file 
void mcmc(double T, int L, int numberCycles, bool ordered, std::string filename);