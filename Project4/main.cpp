#include "include/lattice.hpp"
#include "include/mcmc.hpp"
#include "include/utils.hpp"
#include <vector>
#include "omp.h"  //OpenMP header
using namespace std;
using namespace arma; 


int main(int argc, char* argv[]){
    srand(time(NULL)); //set random seed
    string testString  = argv[1];

    if (testString == "validation2x2"){
        int N = 1e6; //number of MC cycles
        const int L = 2; //size of lattice
        double T = 1.0;
        string filename = "validation2x2.csv";
        mcmc(T, L, N, false, filename);
    }

    else if (testString == "burnIn"){
        int N = 1e6; //number of MC cycles
        const int L = 20; //size of lattice
        double T1 = 1.0;
        double T2 = 2.4;
        string filename1 = "ordered1.csv";
        string filename2 = "ordered2.csv";
        string filename3 = "unordered1.csv";
        string filename4 = "unordered2.csv";

        vector<double> temp = {T1, T2, T1, T2}; 
        vector<string> filenames = {filename1, filename2, filename3, filename4};
        vector<bool> order = {true, true, false, false};

        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < 4; i++){
                mcmc(temp[i], L, N, order[i], filenames[i]); 
            } //end parallelized loop
        } //end parallelized region
    }
    return 0;
}