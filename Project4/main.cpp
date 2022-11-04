#include "include/lattice.hpp"
#include "include/mcmc.hpp"
#include "include/utils.hpp"
#include <vector>
#include <time.h>
#include "omp.h"  //OpenMP header
using namespace std;
using namespace arma; 
using namespace chrono;


int main(int argc, char* argv[]){
    srand(time(NULL)); //set random seed
    string testString  = argv[1];

    if (testString == "validation2x2"){
        //N samples from a 2x2 grid in order to compare with analytical results 
        int N = 1e6; //number of MC cycles
        const int L = 2; //size of lattice
        double T = 1.0;
        string filename = "output/validation2x2.csv";
        mcmc(T, L, N, false, filename);
    }

    else if (testString == "burnIn"){
        //N samples from a 20x20 grid at different initial conditions and at different temperatures to investigate burn-in time
        int N = 1e6; //number of MC cycles
        const int L = 20; //size of lattice
        double T1 = 1.0;
        double T2 = 2.4;
        string filename1 = "output/ordered1.csv";
        string filename2 = "output/ordered2.csv";
        string filename3 = "output/unordered1.csv";
        string filename4 = "output/unordered2.csv";

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

    else if (testString == "histograms"){
        //N samples from a 20x20 grid at 5 different temperatures to compare distributions of epsilon
        int N = 1e6; //number of MC cycles
        const int L = 20; //size of lattice
        vec temp = arma::linspace<vec>(1.0, 2.4, 5);
        string filename;

        #pragma omp parallel
        {
            #pragma omp for
            for (double t : temp){
                filename = "output/L=20_" + to_string(omp_get_thread_num()) + ".csv";
                mcmc(t, L, N, false, filename); 
            } //end parallelized loop
        } //end parallelized region
    }

    else if (testString == "par"){
        // 1 runs with a total of 1e7 samples from 2x2 grid with different number of threads (1 - 4) to find speed-up factor 
        // the total number of samples is distributed evenly across the different threads
        double T = 1.0; //temperature
        int L = 20;
        string filename = "output/test.csv";
        int numberThreads;
        auto t1 = high_resolution_clock::now();
        #pragma omp parallel
        {
            numberThreads = omp_get_num_threads();
            int N = (double) 1e6/numberThreads;  //number steps per cycle
            mcmc(T, L, N, false, filename);
        } //end parallelized region  
        auto t2 = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(t2 - t1);
        cout << scientificFormat(numberThreads) << scientificFormat(duration.count()) << endl; 
        
    }
    else if (testString == "problem8"){
        int N = 1e6; //number of MC cycles
        const int L = atoi(argv[2]); //size of lattice
        vec temp = arma::linspace<vec>(2.1, 2.4, 10);
        string filename;
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < 10; i++){
                filename = "output/L=" + to_string(L) + "_" + to_string(i) + ".csv";
                mcmc(temp[i], L, N, true, filename); 
            } //end parallelized loop
        } //end parallelized region
    }
        
    return 0;
}