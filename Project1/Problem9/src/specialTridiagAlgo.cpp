#include "../include/specialTridiagAlgo.hpp"

std::vector <double> specialTridiagAlgo(std::vector<double> g){
    int N = g.size();
    std::vector <double> gPrime;
    gPrime.push_back(g.at(0)); 
    for (int i = 1; i < N; i++){
        double dev = i + 1.0;
        double gP = g.at(i) + (i)/dev*gPrime.at(i - 1);
        gPrime.push_back(gP);
    }

    std::vector <double> v(N);
    double dev = 1.0*N;
    v.at(N - 1) =  gPrime[N - 1]*(N - 1)/dev;
    for (int ip = 2; ip < N + 1; ip++){
        int i = N - ip;
        double dev = i + 0.0;
        v.at(i) = (gPrime.at(i) + v.at(i + 1))*(i + 1)/(dev + 2);
    }
    return v;
}