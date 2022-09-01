#include "tridiagAlgo.hpp"

std::vector <double> tridiagAlgo(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> g){
    int N = b.size() + 1;
    std::vector <double> bPrime; 
    std::vector <double> gPrime;
    bPrime.push_back(b.at(0)); //initial values 
    gPrime.push_back(g.at(0)); 
    for (int i = 1; i < N - 1; i++){
        double bP = b.at(i) - a.at(i)*c.at(i - 1)/bPrime.at(i - 1);
        double gP = g.at(i) - a.at(i)*gPrime.at(i - 1)/bPrime.at(i - 1);
        bPrime.push_back(bP);
        gPrime.push_back(gP);
    }

    std::vector <double> v(N - 1);
    v.at(N - 2) =  gPrime[N - 2]/bPrime[N - 2];
    for (int ip = 3; ip < N + 1; ip++){
        int i = N - ip;
        v.at(i) = (gPrime.at(i) - c.at(i)*v.at(i + 1))/bPrime.at(i);
    }

    return v;
}