#include "tridiagAlgo.hpp"
#include <fstream>
#include <string>
#include <math.h>
#include "utils.hpp"

int main(){
    int N = 10000; //number of steps 
    double h = 1.0/N; //step size
    double v0 = 0.0;  
    double vN = 0.0;
    
    std::vector<double> x;
    std::vector<double> g;

    for (double i = 0.0; i <= N; i++){
        //x takes steps values from 0 to 1. In total N + 1 points
        x.push_back(i/N);
    }
    for (int i = 1; i < N; i++){
        double value = h*h*100*exp(-10*x.at(i));
        g.push_back(value);
    }

    std::vector<double> a(N - 1, -1.0);
    std::vector<double> b(N - 1, 2);
    std::vector<double> c(N - 1, -1.0);

    std::vector <double> v = tridiagAlgo(a, b, c, g);

    std::string filename = "tridiagAlgo10000.txt";
    std::ofstream ofile;
    ofile.open(filename);

    ofile << scientificFormat(x) << std::endl;
    ofile << scientificFormat(v0);
    ofile << scientificFormat(v);
    ofile << scientificFormat(vN) << std::endl;

    ofile.close();

    return 0;
}