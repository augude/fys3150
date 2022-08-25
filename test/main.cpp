#include "exercise1.hpp"
#include <fstream>
#include <string>

int main(){
    double n = 100.0;
    std::vector<double> x;
    for (int i = 0; i < n; i ++){
        x.push_back(i/n);
    }
    
    std::vector<double> y = expo(x);

    std::string filename = "output.txt";

    std::ofstream ofile;
    ofile.open(filename);

    for (int i = 0; i < n; i++){
        ofile << x[i] << " " << y[i] << "\n"  << std::endl; 
    }
    ofile.close();

    return 0;
}