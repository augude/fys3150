#include "excactPlot.hpp"

std::tuple<std::vector<double>, std::vector<double>> calculateExact(int steps, int prec){
    std::vector<double> x;
    std::vector<double> u;

    for (double i = 0.0; i <= steps; i++){
        //x takes steps values from 0 to 1
        x.push_back(i/steps);
        //x.back() returns last element in the vector.  
        double value = 1 - (1 - exp(-10))*x.back() - exp(-10*x.back());
        u.push_back(value);
    }

    return {x, u};
}

void printToFile(int steps, int prec){
    
    auto [x, u] = calculateExact(steps, prec);    
    std::string filename = "exactPlot.txt";
    std::ofstream ofile;
    ofile.open(filename);

    for (int i = 0; i < steps; i++){
        ofile << std::setprecision(prec) << std::scientific << x[i] << " " << u[i] << "\n"  << std::endl; 
    }

    ofile.close();
}