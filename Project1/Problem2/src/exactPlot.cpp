#include "exactPlot.hpp"
#include "utils.hpp"

std::tuple<std::vector<double>, std::vector<double>> calculateExact(int steps){
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

void printToFile(int steps){
    
    auto [x, u] = calculateExact(steps);    
    std::string filename = "exactPlot.txt";
    std::ofstream ofile;
    ofile.open(filename);

    ofile << scientificFormat(x) << std::endl;
    ofile << scientificFormat(u) << std::endl;

    ofile.close();
}