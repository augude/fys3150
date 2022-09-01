#include "utils.hpp"


std::string scientificFormat(double d, const int width, const int prec){
    std::stringstream ss; 
    ss << std::setw(width) << std::setprecision(prec) << std::scientific << d;
    return ss.str();
}

std::string scientificFormat(const std::vector<double>&v, const int width, const int prec){
    std::stringstream ss;
    for (double elem : v){
        ss << scientificFormat(elem, width, prec);
    }   
    return ss.str();
}

