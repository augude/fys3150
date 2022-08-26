#include "utils.hpp"


std::string scientific_format(const double d, const int width = 20, const int prec = 10){
    std::stringstream ss; 
    ss << std::setw(width) << std::setprecision(prec) << std::scientific;
    return ss.str();

}
