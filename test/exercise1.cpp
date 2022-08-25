#include "exercise1.hpp"

std::vector<double> expo(std::vector<double> x){
    std::vector<double> y;
    for (int i = 0; i < x.size(); i++){
        y.push_back(exp(x[i]));
    }

    return y;
}