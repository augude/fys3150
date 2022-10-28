#include "../include/utils.hpp"


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

std::string scientificFormat(const std::vector<int>& v, const int width){
    std::stringstream ss;
    for (int elem : v){
        ss << std::setw(width) << std::scientific << elem;
    }   
    return ss.str();
}

std::string scientificFormat(const arma::vec v, const int width, const int prec){
    std::stringstream ss;
    for (double elem : v){
        ss << scientificFormat(elem, width, prec);
    }   
    return ss.str();
}


std::string scientificFormat(const std::vector<std::vector<double>>& v, const int width, const int prec){
    std::stringstream ss;
    for (std::vector<double> elem : v){
        ss << scientificFormat(elem, width, prec) << std::endl;
    }   
    return ss.str();
}

void setupTridiag(arma::mat & A){

    int N = A.n_cols;

    double h = (1.0)/(N + 1);
    
    A(0, 0) = 2.0/(h*h);

    for (int i = 1; i < N; i ++){
            A(i, i) = 2.0/(h*h);
            A(i - 1, i) = -1.0/(h*h); 
            A(i, i - 1) = -1.0/(h*h);
    }
}