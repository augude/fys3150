#include "../include/setupTridiag.hpp"
#include "../../include/utils.hpp"

std::map <double, arma::vec> evalsArma(const arma::mat & A){
    /*
        takes in a reference to the matrix. 
    */
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, A);
    
    std::map <double, arma::vec> evalsEvec;

    for (int i = 0; i < eigval.size(); i++){
        evalsEvec[eigval(i)] = -1*eigvec.col(i)*(eigvec.col(i)(0) < 0) + eigvec.col(i)*(eigvec.col(i)(0) > 0);
    }

    return evalsEvec;
}

std::map <double, arma::vec> evalsExact(const arma::mat & A){
    int N = A.n_cols;
    double h = 1.0/(N + 1); //h is defined by the number of steps, not the matrix dim
    std::map <double, arma::vec> evalsEvec;

    for (int i = 0; i < N; i++){
        int indexI = i + 1;
        double eval = 2.0/(h*h) - 2.0/(h*h)*cos((indexI*M_PI)/(N + 1));
        arma::vec evec = arma::vec(N).fill(0.0);
        for (int j = 0; j < N; j++){
            int indexJ = j + 1;
            evec(j) = sin((indexI*indexJ*M_PI)/(N + 1));
        }
        //Want that the first element is positive
        evalsEvec[eval] = normalise(-1*evec*(evec(0) < 0) + evec*(evec(0) > 0)) ;
        
    }

    return evalsEvec;
}

void testSetup(std::map <double, arma::vec> exact, std::map <double, arma::vec> approx){
    int N = exact.size();
    arma::vec approxEvals = arma::vec(N);
    arma::mat approxEvecs = arma::mat(N, N).fill(0.0);
    arma::vec exactEvals = arma::vec(N);
    arma::mat exactEvecs = arma::mat(N, N).fill(0.0);

    for(auto it = exact.cbegin(); it != exact.cend(); ++it){
        int i = 0;
        exactEvals(i) = it -> first;
        exactEvecs.col(i) = it -> second;
        i += 1;
    }

    for(auto it = approx.cbegin(); it != approx.cend(); ++it){
        int i = 0;
        approxEvals(i) = it -> first;
        approxEvecs.col(i) = it -> second;
        i += 1;
    }
    double s = 0.0;
    for (int i = 0; i < N; i++){
        s += (approxEvals(i) - exactEvals(i));
        s += sum(approxEvecs.col(i) - exactEvecs.col(i));
    }
    if (abs(s) < 1e-6){
        std::cout << "The test passed. The difference between the analytical and the numerical values is: " << scientificFormat(abs(s)) << std::endl;
    }
    else{
        std::cout << "The test failed. The difference between the analytical and the numerical values is: " << scientificFormat(abs(s)) << std::endl;
    }
}