#include "../include/setupTridiag.hpp"
#include "../../include/utils.hpp"

std::map <double, arma::vec> evalsArma(int N){
    /*
        takes in the dimension of the matrix
        the dimension of A is two less than the number of points in v
        and one less than the number of steps
    */
    double h = 1.0/(N + 1); //h is defined by the number of steps, not the matrix dim
    
    arma::mat A = arma::mat(N, N).fill(0.0);
    A(0, 0) = 2.0/(h*h);

    for (int i = 1; i < N; i ++){
            A(i, i) = 2.0/(h*h);
            A(i - 1, i) = -1.0/(h*h); 
            A(i, i - 1) = -1.0/(h*h);
    }
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, A);
    
    std::map <double, arma::vec> evalsEvec;

    for (int i = 0; i < eigval.size(); i++){
        evalsEvec[eigval(i)] = -1*eigvec.col(i)*(eigvec.col(i)(0) < 0) + eigvec.col(i)*(eigvec.col(i)(0) > 0);
    }

    return evalsEvec;
}

std::map <double, arma::vec> evalsExact(int N){
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

void testSetup(int N){
    std::map <double, arma::vec> armaEvalsEvec = evalsArma(N);
    std::map <double, arma::vec> exactEvalsEvec = evalsExact(N);
    arma::vec armaEvals = arma::vec(N);
    arma::mat armaEvecs = arma::mat(N, N).fill(0.0);
    arma::vec exactEvals = arma::vec(N);
    arma::mat exactEvecs = arma::mat(N, N).fill(0.0);

    for(auto it = exactEvalsEvec.cbegin(); it != exactEvalsEvec.cend(); ++it){
        int i = 0;
        exactEvals(i) = it -> first;
        exactEvecs.col(i) = it -> second;
        i += 1;
    }

    for(auto it = armaEvalsEvec.cbegin(); it != armaEvalsEvec.cend(); ++it){
        int i = 0;
        armaEvals(i) = it -> first;
        armaEvecs.col(i) = it -> second;
        i += 1;
    }
    double s = 0.0;
    for (int i = 0; i < N; i++){
        s += (armaEvals(i) - exactEvals(i));
        s += sum(armaEvecs.col(i) - exactEvecs.col(i));
    }
    if (abs(s) < 1e-6){
        std::cout << "The test passed. The difference between the analytical and the numerical values is: " << scientificFormat(abs(s)) << std::endl;
    }
    else{
        std::cout << "The test failed. The difference between the analytical and the numerical values is: " << scientificFormat(abs(s)) << std::endl;
    }
}