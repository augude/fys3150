#include "../include/maxSymMatrix.hpp"

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){
    /*
    Loops through the upper triagonal of the symmetric matrix to find
    the maximum value.
    */
    //Determining the number of rows (symmetric matrix)
    int matsize=A.n_rows;

    //Consistency checks:
    assert(A.is_symmetric());
    assert(A.n_rows!=1);

    //setting initial comparison value
    double maxval=A(0,1); k=0;l=1;

    //Looping through the upper triagonal of the matrix
    for(int j=1; j<matsize; j++){
        for (int i=0;i<j;i++){
            if(maxval<A(i,j)){
                maxval=A(i,j); 
                k=i;
                l=j;
            }
        }
    }

    return maxval;
}