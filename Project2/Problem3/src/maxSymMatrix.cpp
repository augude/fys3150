#include "../include/setupTridiag.hpp"

double max_offdiag_symmetric(const arma::mat& A, int& row, int &col){
    /*
    Loops through the upper triagonal of the symmetric matrix to find
    the maximum value.
    */
    //Determining the number of rows (symmetric matrix)
    int matsize=A.n_rows;

    //Consistency checks:

    //setting initial comparison value
    double maxval=A(0,1); row=0;col=1;
    

    for(int i =0; i<matsize-1; i++){
        for (int j=i+1;j<matsize;j++){
            if(maxval<A(i,j)){
                maxval=A(i,j); row=i;col=j;
            }
        }
    }

    return maxval;
}