#include "include/countSimilarity.hpp"
#include "../include/utils.hpp"
#include <math.h>

int main(){

    double eps = 1e-8;
    std::vector<int> listOfN;
    for (int i = 2; i < 101; i++){
        listOfN.push_back(i);
    }
    // vector with values (2, 3, ..., 100)
    std::vector<int> listOfIter = countSimilarities(listOfN, eps);

    std::string filename = "countSimilarities.txt";
    std::ofstream ofile;
    ofile.open(filename);

    ofile << scientificFormat(listOfN) << std::endl;
    ofile << scientificFormat(listOfIter) << std::endl;

    ofile.close();

    return 0;
}