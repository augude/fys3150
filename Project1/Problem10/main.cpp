#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <chrono>
#include "../Problem7/include/tridiagAlgo.hpp"
#include "../Problem2/include/utils.hpp"
#include "include/timers.hpp"

int main(){
    int N = 50;

    std::vector<std::vector<double>> general;
    std::vector<std::vector<double>> special;


    std::vector<int> ns = std::vector<int>{10, 100, 1000, 10000, 100000, 1000000};

    std::vector<double> ts(N);
    std::vector<double> tg(N);

    for (int n: ns){
        for (int i=0; i<N; i++){
            tg.at(i) = time_general(n);
            ts.at(i) = time_special(n);
        }
        general.push_back(tg);
        special.push_back(ts);
    }

    std::string filenameG = "timingsGeneral.txt";  
    std::string filenameS = "timingsSpecial.txt";

    std::ofstream ofileS;
    std::ofstream ofileG;
    ofileG.open(filenameG);
    ofileS.open(filenameS);

    ofileG << scientificFormat(general) << std::endl;
    ofileS << scientificFormat(special) << std::endl;

    ofileG.close();
    ofileS.close();

    return 0;
}
