#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <chrono>
#include "../Problem7/include/tridiagAlgo.hpp"
#include "../Problem2/include/utils.hpp"
#include "include/timers.hpp"

int main(){
    int N = 10;

    std::vector<double> general;
    std::vector<double> special;

    std::vector<int> ns = std::vector<int>{10, 100, 1000, 10000, 100000, 1000000};
    double tg;
    double ts;
    for (int n: ns){
        tg = 0.0;
        ts = 0.0;
        for (int i=0; i<N; i++){
            tg += time_general(n);
            ts += time_special(n);
        }
        general.push_back(tg / N);
        special.push_back(ts / N);
    }

    std::string filename = "timings.txt";  

    std::ofstream ofile;
    ofile.open(filename);

    ofile << scientificFormat(general) << std::endl;
    ofile << scientificFormat(special) << std::endl;

    ofile.close();

    return 0;
}
