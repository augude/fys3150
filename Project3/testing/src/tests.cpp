#include "../../include/PenningTrap.hpp"
#include "../../../Project2/include/utils.hpp"
#include <fstream>
#include <string>


void testPenningSetup(){

    double B0 = 9.51e1;
    double V0 = 9.65e8;
    double d = 1e4;

    std::string filename = "ElectricField.txt";
    std::ofstream ofile;
    ofile.open(filename);

    arma::vec pos(3);
    arma::vec eField(3);
    PenningTrap trap = PenningTrap(B0, V0, d);
    int N = 100;
    for (int z = 0; z < N; z++){
        for (int y = 0; y < N; y++){
            double zpos = 2*d*z/N - d;
            double ypos = 2*d*y/N - d;
            pos(0) = 0;
            pos(1) = ypos;
            pos(2) = zpos;
            double V = V0/(2*d*d)*(2*zpos*zpos - ypos*ypos - 0*0);
            eField = trap.electricField(pos);  
            ofile << scientificFormat(z) << scientificFormat(y);
            ofile << scientificFormat(pos); 
            ofile << scientificFormat(eField) << scientificFormat(V) << std::endl;
        }
    }
    ofile.close();
}