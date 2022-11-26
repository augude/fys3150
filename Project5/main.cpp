#include "include/Solver.hpp"
#include <string>
using namespace std; 
using namespace arma;

int main(int argc, char* argv[])
{
    if (argc != 11){
        string executable_name = argv[0];
        cerr << "Error: Wrong number of input arguments." << std::endl;
        cerr << "Usage: " << executable_name << " <double: h> <double: dt> <double: T> <double: xc> <double: sigx> <double: px> <double: yc> <double: sigy> <double: py> <string: type>" << std::endl;
        return 1;
    }
    else{
        double h = atof(argv[1]); 
        double dt = atof(argv[2]); 
        double T = atof(argv[3]);
        double xc = atof(argv[4]); 
        double sigx = atof(argv[5]);
        double px = atof(argv[6]); 
        double yc = atof(argv[7]);   
        double sigy = atof(argv[8]); 
        double py = atof(argv[9]);
        string type = argv[10];

        Solver obj(h, dt, T);
        obj.set_initial_state(xc, sigx, px, yc, sigy, py);
        obj.set_potential("input/" + type + ".dat");
        obj.fill_matrices();
        obj.solve();
        if (T > 5e-3){
            obj.states.save("output/" + type + "_long_evolution.bin");
        }
        else{
            obj.states.save("output/" + type + "_evolution.bin");
        }
        
        return 0;
    }
}