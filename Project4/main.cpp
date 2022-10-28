#include "include/lattice.hpp"

int main(){
    for (int i = 0; i < 10; i++){
    int N = 2;
    double T = 1.0;
    Lattice test = Lattice(N, T, false);
    vec EM = test.energyMagnetization();
    cout << "Energy: " << EM(0) << " Magnetizaion: " << EM(1) << endl;
    }
    return 0;
}