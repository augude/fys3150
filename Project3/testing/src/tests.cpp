#include "../../include/PenningTrap.hpp"
#include "../../../Project2/include/utils.hpp"
#include "../../include/Particle.hpp"
#include <fstream>
#include <string>


void testPenningSetup(){

    double B0 = 9.51e1;
    double V0 = 9.65e8;
    double d = 1e4;

    arma::vec particlePos;
    arma::vec particleVec;
    double charge;
    double mass;

    Particle part = Particle(charge, mass, particlePos, particleVec);

    std::vector<Particle> parts = {part};
    std::string filename = "ElectricField.txt";
    std::ofstream ofile;
    ofile.open(filename);

    arma::vec pos(3);
    arma::vec eField(3);
    PenningTrap trap = PenningTrap(B0, V0, d, parts);
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

void testOneParticleFE(){
    //trap parameters
    double B0 = 9.51e1;
    double V0 = 9.65e8;
    double d = 1e4;
    //evolution parameter
    double dt = 1000; //microseconds
    int N = 100*dt; //number of timesteps
    //init conditions
    double x0 = d/2;
    double z0 = d/2;
    double vy0 = d/dt;
    arma::vec initPos(3);
    arma::vec initVel(3);
    initPos(0) = x0; initPos(1) = 0; initPos(2) = z0;
    initVel(0) = 0; initVel(1) = vy0; initVel(2) = 0; 
    //particle parameters
    double mass = 40.078;
    double charge = 1; 
    Particle calsium = Particle(charge, mass, initPos, initVel);
    std::vector<Particle> particles;
    PenningTrap Trap = PenningTrap(B0, V0, d, particles);
    Trap.addParticle(calsium);

    for (int t = 0; t < N; t++){
        std::cout << scientificFormat(t);
        Trap.particles_.at(0).printCurrentPos();
        Trap.particles_.at(0).printCurrentVel();
        std::cout << "" << std::endl;
        Trap.evolveForwardEuler(dt);
        t += dt;
    }
}