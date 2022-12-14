#include "../../include/PenningTrap.hpp"
#include "../../../Project2/include/utils.hpp"
#include "../../include/Particle.hpp"
#include <fstream>
#include <string>
#include <time.h>

void testPenningSetup(){

    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500;

    arma::vec particlePos;
    arma::vec particleVec;
    double charge;
    double mass;

    Particle part = Particle(charge, mass, particlePos, particleVec);

    std::vector<Particle> parts = {part};
    std::string filename = "ElectricField.txt";
    arma::vec pos(3);
    arma::vec eField(3);
    PenningTrap trap = PenningTrap(B0, V0, d, parts);
    int N = 10;
    for (int z = 0; z < N; z++){
        for (int y = 0; y < N; y++){
            double zpos = 2*d*z/N - d;
            double ypos = 2*d*y/N - d;
            pos(0) = 0;
            pos(1) = ypos;
            pos(2) = zpos;
            double V = V0/(2*d*d)*(2*zpos*zpos - ypos*ypos - 0*0);
            eField = trap.electricField(pos);  
            std::cout << scientificFormat(z) << scientificFormat(y);
            std::cout << scientificFormat(pos); 
            std::cout << scientificFormat(eField) << scientificFormat(V) << std::endl;
        }
    }
}

void testOneParticleFE(){
    //trap parameters
    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500;
    //evolution parameter
    double dt = 0.001; //microseconds
    double T = 50.0; //end time
    int N = T/dt; //number of timesteps
    //init conditions
    double x0 = d/2;
    double z0 = d/2;
    double vy0 = 10;
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

    for (int i = 0; i < N; i++){
        double t = dt*i;
        std::cout << scientificFormat(t);
        Trap.particles.at(0).printCurrentPos();
        Trap.particles.at(0).printCurrentVel();
        std::cout << "" << std::endl;
        Trap.evolveForwardEuler(dt);
    }   
}

void testOneParticleRK4(){
    //trap parameters
    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500;
    //evolution parameters
    double dt = 0.001; //microseconds
    double T = 50.0; //end time
    int N = T/dt; //number of timesteps
    //init conditions
    double x0 = d/2;
    double z0 = d/2;
    double vy0 = 10;
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

    for (int i = 0; i < N; i++){
        double t = dt*i;
        std::cout << scientificFormat(t);
        Trap.particles.at(0).printCurrentPos();
        Trap.particles.at(0).printCurrentVel();
        std::cout << "" << std::endl;
        Trap.evolveRK4(dt);
    }   
}

void testDoubleSetup(bool internalForces = true){
    //trap parameters
    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500;
    //evolution parameters
    double dt = 0.001; //microseconds
    int T = 50; //end time
    int N = T/dt; //number of timesteps
    //init conditions, must have different init pos for the particles to avoid infinite Coulomb-forces
    arma::vec initPos1(3); arma::vec initVel1(3);
    arma::vec initPos2(3); arma::vec initVel2(3);
    initPos1(0) = 20.0; initPos1(1) = 0.0; initPos1(2) = 20.0;
    initVel1(0) = 0.0; initVel1(1) = 25.0; initVel1(2) = 0.0; 
    initPos2(0) = 25.0; initPos2(1) = 25.0; initPos2(2) = 0.0;
    initVel2(0) = 0.0; initVel2(1) = 40.0; initVel2(2) = 5.0; 
    //particle parameters
    double mass = 40.078;
    double charge = 1; 
    
    if (internalForces == true){
        //add particles to the same trap if you want internal forces
        Particle calsium1 = Particle(charge, mass, initPos1, initVel1);
        Particle calsium2 = Particle(charge, mass, initPos2, initVel2);
        std::vector<Particle> particles;
        PenningTrap Trap = PenningTrap(B0, V0, d, particles);
        Trap.addParticle(calsium1);
        Trap.addParticle(calsium2);

        for (int i = 0; i < N; i++){
            double t = dt*i;
            std::cout << scientificFormat(t);
            for (int n = 0; n < size(Trap.particles); n++){
                Trap.particles.at(n).printCurrentPos();
                Trap.particles.at(n).printCurrentVel();
            }
            std::cout << "" << std::endl;
            Trap.evolveRK4(dt);
        }   
    }
    else{
        //add particles to identical, but different traps if you don't want internal forces
        Particle calsium1 = Particle(charge, mass, initPos1, initVel1);
        Particle calsium2 = Particle(charge, mass, initPos2, initVel2);
        std::vector<Particle> particles;
        PenningTrap Trap1 = PenningTrap(B0, V0, d, particles);
        PenningTrap Trap2 = PenningTrap(B0, V0, d, particles);
        Trap1.addParticle(calsium1);
        Trap2.addParticle(calsium2);

        for (int i = 0; i < N; i++){
            double t = dt*i;
            std::cout << scientificFormat(t);            
            Trap1.particles.at(0).printCurrentPos();
            Trap1.particles.at(0).printCurrentVel();
            Trap2.particles.at(0).printCurrentPos();
            Trap2.particles.at(0).printCurrentVel();
            std::cout << "" << std::endl;
            Trap1.evolveRK4(dt);
            Trap2.evolveRK4(dt);
        }   
    } 
}

void compareStepsize(int N){
    //trap parameters
    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500;
    //evolution parameter
    double T = 50.0; //end time
    double stepSize = T/N; //number of timesteps
    //init conditions
    double x0 = d/2;
    double z0 = d/2;
    double vy0 = 10;
    arma::vec initPos(3);
    arma::vec initVel(3);
    initPos(0) = x0; initPos(1) = 0; initPos(2) = z0;
    initVel(0) = 0; initVel(1) = vy0; initVel(2) = 0; 
    //particle parameters
    double mass = 40.078;
    double charge = 1; 
    Particle calsium = Particle(charge, mass, initPos, initVel);
    std::vector<Particle> particles;
    PenningTrap TrapFE = PenningTrap(B0, V0, d, particles);
    PenningTrap TrapRK4 = PenningTrap(B0, V0, d, particles);
    TrapFE.addParticle(calsium);
    TrapRK4.addParticle(calsium);

    for (int i = 0; i < N; i++){
        double t = stepSize*i;
        std::cout << scientificFormat(t);
        TrapFE.particles.at(0).printCurrentPos();
        TrapFE.particles.at(0).printCurrentVel();
        TrapRK4.particles.at(0).printCurrentPos();
        TrapRK4.particles.at(0).printCurrentVel();
        std::cout << "" << std::endl;
        TrapFE.evolveForwardEuler(stepSize);
        TrapRK4.evolveRK4(stepSize);
    }   
}


arma::vec fractionWithin(double f){
    //trap parameters
    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500;
    int numberParticles = 100;
    //evolution parameter
    double T = 500.0; //end time
    double stepSize = 0.01; //number of timesteps
    int numberSteps = T/stepSize;
    std::vector<Particle> particles;
    arma::vec omega = arma::linspace(0.2, 2.5, 115); 
    arma::vec fraction = arma::vec(omega.size());
    int index = 0;
    for (double o : omega){
        PenningTrap Trap = PenningTrap(B0, V0, d, particles, f, o);
        Trap.fillTrap(numberParticles); //fill trap with 100 particles
        for (int i = 0; i < numberSteps; i++){
            double time = stepSize*i;
            Trap.evolveRK4(stepSize, time, false);
        }
        fraction(index) = Trap.numberWithin();
        index += 1;
    }
    return fraction;
}

arma::vec fractionWithinZoom(bool internal){
    //trap parameters
    double B0 = 9.65e1;
    double V0 = 2.41e6;
    double d = 500;
    double f = 0.7;
    int numberParticles = 100;
    //evolution parameter
    double T = 500.0; //end time
    double stepSize = 0.01; //number of timesteps
    int numberSteps = T/stepSize;
    std::vector<Particle> particles;
    arma::vec omega = arma::linspace(0.4, 0.5, 20); 
    arma::vec fraction = arma::vec(omega.size());
    int index = 0;
    for (double o : omega){
        PenningTrap Trap = PenningTrap(B0, V0, d, particles, f, o);
        Trap.fillTrap(numberParticles); //fill trap with 100 particles
        clock_t t1 = clock();
        for (int i = 0; i < numberSteps; i++){
            double time = stepSize*i;
            Trap.evolveRK4(stepSize, time, internal);
        }
        clock_t t2 = clock();
        std::cout << "Iteration " << index +1 << " of " << omega.size() << " used " << ((double) (t2 - t1)) / CLOCKS_PER_SEC << " seconds" << std::endl;
        fraction(index) = Trap.numberWithin();
        index += 1;
    }
    return fraction;

}