#pragma once 
#include <iostream>
#include <armadillo>
#include <vector>
#include <cassert>
#include <math.h>
#include <armadillo>
#include "Particle.hpp"

class PenningTrap{
    public:
        double B0;
        double V0;
        double d;
        double f;
        double omega;
        std::vector<Particle> particles;

        PenningTrap(double B0In, double V0In, double dIn, std::vector<Particle> particlesIn, double f = 0, double omega = 0);

        // Add a particle to the trap
        void addParticle(Particle pIn);

        // External electric field at point r=(x,y,z) and time=t
        arma::vec electricField(arma::vec position, double time = 0);

        // External magnetic field at point r=(x,y,z) and time = t
        arma::vec magneticField(arma::vec position, double time = 0);

        // Force on particle i from particle j
        arma::vec forceParticle(int i, int j, double time = 0);
        
        // The total force on particle i from the external fields
        arma::vec totalForceExternal(int i, double time = 0);

        // The total force on particle i from the other particles
        arma::vec totalForceParticles(int i, double time = 0);

        // The total force on particle i from both external fields and other particles
        arma::vec totalForce(int i, double time = 0, bool internal = true);

        // Evolve the system one time step (dt) using Forward Euler
        void evolveForwardEuler(double dt, double time = 0, bool internal = true);

        // Evolve the system one time step (dt) using RK4
        void evolveRK4(double dt, double time = 0, bool internal = true);

        //count the number of particles within the trap
        int numberWithin();

        //fill the trap with N Ca+ ions
        void fillTrap(int N);
};