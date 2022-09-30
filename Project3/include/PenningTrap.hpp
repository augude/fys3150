#pragma once 
#include <iostream>
#include <armadillo>
#include <vector>
#include <cassert>
#include "Particle.hpp"

class PenningTrap{
    public:
        double B0_;
        double V0_;
        double d_;
        std::vector<Particle> particles_;

        PenningTrap(double B0In, double V0In, double dIn, std::vector<Particle> particlesIn);

        // Add a particle to the trap
        void addParticle(Particle pIn);

        // External electric field at point r=(x,y,z)
        arma::vec electricField(arma::vec position);

        // External magnetic field at point r=(x,y,z)
        arma::vec magneticField(arma::vec position);

        // Force on particle_i from particle_j
        arma::vec forceParticle(int i, int j);
        
        // The total force on particle_i from the external fields
        arma::vec totalForceExternal(int i);

        // The total force on particle_i from the other particles
        arma::vec totalForceParticles(int i);

        // The total force on particle_i from both external fields and other particles
        arma::vec totalForce(int i);

        // Evolve the system one time step (dt) using Forward Euler
        void evolveForwardEuler(double dt);
};