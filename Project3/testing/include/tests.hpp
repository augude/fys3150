#pragma once

//testing that the electric field looks as expected
void testPenningSetup();

//testing the time evolution of one particle with Forward Euler
void testOneParticleFE();

//testing the time evolution of one particle with RK4
void testOneParticleRK4();

//testing the movement of two particles in the trap
void testDoubleSetup(bool internalForces = true);

//calculate time evolution for different stepsizes using both FE and RK4
void compareStepsize(double stepSize);

//Calculates fraction of particles that are still trapped after 500mus as a function of angular frequency
void comparefValue(double f,double tStepSize,double wStepSize=0.02,double wStart=0.2,double wEnd=2.5, int nrParticles=100,bool coulombInteractions=true);
