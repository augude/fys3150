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

