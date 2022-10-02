#pragma once

//testing that the electric field looks as expected
void testPenningSetup();

//testing the time evolution of one particle with Forward Euler
void testOneParticleFE();

//testing the time evolution of one particle with RK4
void testOneParticleRK4();

//testing the movement of two particles in the trap
//NB: The code does not yet support the option of turning off the internal fields
void testDoubleSetup(bool internalForces = true);