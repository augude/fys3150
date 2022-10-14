#include "../include/PenningTrap.hpp"

PenningTrap::PenningTrap(double B0In, double V0In, double dIn, std::vector<Particle> particlesIn, double fIn, double omegaIn){
    B0 = B0In;
    V0 = V0In;
    d = dIn;
    particles = particlesIn;
    f = fIn;
    omega = omegaIn;
}

void PenningTrap::addParticle(Particle pIn){

    particles.push_back(pIn);
}

arma::vec PenningTrap::electricField(arma::vec position, double time){
    arma::vec E = arma::vec(3);
    // check if particle is in the trap
    if (norm(position) > d){
        E(0) = 0; E(1) = 0; E(2) = 0;
    }
    else{
        E(0) = position(0);
        E(1) = position(1);
        E(2) = -2*position(2);
        double V = V0*(1 + f*cos(omega*time));
        E = E*V/(d*d);
    }
    return E;
}

arma::vec PenningTrap::magneticField(arma::vec position, double time){
    arma::vec B = arma::vec(3);
    // check if particle is in the trap
    if (norm(position) > d){
        B(0) = 0; B(1) = 0; B(2) = 0;
    }
    else{
        B(0) = 0;
        B(1) = 0;
        B(2) = B0; 
    }
    return B;
}

arma::vec PenningTrap::forceParticle(int i, int j, double time){

    //thcol an assertion if you try to calculate the force from the particle itself 
    assert(i != j);

    Particle particleOn = particles.at(i);
    Particle particleFrom = particles.at(j);

    double ke = 1.38935333e5;
    arma::vec distance = particleOn.position - particleFrom.position;
    arma::vec force = ke*particleFrom.charge*particleOn.charge*distance/pow(norm(distance), 3);

    return force;
}

arma::vec PenningTrap::totalForceExternal(int i, double time){
    Particle particleOn = particles.at(i);
    
    arma::vec EField = electricField(particleOn.position, time);
    arma::vec BField = magneticField(particleOn.position);
    arma::vec externalForce = particleOn.charge*(EField + cross(particleOn.velocity, BField));

    return externalForce;
}

arma::vec PenningTrap::totalForceParticles(int i, double time){
    arma::vec internalForce(3);
    int n = particles.size();
    for (int j = 0; j < n; j++){
        if (j != i){
            internalForce += forceParticle(i, j, time);
        }
    }

    return internalForce;
}

arma::vec PenningTrap::totalForce(int i, double time, bool internal){
    arma::vec totalForce;
    if (internal){
        totalForce = totalForceExternal(i, time) + totalForceParticles(i, time);
    }
    else{
        totalForce = totalForceExternal(i, time);
    }
    

    return totalForce;
}

void PenningTrap::evolveForwardEuler(double dt, double time, bool internal){
    int n = size(particles);

    //Creating new matrices to store new values temporarily
    arma::mat newPos(3, n);
    arma::mat newVel(3, n);

    for (int i = 0; i < n; i++){
        arma::vec force = totalForce(i, time, internal)/particles.at(i).mass;
        //need to calculate change in position before change in velocity to 'aviod' Euler-Cromer

        newPos.col(i) = particles.at(i).position + dt*particles.at(i).velocity;
        newVel.col(i) = particles.at(i).velocity + dt*force;
    }

    for (int i = 0; i < n; i++){
        particles.at(i).position = newPos.col(i);
        particles.at(i).velocity = newVel.col(i);
    }

}

void PenningTrap::evolveRK4(double dt, double time, bool internal){
    int n = size(particles);

    arma::mat posI(3, n); //to store initial pos
    arma::mat velI(3, n); //to store initial vel

    arma::mat k1vel(3, n); //to store k's
    arma::mat k1pos(3, n);
    arma::mat k2vel(3, n);
    arma::mat k2pos(3, n);
    arma::mat k3vel(3, n);
    arma::mat k3pos(3, n);
    arma::mat k4vel(3, n);
    arma::mat k4pos(3, n);
    //update all particles using k1
    for (int i = 0; i < n; i++){
        posI.col(i) = particles.at(i).position; //store old initial positions
        velI.col(i) = particles.at(i).velocity; //store old initial velocities
        k1vel.col(i) = totalForce(i, time, internal)/particles.at(i).mass;
        k1pos.col(i) = particles.at(i).velocity;
        particles.at(i).velocity = velI.col(i) + 0.5*dt*k1vel.col(i); //calculate midpoint using k1
        particles.at(i).position = posI.col(i) + 0.5*dt*k1pos.col(i); //calculate midpoint using k1
        
    }
    //update all particles using k2
    for (int i = 0; i < n; i++){
        k2vel.col(i) = totalForce(i, time, internal)/particles.at(i).mass;
        k2pos.col(i) = particles.at(i).velocity;
        particles.at(i).velocity = velI.col(i) + 0.5*dt*k2vel.col(i); //calculate midpoint using k2
        particles.at(i).position = posI.col(i) + 0.5*dt*k2pos.col(i); //calculate midpoint using k2
    }
    //update all particles using k3
    for (int i = 0; i < n; i++){
        k3vel.col(i) = totalForce(i, time, internal)/particles.at(i).mass;
        k3pos.col(i) = particles.at(i).velocity;
        particles.at(i).velocity = velI.col(i) + dt*k3vel.col(i); //calculate endpoint using k3
        particles.at(i).position = posI.col(i) + dt*k3pos.col(i); //calculate endpoint using k3
    } 
    //update all particles using weigthed sum
    for (int i = 0; i < n; i++){
        k4vel.col(i) = totalForce(i, time, internal)/particles.at(i).mass;
        k4pos.col(i) = particles.at(i).velocity;
        particles.at(i).velocity = velI.col(i) + dt/6*(k1vel.col(i) + 2*k2vel.col(i) + 2*k3vel.col(i) + k4vel.col(i)); //calculate endpoint using weigthed sum
        particles.at(i).position = posI.col(i) + dt/6*(k1pos.col(i) + 2*k2pos.col(i) + 2*k3pos.col(i) + k4pos.col(i)); //calculate endpoint using weigthed sum
    }   
}

int PenningTrap::numberWithin(){
    int n = particles.size();
    int number = 0;
    for (int i = 0; i < n; i++){
        number += (norm(particles.at(i).position) < d);
    }
    return number;
}

void PenningTrap::fillTrap(int N){
    arma::arma_rng::set_seed_random();
    double mass = 40.078;
    double charge = 1; 
    for (int i = 0; i < N; i++){
        arma::vec pos = arma::vec(3).randn() * 0.1 * d;  // random initial position
        arma::vec vel = arma::vec(3).randn() * 0.1 * d;  // random initial velocity
        Particle ion = Particle(charge, mass, pos, vel);
        addParticle(ion);
    }
}