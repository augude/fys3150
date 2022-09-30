#include "../include/PenningTrap.hpp"

PenningTrap::PenningTrap(double B0In, double V0In, double dIn, std::vector<Particle> particlesIn){
    B0_ = B0In;
    V0_ = V0In;
    d_ = dIn;
    particles_ = particlesIn;
}

void PenningTrap::addParticle(Particle pIn){

    particles_.push_back(pIn);
}

arma::vec PenningTrap::electricField(arma::vec & position){
    arma::vec E = arma::vec(3);
    E(0) = position(0);
    E(1) = position(1);
    E(2) = -2*position(2);
    E = E*V0_/(2*d_*d_);
    return E;
}

arma::vec PenningTrap::magneticField(arma::vec & postion){
    arma::vec B = arma::vec(3);
    B(0) = 0;
    B(1) = 0;
    B(2) = B0_;
    return B;
}

arma::vec PenningTrap::forceParticle(int i, int j){

    //throw an assertion if you try to calculate the force form the particle itself 
    assert(i != j);

    Particle particleOn = particles_.at(i);
    Particle particleFrom = particles_.at(j);

    double ke = 1.38935333e5;
    arma::vec distance = particleOn.postion_ - particleFrom.postion_;
    arma::vec force = ke*particleFrom.charge_*particleOn.charge_*distance/pow(norm(distance), 3);

    return force;
}

arma::vec PenningTrap::totalForceExternal(int i){
    Particle particleOn = particles_.at(i);
    
    arma::vec EField = electricField(particleOn.postion_);
    arma::vec BField = magneticField(particleOn.postion_);
    arma::vec externalForce = particleOn.charge_*(EField + cross(particleOn.velocity_, BField));

    return externalForce;
}

arma::vec PenningTrap::totalForceParticles(int i){
    arma::vec internalForce(3);
    int n = particles_.size();
    for (int j = 0; j < n; j++){
        if (j != i){
            internalForce += forceParticle(i, j);
        }
    }

    return internalForce;
}

arma::vec PenningTrap::totalForce(int i){
    arma::vec totalForce = totalForceExternal(i) + totalForceParticles(i);

    return totalForce;
}

void PenningTrap::evolveForwardEuler(double dt){
    int n = size(particles_);

    for (int i = 0; i < n; i++){
        Particle particleOn = particles_.at(i);
        arma::vec force = totalForce(i)/particleOn.mass_;
        //need to calculate change in postion before change in velocity to 'aviod' Euler-Cromer
        particleOn.postion_ += dt*particleOn.velocity_;
        particleOn.velocity_ += dt*force;
    }

}