#include "../include/PenningTrap.hpp"

PenningTrap::PenningTrap(double B0In, double V0In, double dIn, std::vector<Particle> particlesIn, double fIn, double wIn){
    B0 = B0In;
    V0 = V0In;
    d = dIn;
    f=fIn;
    w=wIn;

    particles = particlesIn;
    //TODO add changes to t into the functins
    double t = 0;
}

void PenningTrap::addParticle(Particle pIn){

    particles.push_back(pIn);
}

void PenningTrap::updateV0(){
    V0 = V0*(1+f*cos(w*t));
}

arma::vec PenningTrap::electricField(arma::vec position){
    arma::vec E = arma::vec(3);

    if (norm(position) > d){
        E.zeros();
    } else {
        E(0) = position(0);
        E(1) = position(1);
        E(2) = -2*position(2);
        E = E*V0/(d*d);
    }

    return E;
}

arma::vec PenningTrap::magneticField(arma::vec position){
    arma::vec B = arma::vec(3);

    if (norm(position) > d){
        B.zeros();
    } else {
        B(0) = 0;
        B(1) = 0;
        B(2) = B0;
    }

    return B;
}

arma::vec PenningTrap::forceParticle(int i, int j){

    //throw an assertion if you try to calculate the force from the particle itself 
    assert(i != j);

    Particle particleOn = particles.at(i);
    Particle particleFrom = particles.at(j);

    double ke = 1.38935333e5;
    arma::vec distance = particleOn.position - particleFrom.position;
    arma::vec force = ke*particleFrom.charge*particleOn.charge*distance/pow(norm(distance), 3);

    return force;
}

arma::vec PenningTrap::totalForceExternal(int i){
    Particle particleOn = particles.at(i);
    
    arma::vec EField = electricField(particleOn.position);
    arma::vec BField = magneticField(particleOn.position);
    arma::vec externalForce = particleOn.charge*(EField + cross(particleOn.velocity, BField));

    return externalForce;
}

arma::vec PenningTrap::totalForceParticles(int i){
    arma::vec internalForce(3);
    int n = particles.size();
    for (int j = 0; j < n; j++){
        if (j != i){
            internalForce += forceParticle(i, j);
        }
    }

    return internalForce;
}

arma::vec PenningTrap::totalForce(int i,bool coulombInteractions){
    arma::vec totalForce;
    if(coulombInteractions){
        totalForce = totalForceExternal(i) + totalForceParticles(i);
    } else{
        totalForce = totalForceExternal(i);
    }
    
    //std::cout << totalForceExternal(i) << std::endl;
    //std::cout << totalForceParticles(i) << std::endl;

    return totalForce;
}

int PenningTrap::countParticlesInside(){
    int nrInside=0;
    for(Particle p : particles){
        if(norm(p.position)<d){
            nrInside+=1;
        }
    }
    return nrInside;
}

void PenningTrap::evolveForwardEuler(double dt){
    int n = size(particles);

    //Creating new matrices to store new values temporarily
    arma::mat newPos(3, n);
    arma::mat newVel(3, n);

    for (int i = 0; i < n; i++){
        arma::vec force = totalForce(i)/particles.at(i).mass;
        //need to calculate change in position before change in velocity to 'aviod' Euler-Cromer

        newPos.col(i)=particles.at(i).position + dt*particles.at(i).velocity;
        newVel.col(i)=particles.at(i).velocity + dt*force;
    }

    for (int i = 0; i < n; i++){
        particles.at(i).position=newPos.col(i);
        particles.at(i).velocity=newVel.col(i);
    }

}

void PenningTrap::evolveRK4(double dt){
    int n = size(particles);

    arma::mat newPos(3, n);
    arma::mat newVel(3, n);

    for (int i = 0; i < n; i++){
        arma::vec posI = particles.at(i).position; //store pos at time i 
        arma::vec velI = particles.at(i).velocity; //store vel at time i

        arma::vec k1vel = dt*totalForce(i)/particles.at(i).mass;
        arma::vec k1pos = dt*particles.at(i).velocity;
        newVel.col(i) = velI + 0.5*k1vel; //calculate midpoint using k1
        newPos.col(i) = posI + 0.5*k1pos; //calculate midpoint using k1
        
        arma::vec k2vel = dt*totalForce(i)/particles.at(i).mass;
        arma::vec k2pos = dt*particles.at(i).velocity;
        newVel.col(i) = velI + 0.5*k2vel; //calculate midpoint using k2
        newPos.col(i) = posI + 0.5*k2pos; //calculate midpoint using k2
        
        arma::vec k3vel = dt*totalForce(i)/particles.at(i).mass;
        arma::vec k3pos = dt*particles.at(i).velocity;
        newVel.col(i) = velI + k3vel; //calculate endpoint using k3
        newPos.col(i) = posI + k3pos; //calculate endpoint using k3

        arma::vec k4vel = dt*totalForce(i)/particles.at(i).mass;
        arma::vec k4pos = dt*particles.at(i).velocity;
        
        newVel.col(i) = velI + 1.0/6*(k1vel + 2*k2vel + 2*k3vel + k4vel); //calculate endpoint using weighted sum
        newPos.col(i) = posI + 1.0/6*(k1pos + 2*k2pos + 2*k3pos + k4pos); //calculate endpoint using weighted sum
    }

    for (int i = 0; i < n; i++){
        particles.at(i).velocity = newVel.col(i);
        particles.at(i).position = newPos.col(i);
    }
}

void PenningTrap::evolveRK4witht(double dt,bool coulombInteractions){
    int n = size(particles);

    arma::mat newPos(3, n);
    arma::mat newVel(3, n);
    //std::cout << V0 << std::endl;
    //std::cout << t << std::endl

    for (int i = 0; i < n; i++){
        arma::vec posI = particles.at(i).position; //store pos at time i 
        arma::vec velI = particles.at(i).velocity; //store vel at time i

        arma::vec k1vel = dt*totalForce(i,coulombInteractions)/particles.at(i).mass;
        arma::vec k1pos = dt*particles.at(i).velocity;
        newVel.col(i) = velI + 0.5*k1vel; //calculate midpoint using k1
        newPos.col(i) = posI + 0.5*k1pos; //calculate midpoint using k1
        
        arma::vec k2vel = dt*totalForce(i,coulombInteractions)/particles.at(i).mass;
        arma::vec k2pos = dt*particles.at(i).velocity;
        newVel.col(i) = velI + 0.5*k2vel; //calculate midpoint using k2
        newPos.col(i) = posI + 0.5*k2pos; //calculate midpoint using k2
        
        arma::vec k3vel = dt*totalForce(i,coulombInteractions)/particles.at(i).mass;
        arma::vec k3pos = dt*particles.at(i).velocity;
        newVel.col(i) = velI + k3vel; //calculate endpoint using k3
        newPos.col(i) = posI + k3pos; //calculate endpoint using k3

        arma::vec k4vel = dt*totalForce(i,coulombInteractions)/particles.at(i).mass;
        arma::vec k4pos = dt*particles.at(i).velocity;
        
        newVel.col(i) = velI + 1.0/6*(k1vel + 2*k2vel + 2*k3vel + k4vel); //calculate endpoint using weighted sum
        newPos.col(i) = posI + 1.0/6*(k1pos + 2*k2pos + 2*k3pos + k4pos); //calculate endpoint using weighted sum
    }

    for (int i = 0; i < n; i++){
        particles.at(i).velocity = newVel.col(i);
        particles.at(i).position = newPos.col(i);
    }

    //Update t and then update V0 for next evolution
    t+=dt;
    updateV0();
}