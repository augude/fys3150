#include "../include/mcmc.hpp"
using namespace std;
using namespace arma;

void mcCycle(Lattice& s, vector<int> randomi, vector<int> randomj,
             vector<double> r) {
    int L = s.L;
    for (int i = 0; i < L * L; i++) {
        int ed =
            -2 *
            s.energyij(randomi[i],
                       randomj[i]);  // energy difference between the two states
        double A = s.energyDiff[ed];  // acceptance probability from the energy
                                      // difference using prerecorded values

        if (r[i] <= A) {
            // update the state, energy and magnetization
            int md = -2 * s.spins(randomi[i],
                                  randomj[i]);      // magnatization difference
                                                    // between the two states
            s.spins(randomi[i], randomj[i]) *= -1;  // update the spin ensemble
            s.E += ed;                              // update the energy
            s.M += md;                              // update the magnetization
        }
    }
}

void mcmc(double T, int L, int numberCycles, bool ordered, string filename) {
    mt19937 generator;
    // generator of random doubles on [0.0, 1.0)
    uniform_real_distribution<double> uniformDouble =
        uniform_real_distribution<double>(0.0, 1.0);
    // generator of random ints on [0, L - 1]
    uniform_int_distribution<> uniformInt(0, L - 1);

    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed);

    vector<int> randomi(L * L);
    vector<int> randomj(L * L);
    vector<double> r(L * L);

    // functions for filling vectors with random numbers
    auto genInt = [&uniformInt, &generator]() { return uniformInt(generator); };

    auto genDouble = [&uniformDouble, &generator]() {
        return uniformDouble(generator);
    };

    // initialize a lattice of size LxL at temperature T
    Lattice s = Lattice(L, T, ordered);

    mat EM = mat(numberCycles + 1, 2);
    EM(0, 0) = L;
    EM(0, 1) = T;

    s.energyMagnetization();  // initialize E and M

    for (int i = 0; i < numberCycles; i++) {
        generate(begin(randomi), end(randomi), genInt);
        generate(begin(randomj), end(randomj), genInt);
        generate(begin(r), end(r), genDouble);

        mcCycle(s, randomi, randomj, r);  // do one whole mc cycle and update
                                          // the state, energy and magnetization
        // write values to file
        EM(i + 1, 0) = s.E;
        EM(i + 1, 1) = s.M;
    }
    EM.save(filename);  // saving as binaryfile
}

