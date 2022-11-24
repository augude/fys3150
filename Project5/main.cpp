#include "include/Solver.hpp"

int main()
{
    Solver obj(0.005, 2.5e-5, 0.008);
    obj.fill_matrices();
    obj.set_initial_state(0.25, 0.05, 200, 0.5, 0.05, 0);
    obj.set_potential("potential.txt", 0);
    obj.solve();
    obj.write_to_file(); 

    return 0;
}