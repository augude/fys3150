#include "include/Solver.hpp"

int main()
{
    Solver obj(5, 5, 10);
    obj.fill_matrices();

    //Testing forward
    /* arma::cx_mat u = obj.set_initial_state(0.1, 0.1, 0.1, 0.2, 0.2, 0.2);
    u.reshape((u.n_rows*u.n_rows), 1);
    std::cout << obj.forward(u) << std::endl; */

    obj.set_initial_state(0.1, 0.1, 0.1, 0.2, 0.2, 0.2);
    obj.forward();

    return 0;
}