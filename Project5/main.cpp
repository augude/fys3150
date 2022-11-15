#include "include/Solver.hpp"

int main()
{
    Solver obj(5, 5, 10);
    obj.set_initial_state(0.1, 0.1, 0.1, 0.2, 0.2, 0.2);

    return 0;
}