#include <iostream>
#include <armadillo>
#include <cmath>
#include "solver.h"
#include <time.h>

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
    //                      n   t_min  t_max    N      R0    mu     eps
    Solver object = Solver(100,   0,    5,     5000,   20,   10,    0.15);
    object.init_cluster();
    object.system_state_vector();

    // decide method (1 = RK4, 2 = Leap-Frog)
    int method = atoi(argv[1]);

    if (method == 1) object.RK4();
    if (method == 2) object.Leap_Frog();

    return 0;
}
