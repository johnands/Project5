#include <iostream>
#include <armadillo>
#include <cmath>
#include "solver.h"
#include <time.h>
#include <fstream>

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
    ofstream outFile5;
    outFile5.open("ejection.dat");

    // n values
    vec n_part = linspace<vec>(50, 150, 10);

    // find ejection as function of n
    //for (int i=0; i<n_part.n_elem; i++) {

        //                        n       t_min  t_max    N      R0    mu     eps
        Solver object = Solver(100,   0,    10,     5000,   20,   10,    0.15);
        object.init_cluster();
        object.system_state_vector();

        // decide method (1 = RK4, 2 = Leap-Frog)
        int method = atoi(argv[1]);

        if (method == 1) object.RK4();
        if (method == 2) object.Leap_Frog();

        double Etot = sum(object.K) + sum(object.U);
        outFile5 << object.ejected_particles << endl;
        outFile5 << object.ejected_energy/Etot << endl;
    //}

    return 0;
}
