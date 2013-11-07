#include <iostream>
#include <armadillo>
#include <cmath>
#include "constants.h"
#include "solver.h"
#include <fstream>
#include <time.h>

using namespace std;
using namespace arma;

// function declarations
//double kin_energy(vec A, Planet **B);
//double pot_energy(vec A, Planet **B);

int main(int argc, char* argv[])
{

    //                      n   t_min  t_max     N     R0    mu
    Solver object = Solver(100,   0,    5,     500,   20,   10);
    object.init_cluster();
    object.system_state_vector();

    // decide method (1 = RK4, 2 = Leap-Frog)
    int method = atoi(argv[1]);

    if (method == 1) object.RK4();
    if (method == 2) object.Leap_Frog();


    /*
    // calculate initial energy and angular momentum
    double K, U, L;
    K = kin_energy(A, B);
    U = pot_energy(A, B);
    L = ang_mom(A, B);
    */

    //outFile2 << t_max << endl;
    //outFile2 << K << " " << U << " " << L << endl;

    return 0;
}

/*
// function to find Earth's kinetic energy
double kin_energy(vec A, Planet **B)
{
    double K, vx, vy, vz;
    if (n == 2 || n == 3) {
        vx = A(8); vy = A(9); vz = A(10);
    }
    else {
        vx = A(20); vy = A(21); vz = A(22);
    }
    K = 0.5*B[1]->M*(pow(vx,2) + pow(vy,2) + pow(vz,2));
    return K;
}


// function to find Earth's potential energy
double pot_energy(vec A, Planet **B)
{
    double r, U;
    if (n == 2 || n == 3) {
        r = sqrt(pow(A(6)-A(0),2) + pow(A(7)-A(1),2) + pow(A(8)-A(2),2));
        U = -(G*B[0]->M*B[1]->M)/r;
    }
    else {
        r = sqrt(pow(A(18)-A(0),2) + pow(A(19)-A(1),2) + pow(A(20)-A(2),2));
        U = -(G*B[0]->M*B[3]->M)/r;
    }
    return U;
}
*/
