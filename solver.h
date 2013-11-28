#include <iostream>
#include <cmath>
#include "planet.h"
#include <armadillo>
#include <time.h>
#include <fstream>

using namespace std;
using namespace arma;

class Solver {
public:
    int n;
    double t_min;
    double t_max;
    int N;
    double R0;
    double mu;
    double eps2;

    double h;
    double h2;
    double G;
    Planet **B;
    vec A;
    vec K;
    vec U;

    // constructor
    Solver(int, double, double, int, double, double, double);

    // initialization functions
    void init_cluster();
    void system_state_vector();

    // solver functions
    void RK4();
    void Leap_Frog();
    vec f(vec, int check = 1);

    // particle energies
    void energy();
};
