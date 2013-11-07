#include <iostream>
#include <cmath>
#include "constants.h"
#include "planet.h"
#include <armadillo>
#include <time.h>
#include <fstream>

using namespace std;
using namespace arma;

// class to make planet objects containing data of each planet: mass, x0, y0, vx0, vy0
class Solver {
public:
    int n;
    double t_min;
    double t_max;
    int N;
    double R0;
    double mu;

    double h;
    double h2;
    double G;
    Planet **B;
    vec A;

    // file variables
    fstream *outFile1;
    //fstream outFile2;

    // constructor
    Solver(int, double, double, int, double, double);

    // initialization functions
    void init_cluster();
    void system_state_vector();

    // solver functions
    void RK4();
    void Leap_Frog();
    vec f(vec, int check = 1);
    void write_to_file();
};
