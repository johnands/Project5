#include "planet.h"

// constructor
Planet::Planet(double dX0, double dY0, double dZ0, double dvx0, double dvy0, double dvz0, double dM) {
    X0 = dX0;
    Y0 = dY0;
    Z0 = dZ0;
    vx0 = dvx0;
    vy0 = dvy0;
    vz0 = dvz0;

    M = dM/MSun;
}
