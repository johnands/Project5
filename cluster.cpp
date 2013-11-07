#include <iostream>
#include <armadillo>
#include <cmath>
#include "constants.h"
#include "planet.h"
#include <fstream>
#include <time.h>
#include "gaussiandeviate.cpp"

using namespace std;
using namespace arma;

// function declarations
vec    f(vec A, Planet **B, int check = 1);
double kin_energy(vec A, Planet **B);
double pot_energy(vec A, Planet **B);
double ang_mom(vec A, Planet **B);

// number of objects
int n = 100;

int main(int argc, char* argv[])
{
    // initialization
    double t_min = 0.0;
    double t_max = 50;
    int    N = atoi(argv[1]);
    double h = (t_max - t_min)/N;
    double h2 = h/2.0;
    double R0 = 20;                 // radius of sphere [l.y.]

    // system state vector
    vec A(6*n);

    // system vector: array of pointers to planet objects
    Planet **B = new Planet*[n];

    // initialize cluster
    long idum = -1;                 // seed
    double u, v, w;                 // uniformly dist. coordinates
    double r, theta, phi;           // spherical coordinates
    double x, y, z;                 // cartesian coordinates
    double sigma = 1, mu = 10;      // standard deviation and average
    double M;                       // mass in solar mass units
    for (int i=0; i<n; i++) {
        // fetch uniformly distributed random numbers between 0 and 1
        u = ran2(&idum); v = ran2(&idum); w = ran2(&idum);

        // convert to cartesian coordinates
        r = R0*pow(u,1.0/3); theta = acos(1-2*v); phi = 2*pi*w;
        x = r*sin(theta)*cos(phi);
        y = r*sin(theta)*sin(phi);
        z = r*cos(theta);

        // gaussian distributed masses
        M = sigma*gaussian_deviate(&idum) + mu;

        // make pointer to class object (zero initial velocity)
        B[i] = new Planet(x, y, z, 0.0, 0.0, 0.0, M);
    }

    // filling system state vector with initial values
    for (int i=0; i<n; i++) {
        A(i*6)   = B[i]->X0;
        A(i*6+1) = B[i]->Y0;
        A(i*6+2) = B[i]->Z0;
        A(i*6+3) = B[i]->vx0;
        A(i*6+4) = B[i]->vy0;
        A(i*6+5) = B[i]->vz0;
    }
    /*
    // calculate initial energy and angular momentum
    double K, U, L;
    K = kin_energy(A, B);
    U = pot_energy(A, B);
    L = ang_mom(A, B);
    */
    // write inital positions of all objects and Earth's energymom to file
    fstream outFile1;
    outFile1.open("positions.dat", ios::out);
    //outFile2.open("energymom.dat", ios::out);

    for (int i=0; i<n; i++) outFile1 << A(6*i) << " " << A(6*i+1) << " " << A(6*i+2) << " ";
    outFile1 << endl;

    //outFile2 << t_max << endl;
    //outFile2 << K << " " << U << " " << L << endl;

    // deciding which method to be used
    int method = atoi(argv[2]);

    // vectors to use in the algorithms
    vec k1(6*n), k2(6*n), k3(6*n), k4(6*n);
    vec v1(6*n), v2(6*n);

    // time integration
    for (int i=0; i<N; i++) {

        // Runge-Kutta 4
        if (method == 1) {
            // compute time usage
            clock_t start, finish;
            start = clock();

            k1 = f(A, B);
            k2 = f(A + k1*h2, B);
            k3 = f(A + k2*h2, B);
            k4 = f(A + k3*h, B);
            A += (h/6.0)*(k1 + 2*k2 + 2*k3 + k4);

            finish = clock();

            // write out time usage
            if (i == 0.5*N) cout << "Time elapsed RK4: " << ((finish-start)/CLOCKS_PER_SEC) << endl;
        }

        // Leap-Frog
        if (method == 2) {
            // compute time usage
            clock_t start, finish;
            start = clock();

            k1 = f(A, B, 3);                // ai, zeros elsewhere
            v1 = A + k1*h2;                 // v1 in final update equation
            k2 = f(A + k1*h2, B, 2);        // vi+1/2, zeros elsewhere; one of two vecs to find v2
            k3 = f(A + k2*h, B, 3);         // ai+1, zeros elsewhere; one of two vecs to find v2
            v2 = k2*h + k3*h2;              // v2 in final update equation
            A  = v1 + v2;                   // compute new A

            finish = clock();

            // write out time usage
            if (i == 0.5*N) cout << "Time elapsed Leap-Frog: " << ((finish-start)/CLOCKS_PER_SEC) << endl;
        }
        /*
        // calculate energy and angular momentum
        K = kin_energy(A, B);
        U = pot_energy(A, B);
        L = ang_mom(A, B);
        */
        // write position values of all objects to file
        for (int j=0; j<n; j++) outFile1 << A(6*j) << " " << A(6*j+1) << " " << A(6*j+2) << " ";
        outFile1 << endl;

        // write energymom of Earth to file
        //outFile2 << K << " " << U << " " << L << endl;
    }
    outFile1.close();
    //outFile2.close();

    return 0;
}


// function that computes derivate of vector A
vec f(vec A, Planet **B, int check)
{
    vec dAdt(6*n);  // derivative vector
    // setting dxdt and dydt for each planet
    // which in each case is just the velocity of that planet
    for (int i=0; i<n; i++) {
        // if check is equal to 3, only accelerations are computed
        if (check == 3) {
            dAdt(6*i)   = 0.0;      // dx/dt = vx
            dAdt(6*i+1) = 0.0;      // dy/dt = vy
            dAdt(6*i+2) = 0.0;      // dz/dt = vz
        }
        else {
            dAdt(6*i)   = A(6*i+3); // dx/dt = vx
            dAdt(6*i+1) = A(6*i+4); // dy/dt = vy
            dAdt(6*i+2) = A(6*i+5); // dz/dt = vz
        }

        // if check is equal to 2, only velocities are computed
        if (check == 2) {
            dAdt(6*i+3) = 0.0;
            dAdt(6*i+4) = 0.0;
            dAdt(6*i+5) = 0.0;
        }
        else {
            // finding the accelerations
            double ax, ay, az, r;
            ax = ay = az = 0.0;
            // i is current planet, j is the other planets
            // must find the force on planet i from the other planets j
            for (int j=0; j<n; j++) {
                if (i != j) {
                    r = sqrt(pow(A(6*j)-A(6*i),2) + pow(A(6*j+1)-A(6*i+1),2) + pow(A(6*j+2)-A(6*i+2),2));
                    ax += -(G*B[j]->M*(A(6*i)-A(6*j)))/(r*r*r);
                    ay += -(G*B[j]->M*(A(6*i+1)-A(6*j+1)))/(r*r*r);
                    az += -(G*B[j]->M*(A(6*i+2)-A(6*j+2))/(r*r*r));
                }
            }
            dAdt(6*i+3) = ax;
            dAdt(6*i+4) = ay;
            dAdt(6*i+5) = az;
        }
    }

return dAdt;
}


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


// function to find Earth's angular momentum in z-direction
double ang_mom(vec A, Planet **B)
{
    double Lx, Ly, Lz, L;
    if (n == 2 || n == 3) {
        Lx = A(7)*A(11) - A(8)*A(10);
        Ly = A(6)*A(11) - A(8)*A(9);
        Lz = A(6)*A(10) - A(7)*A(9);
    }
    L = B[1]->M*sqrt(pow(Lx,2) + pow(Ly,2) + pow(Lz,2));

    return L;

}
