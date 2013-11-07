#include <iostream>
#include <armadillo>
#include <cmath>
#include "constants.h"
#include "planet.h"
#include <fstream>
#include <time.h>

using namespace std;
using namespace arma;

// function declarations
vec    f(vec A, Planet **B, int check = 1);
double kin_energy(vec A, Planet **B);
double pot_energy(vec A, Planet **B);
double ang_mom(vec A, Planet **B);
double lin_mom(Planet **B);

// number of objects
int n = 2;

int main(int argc, char* argv[])
{
    // initialization
    double t_min = 0.0;
    double t_max = 100;
    int    N = atoi(argv[1]);
    double h = (t_max - t_min)/N;
    double h2 = h/2.0;

    // make pointers to planet objects
    // data from: http://nssdc.gsfc.nasa.gov/planetary/factsheet/
    Planet *Sun     = new Planet(0.0,   0.0, 0.0, 0.0, 0.0,      0.0, MSun);
    Planet *Mercury = new Planet(0.387, 0.0, 0.0, 0.0, 3.215*pi, 0.0, 2.4e23);
    Planet *Venus   = new Planet(0.723, 0.0, 0.0, 0.0, 2.345*pi, 0.0, 4.87e24);
    Planet *Earth   = new Planet(1.0,   0.0, 0.0, 0.0, 2.0*pi,   0.0, 5.97e24);
    Planet *Mars    = new Planet(1.523, 0.0, 0.0, 0.0, 1.617*pi, 0.0, 6.6e23);
    Planet *Jupiter = new Planet(5.205, 0.0, 0.0, 0.0, 0.879*pi, 0.0, 1.898e27);
    Planet *Saturn  = new Planet(9.54,  0.0, 0.0, 0.0, 0.651*pi, 0.0, 5.5e26);
    Planet *Uranus  = new Planet(19.19, 0.0, 0.0, 0.0, 0.456*pi, 0.0, 8.8e25);
    Planet *Neptune = new Planet(30.06, 0.0, 0.0, 0.0, 0.362*pi, 0.0, 1.03e26);
    Planet *Pluto   = new Planet(39.53, 0.0, 0.0, 0.0, 0.315*pi, 0.0, 1.31e22);

    // system state vector
    vec A(6*n);

    // system vector: array of pointers to planet objects
    Planet **B = new Planet*[n];

    switch(n) {
    // two-body system
    case 2:
        B[0] = Sun;
        B[1] = Earth;
        break;

    // three-body system
    case 3:
        B[0] = Sun;
        B[1] = Earth;
        B[2] = Jupiter;
        break;

    // solar system
    case 10:
        B[0] = Sun;
        B[1] = Mercury;
        B[2] = Venus;
        B[3] = Earth;
        B[4] = Mars;
        B[5] = Jupiter;
        B[6] = Saturn;
        B[7] = Uranus;
        B[8] = Neptune;
        B[9] = Pluto;
        break;
    }

    // fixing center of mass at origin
    double p_tot = lin_mom(B);
    Sun->vy0 = -p_tot/Sun->M;

    // filling system state vector with initial values
    for (int i=0; i<n; i++) {
        A(i*6)   = B[i]->X0;
        A(i*6+1) = B[i]->Y0;
        A(i*6+2) = B[i]->Z0;
        A(i*6+3) = B[i]->vx0;
        A(i*6+4) = B[i]->vy0;
        A(i*6+5) = B[i]->vz0;
    }

    // calculate initial energy and angular momentum
    double K, U, L;
    K = kin_energy(A, B);
    U = pot_energy(A, B);
    L = ang_mom(A, B);

    // write inital positions of all objects and Earth's energymom to file
    fstream outFile1, outFile2, outFile3;
    outFile1.open("positions.dat", ios::out);
    outFile2.open("energymom.dat", ios::out);
    outFile3.open("radius.dat", ios::out);

    outFile1 << n << endl;
    for (int i=0; i<n; i++) outFile1 << A(6*i) << " " << A(6*i+1) << " " << A(6*i+2) << " ";
    outFile1 << endl;

    outFile2 << t_max << endl;
    outFile2 << K << " " << U << " " << L << endl;

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

        // calculate energy and angular momentum
        K = kin_energy(A, B);
        U = pot_energy(A, B);
        L = ang_mom(A, B);

        // write position values of all objects to file
        for (int j=0; j<n; j++) outFile1 << A(6*j) << " " << A(6*j+1) << " " << A(6*j+2) << " ";
        outFile1 << endl;

        // write energymom of Earth to file
        outFile2 << K << " " << U << " " << L << endl;

        // write out last radius of Earth
        if (i == N-1) {
            double r = sqrt(pow(A(6),2) + pow(A(7),2) + pow(A(8),2));
            outFile3 << t_max << endl << r << endl;
        }
    }
    outFile1.close();
    outFile2.close();
    outFile3.close();

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


// function to find the total linear momentum of the planets
double lin_mom(Planet **B)
{
    double p_tot = 0.0;
    for (int i=1; i<n; i++) p_tot += B[i]->M*sqrt(pow(B[i]->vx0,2) + pow(B[i]->vy0,2) + pow(B[i]->vz0,2));
    return p_tot;
}






