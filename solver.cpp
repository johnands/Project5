#include "solver.h"
#include "gaussiandeviate.cpp"

#define pi 3.141592653589793

// constructor
Solver::Solver(int dn, double dt_min, double dt_max, int dN, double dR0, double dmu, double eps) {
    n = dn;                             // number of objects
    t_min = dt_min;
    t_max = dt_max;
    N = dN;                             // number of time points
    R0 = dR0;                           // radius of sphere [ly]
    mu = dmu;                           // average mass [Solar mass]
    eps2 = eps*eps;                     // smoothing constant

    h = (t_max - t_min)/N;              // time step
    h2 = h/2.0;
    G = (pi*pi*pow(R0,3))/(8*mu*n);     // gravitational constant [ly^3/tau^2]
    B = new Planet*[n];
    A = zeros(6*n);                     // system state vector
    K = zeros(n);                       // kinetic energies
    U = zeros(n);                       // potential energies
}

void Solver::init_cluster() {
    // initialize cluster
    long idum = -1;                 // seed
    double u, v, w;                 // uniformly dist. coordinates
    double r, theta, phi;           // spherical coordinates
    double x, y, z;                 // cartesian coordinates
    double sigma = 1;               // standard deviation
    double M;                       // mass [Solar mass]
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

}

void Solver::system_state_vector() {
    // filling system state vector with initial values
    for (int i=0; i<n; i++) {
        A(i*6)   = B[i]->X0;
        A(i*6+1) = B[i]->Y0;
        A(i*6+2) = B[i]->Z0;
        A(i*6+3) = B[i]->vx0;
        A(i*6+4) = B[i]->vy0;
        A(i*6+5) = B[i]->vz0;
    }
}


void Solver::RK4() {

    // compute initial energies
    energy();

    // write initial positions and energies to file
    ofstream outFile1;           ofstream outFile3;
    outFile1.open("pos_RK4.dat"); outFile3.open("energies_RK4.dat");
    for (int j=0; j<n; j++) {
        outFile1 << A(6*j) << " " << A(6*j+1) << " " << A(6*j+2) << " ";
        outFile3 << K(j) << " " << U(j) << " ";
    }
    outFile1 << endl;
    outFile3 << endl;

    // vectors to use in the algorithm
    vec k1(6*n), k2(6*n), k3(6*n), k4(6*n);

    // compute time usage
    clock_t start, finish;
    start = clock();

    // time integration
    for (int i=0; i<N; i++) {

        k1 = f(A);
        k2 = f(A + k1*h2);
        k3 = f(A + k2*h2);
        k4 = f(A + k3*h);
        A += (h/6.0)*(k1 + 2*k2 + 2*k3 + k4);

        // calculate new energies
        energy();

        // write positions and energies to files
        for (int j=0; j<n; j++) {
            outFile1 << A(6*j) << " " << A(6*j+1) << " " << A(6*j+2) << " ";
            outFile3 << K(j) << " " << U(j) << " ";
        }
        outFile1 << endl;
        outFile3 << endl;
    }
    finish = clock();
    // write out time usage
    cout << "Time elapsed RK4: " << ((finish-start)/CLOCKS_PER_SEC) << " s" << endl;

    // close files
    outFile1.close();
    outFile3.close();
}

void Solver::Leap_Frog() {

    // compute initial energies
    energy();

    // write initial positions and energies to file
    ofstream outFile2;           ofstream outFile4;
    outFile2.open("pos_LF.dat"); outFile4.open("energies_LF.dat");
    for (int j=0; j<n; j++) {
        outFile2 << A(6*j) << " " << A(6*j+1) << " " << A(6*j+2) << " ";
        outFile4 << K(j) << " " << U(j) << " ";
    }
    outFile2 << endl;
    outFile4 << endl;

    // vectors to use in the algorithm
    vec k1(6*n), k2(6*n), k3(6*n);
    vec v1(6*n), v2(6*n);

    // compute time usage
    clock_t start, finish;
    start = clock();

    // time integration
    for (int i=0; i<N; i++) {

        k1 = f(A, 3);                // ai, zeros elsewhere
        v1 = A + k1*h2;              // v1 in final update equation
        k2 = f(A + k1*h2, 2);        // vi+1/2, zeros elsewhere; one of two vecs to find v2
        k3 = f(A + k2*h, 3);         // ai+1, zeros elsewhere; one of two vecs to find v2
        v2 = k2*h + k3*h2;           // v2 in final update equation
        A  = v1 + v2;                // compute new A

        // calculate new energies
        energy();

        // write positions and energies to files
        for (int j=0; j<n; j++) {
            outFile2 << A(6*j) << " " << A(6*j+1) << " " << A(6*j+2) << " ";
            outFile4 << K(j) << " " << U(j) << " ";
        }
        outFile2 << endl;
        outFile4 << endl;
    }
    finish = clock();
    // write out time usage
    cout << "Time elapsed Leap-Frog: " << ((finish-start)/CLOCKS_PER_SEC) << " s" << endl;

    // close files
    outFile2.close();
    outFile4.close();
}

vec Solver::f(vec A, int check) {
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
                    ax += -(G*B[j]->M*(A(6*i)-A(6*j)))/(r*(r*r + eps2));
                    ay += -(G*B[j]->M*(A(6*i+1)-A(6*j+1)))/(r*(r*r + eps2));
                    az += -(G*B[j]->M*(A(6*i+2)-A(6*j+2)))/(r*(r*r + eps2));
                }
            }
            dAdt(6*i+3) = ax;
            dAdt(6*i+4) = ay;
            dAdt(6*i+5) = az;
        }
    }

return dAdt;
}


// function to find particle energies
void Solver::energy()
{
    double r, dU = 0;
    for (int i=0; i<n; i++) {
        K(i) = 0.5*B[i]->M*(pow(A(6*i+3),2) + pow(A(6*i+4),2) + pow(A(6*i+5),2));
        for (int j=0; j<n; j++) {
            if (i != j) {
                r = sqrt(pow(A(6*j)-A(6*i),2) + pow(A(6*j+1)-A(6*i+1),2) + pow(A(6*j+2)-A(6*i+2),2));
                dU += -(G*B[i]->M*B[j]->M)/r;
            }
        }
        U(i) = dU;

    }
}

