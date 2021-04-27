/*  Dies ist ein Programm zur Implementierung des forward-euler,
    rk4 und des leapfrog-Verfahrens.
*/

//pre-processor instructions
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>
#include <cfloat>
#include <chrono>


//for global "short-hand" notation - need not to write 'std::' in front of most things
using namespace std;

void initialize(int n, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz){
    //Change size of vectors
    x.resize(n);
    y.resize(n);
    z.resize(n);
    vx.resize(n);
    vy.resize(n);
    vz.resize(n);

    //set startvalues
    x[0] = 0.;
    x[1] = 0.5;
    y[0] = 0.;
    y[1] = 0.;
    z[0] = 0.;
    z[1] = 0.;

    //fwd and rk4 method need initialisation of v at the same timestep as x
    vx[0] = 0.;
    vx[1] = 0.;
    vy[0] = 0.;
    vy[1] = 0.;
    vz[0] = 0.;
    vz[1] = 0.;
    
    //mabye initialize natural constants
}

void initialize_symplectic(int n, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz){
    //Change size of vectors
    x.resize(n);
    y.resize(n);
    z.resize(n);
    vx.resize(n);
    vy.resize(n);
    vz.resize(n);

    //set startvalues
    x[0] = 0.;
    x[1] = 0.5;
    y[0] = 0.;
    y[1] = 0.;
    z[0] = 0.;
    z[1] = 0.;

    //contrary the lf method needs initialisation of v a half timestep before the first timestep of x
    vx[0] = 0.999999875;
    vx[1] = 0.999999875;
    vy[0] = 0.;
    vy[1] = 0.;
    vz[0] = 0.;
    vz[1] = 0.;
    
    //mabye initialize natural constants
}

void *testfunction(double t, vector<double> x, vector<double> &u_rhs, int n){
    int i;
    u_rhs.resize(n);
    for(int i=0; i<n; i++) u_rhs[i] = cos(t);
}

void *testsymplectic(double t, vector<double> x, vector<double> &u_rhs, int n){
    int i;
    u_rhs.resize(n);
    for(int i=0; i<n; i++) u_rhs[i] = -sin(t);
}

void fwd_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, void* rhs(double t, vector<double> x, vector<double> &u_rhs, int n), int n){
    //Initialize vectors for the steps - only one step here!
    vector<double> kx, ky, kz;

    for(int i=0; i<n; i++) kx.push_back(0.);
    for(int i=0; i<n; i++) ky.push_back(0.);
    for(int i=0; i<n; i++) kz.push_back(0.);

    //determine the derivative in every direction (calculate kx, ky, kz)
    rhs(t, x, kx, n);
    rhs(t, y, ky, n);
    rhs(t, z, kz, n);

    //do the iteration step (update the positions)
    for(int i=0; i<n; i++) x[i] += dt * kx[i];
    for(int i=0; i<n; i++) y[i] += dt * ky[i];
    for(int i=0; i<n; i++) z[i] += dt * kz[i];
}

void rk4_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, void* rhs(double t, vector<double> x, vector<double> &u_rhs, int n), int n){
    //Initialize vectors for the steps - only one step here!
    vector<double> kx1, kx2, kx3, kx4, tmpx;
    vector<double> ky1, ky2, ky3, ky4, tmpy;
    vector<double> kz1, kz2, kz3, kz4, tmpz;

    for(int i=0; i<n; i++) kx1.push_back(0.);
    for(int i=0; i<n; i++) kx2.push_back(0.);
    for(int i=0; i<n; i++) kx3.push_back(0.);
    for(int i=0; i<n; i++) kx4.push_back(0.);
    for(int i=0; i<n; i++) ky1.push_back(0.);
    for(int i=0; i<n; i++) ky2.push_back(0.);
    for(int i=0; i<n; i++) ky3.push_back(0.);
    for(int i=0; i<n; i++) ky4.push_back(0.);
    for(int i=0; i<n; i++) kz1.push_back(0.);
    for(int i=0; i<n; i++) kz2.push_back(0.);
    for(int i=0; i<n; i++) kz3.push_back(0.);
    for(int i=0; i<n; i++) kz4.push_back(0.);

    for(int i=0; i<n; i++) tmpx.push_back(0.);
    for(int i=0; i<n; i++) tmpy.push_back(0.);
    for(int i=0; i<n; i++) tmpz.push_back(0.);

    //first rk4 step
    rhs(t, x, kx1, n);
    rhs(t, y, ky1, n);
    rhs(t, z, kz1, n);
    for(int i=0; i<n; i++) tmpx[i] = x[i] + (dt/2.) * kx1[i];
    for(int i=0; i<n; i++) tmpy[i] = y[i] + (dt/2.) * ky1[i];
    for(int i=0; i<n; i++) tmpz[i] = z[i] + (dt/2.) * kz1[i];

    //second rk4 step
    rhs(t+dt/2., tmpx, kx2, n);
    rhs(t+dt/2., tmpy, ky2, n);
    rhs(t+dt/2., tmpz, kz2, n);
    for(int i=0; i<n; i++) tmpx[i] = x[i] + (dt/2.) * kx2[i];
    for(int i=0; i<n; i++) tmpy[i] = y[i] + (dt/2.) * ky2[i];
    for(int i=0; i<n; i++) tmpz[i] = z[i] + (dt/2.) * kz2[i];

    //third rk4 step
    rhs(t+dt/2., tmpx, kx3, n);
    rhs(t+dt/2., tmpy, ky3, n);
    rhs(t+dt/2., tmpz, kz3, n);
    for(int i=0; i<n; i++) tmpx[i] = x[i] + dt * kx3[i];
    for(int i=0; i<n; i++) tmpy[i] = y[i] + dt * ky3[i];
    for(int i=0; i<n; i++) tmpz[i] = z[i] + dt * kz3[i];

    //fourth rk4 step
    rhs(t+dt, tmpx, kx4, n);
    rhs(t+dt, tmpy, ky4, n);
    rhs(t+dt, tmpz, kz4, n);

    //do the iteration step (update the positions)
    for(int i=0; i<n; i++) x[i] += (dt/6.) * (kx1[i] + 2.*kx2[i] + 2.*kx3[i] + kx4[i]);
    for(int i=0; i<n; i++) y[i] += (dt/6.) * (ky1[i] + 2.*ky2[i] + 2.*ky3[i] + ky4[i]);
    for(int i=0; i<n; i++) z[i] += (dt/6.) * (kz1[i] + 2.*kz2[i] + 2.*kz3[i] + kz4[i]);
}

void lf_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, void* rhs(double t, vector<double> x, vector<double> &u_rhs, int n), int n){
    //Initialize vectors for the steps - only one step here!
    vector<double> kvx1, kvy1, kvz1;

    for(int i=0; i<n; i++) kvx1.push_back(0.);
    for(int i=0; i<n; i++) kvy1.push_back(0.);
    for(int i=0; i<n; i++) kvz1.push_back(0.);

    //lf step - calculate derivative of v
    rhs(t, x, kvx1, n);         
    rhs(t, y, kvy1, n);
    rhs(t, z, kvz1, n);

    //calculate n+1/2 value of v
    for(int i=0; i<n; i++) vx[i] += dt * kvx1[i];  
    for(int i=0; i<n; i++) vy[i] += dt * kvy1[i];
    for(int i=0; i<n; i++) vz[i] += dt * kvz1[i];

    //do the iteration step (update the positions (n+1))
    for(int i=0; i<n; i++) x[i] += dt * vx[i];
    for(int i=0; i<n; i++) y[i] += dt * vy[i];
    for(int i=0; i<n; i++) z[i] += dt * vz[i];
}

void driver_fwd(double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, int n){
    //Create and open output file
    fstream file;
    file.open("fwd-solution.csv", ios::out);
    file.precision(10);

    //loop that iterates up to a certain chosen time (end)
    while((t_end - t) > DBL_EPSILON){
        //Output current values to file - "; " is needed as delimiter for cells
        //Iterations are needed to generally output for n objects without adjusting anything
        file << t << "; ";
            for(int i=0; i<n; i++) file << x[i] << "; ";
            for(int i=0; i<n; i++) file << y[i] << "; ";
            for(int i=0; i<n; i++) file << z[i] << "; ";
            for(int i=0; i<n; i++) file << vx[i] << "; ";
            for(int i=0; i<n; i++) file << vy[i] << "; ";
            for(int i=0; i<n-1; i++) file << vz[i] << "; ";
        file << vz[n-1] << endl;

        //Calculate next timestep
        fwd_step(t, dt, x, y, z, vx, vy, vz, testfunction, n);

        //update time - so the loop will have a chance to end
        t += dt;
    }

    //close the output file after the iterations are done
    file.close();
}

void driver_rk4(double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, int n){
    //Create and open output file
    fstream file;
    file.open("rk4-solution.csv", ios::out);
    file.precision(10);

    //loop that iterates up to a certain chosen time (end)
    while((t_end - t) > DBL_EPSILON){
        //Output current values to file - "; " is needed as delimiter for cells
        //Iterations are needed to generally output for n objects without adjusting anything
        file << t << "; ";
            for(int i=0; i<n; i++) file << x[i] << "; ";
            for(int i=0; i<n; i++) file << y[i] << "; ";
            for(int i=0; i<n; i++) file << z[i] << "; ";
            for(int i=0; i<n; i++) file << vx[i] << "; ";
            for(int i=0; i<n; i++) file << vy[i] << "; ";
            for(int i=0; i<n-1; i++) file << vz[i] << "; ";
        file << vz[n-1] << endl;

        //Calculate next timestep
        rk4_step(t, dt, x, y, z, vx, vy, vz, testfunction, n);

        //update time - so the loop will have a chance to end
        t += dt;
    }

    //close the output file after the iterations are done
    file.close();   
}

void driver_lf(double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, int n){
    //Create and open output file
    fstream file;
    file.open("lf-solution.csv", ios::out);
    file.precision(10);

    //loop that iterates up to a certain chosen time (end)
    while((t_end - t) > DBL_EPSILON){
        //Output current values to file - "; " is needed as delimiter for cells
        //Iterations are needed to generally output for n objects without adjusting anything
        file << t << "; ";
            for(int i=0; i<n; i++) file << x[i] << "; ";
            for(int i=0; i<n; i++) file << y[i] << "; ";
            for(int i=0; i<n; i++) file << z[i] << "; ";
            for(int i=0; i<n; i++) file << vx[i] << "; ";
            for(int i=0; i<n; i++) file << vy[i] << "; ";
            for(int i=0; i<n-1; i++) file << vz[i] << "; ";
        file << vz[n-1] << endl;

        //Calculate next timestep
        lf_step(t, dt, x, y, z, vx, vy, vz, testsymplectic, n);

        //update time - so the loop will have a chance to end
        t += dt;
    }

    //close the output file after the iterations are done
    file.close();
}

void programmteil(string command, vector<string> &commands, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz){
    commands.resize(3);
    commands[0] = "fwd";    //forward euler
    commands[1] = "rk4";    //runge kutta 4
    commands[2] = "lf";     //leapfrog

    int n = 2;              //Number of particles
    double t_end = 2*M_PI;
    double dt = pow(10,-2);
    double t = 0.;

    if (command == commands[0]){
        initialize(n, x, y, z, vx, vy, vz);
        driver_fwd(t, t_end, dt, x, y, z, vx, vy, vz, n);
    }

    if (command == commands[1]){
        initialize(n, x, y, z, vx, vy, vz);
        driver_rk4(t, t_end, dt, x, y, z, vx, vy, vz, n);
    }

    if (command == commands[2]){
        initialize_symplectic(n, x, y, z, vx, vy, vz);
        driver_lf(t, t_end, dt, x, y, z, vx, vy, vz, n);
    }

    if ( (command != commands[0]) && (command != commands[1]) && (command != commands[2])) cout << "Wrong parameter!" << endl;
}

int main(int argc, char** argv){
    auto t1 = chrono::high_resolution_clock::now();
    //arguments of main used to hand over parameters at the programm start
    //1st parameter is programm name
    if(argc != 2) {
        cout << "usage:\n" << argv[0] << " <parameter> || possible parameters: 'fwd', 'rk4' or 'lf' " << endl;
        return -1;
    }
    else{
        stringstream input{argv[1]};
        string command;
        input >> command;
        vector<string> commands = {};

        //vectors chosen such that n-particles can be realized
        //Only current values are stored and after the output overridden with the new ones
        vector<double> x = {};
        vector<double> y = {};
        vector<double> z = {};

        vector<double> vx = {};
        vector<double> vy = {};
        vector<double> vz = {};

        programmteil(command, commands, x, y, z, vx, vy, vz);
        auto t2 = chrono::high_resolution_clock::now();
        auto time = chrono::duration<float>(t2-t1).count();
        cout << "Reached end of main." << endl;
        cout << "Die Berechnung hat " << time << " Sekunden gedauert." << endl;
        
        return 0;
    }
}

/* To do:
 Programming part:
    -Write .csv/.txt -file with initialvalues for ~ 20-30 objects we want to consider with lines like: x, y, z, vx, vy, vz, m;
    -Write function that reads the initialvalue file and uses n objects (first n lines)
    -Test 1d integrators for different functions and determine order of convergence with log-log-plot
    -update testfunctions to newtons grav. law
    -update stepfunctions to deal with d > 2 dimensions
    -mass-vector in main
    -calculate E_tot, E_pot, E_kin in driver functions and output them for conservation plot
    -calculate L_tot in driver functions and output it for conservation plot
    -calculate Laplace integral in driver functions and output it for conservation plot
    -implement adaptive stepsize control in the driver functions
    -Write a script for plots:
        -Energy
        -L
        -Laplace integral
        -trajectories in position space
        -trajectories in phase space
        (-trajectories in momentum space)
    -mabye unify the driver functions, since they are pretty much the same

 Theoretical foundations:
    -datatypes and their representation
    -numerical errors: their occurances and reduction
    -idea of integrator schemes: from taylor-expansion to final formulas
    -representation of formulas in butcher tables and how to implement code from there
    -determining order of convergence (power laws) from log-log-plots

    -from newtons law / keplers law to final functions for acceleration calculations
    (-post-newtonian/relativistic corrections)
    -transformation of data/initial data in order to avoid poorly scaling prefactors like G (pow(10,-11))

 Other challenges:
    -how to call c++ programms / functions from a python script
    -how to unify the plotting python scrips to a complete script that does programm execution, numerical analysis and plotting
    -how to handle auxiliary functions of symplectic integrators if they appear in the output at all, since they are off by 0.5 time step

 Questions:
    -can we reduce calculations and ensure that the 3 Quantities are conserved by calculating them and adopt one planet
    otherwise we had some redundant information.
    -if we can, would this or a similar procedure be stable since small deviations from all n-1 objects summed up and put to one object could mabye cause problems.
*/