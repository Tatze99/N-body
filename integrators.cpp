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
#include <chrono>
#include <cfloat>


//for global "short-hand" notation - need not to write 'std::' in front of most things
using namespace std;

//global variables
//const int n = 2;

typedef void (* Step_function)(double, double, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, void* (*)(double, vector<double>, vector<double>&, int), int);
//const double PI = 4.*atan(1.);

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

//Functions for read ----------------------------------------------------------------------------

inline bool fileexists (const string& name) {
    if (fstream file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

void read_file(vector<string> &v, string &filename, string &str){
    //read complete file and store to vector v
    //create vector s
    fstream s;
    char cstring[256];
    string tmp;
    string tmp1;
    int tmpascii;

    //open file and read data
    s.open(filename, ios::in);

    while (!s.eof())
    {
        //tmp = "";
        s.getline(cstring, sizeof(cstring));
        //tmp1 = cstring;
        //for(int i=0; i<tmp1.length();i++){
        //    tmpascii = tmp1[i];
        //    //ascii-sign 015 is some sort of newline command in CCLI or songbeamer files
        //    if (tmpascii != 015) tmp += cstring[i];
        //}
        //v.push_back(tmp);
        v.push_back(cstring);
    }

    //set breaking command for printing and line counting functions
    //v.push_back(str);
    s.close();
}

void seperate_to_files(vector<string> &file, vector<string> &sngspecific, vector<string> &data, string &s){
   //create from vector (file) with whole file two seperate files containing songbeamber specific information (sngspecific) or song data (data)
   using namespace std;
   
   int i = 0;
   string sngtag = "#";
   string snginfo;

    //go through all lines of file and sort them
    while (file[i].compare(s) != 0){
        snginfo = file[i];

        //sort songbeamer specific information to sngspecific
        if (snginfo[0] == sngtag[0]){
            sngspecific.push_back(file[i]);
            i++;
            continue;
        }

        //sort song data to data
        data.push_back(file[i]);
        i++;
    }

    //set breaking command for printing and line counting functions
    sngspecific.push_back(s);
    data.push_back(s);
}

//Functions for read ----------------------------------------------------------------------------

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

/*
//Function for speed calculation
	double speed(double t, double v_alt, double a) {
		double v_neu = v_alt + a * t;
			return v_neu;
		}

	//Function for acceleration calculation
	double acc(double x, double y) {
		double a = -4 * pow(M_PI, 2) * x / pow(hypot(x,y), 3);
		return a;
	}

void *speed(){

}


void *acc(double t, vector<double> x, vector<double> y, vector<double> z, vector<double> &u_rhs){
    u_rhs = 0.;
    for(int i=0; i<n; i++) u_rhs[i] += -4 * pow(M_PI, 2) * x[i] / pow(hypot(x,y), 3);
}

*/

void fwd_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, void* rhs(double t, vector<double> x, vector<double> &u_rhs, int n), int n){
    //Initialize vectors for the steps - only one step here!
    vector<double> kx, ky, kz;
    //vector<double> kvx, kvy, kvz;

    for(int i=0; i<n; i++) kx.push_back(0.);
    for(int i=0; i<n; i++) ky.push_back(0.);
    for(int i=0; i<n; i++) kz.push_back(0.);
    //for(int i=0; i<n; i++) kvx.push_back(0.);
    //for(int i=0; i<n; i++) kvy.push_back(0.);
    //for(int i=0; i<n; i++) kvz.push_back(0.);

    //determine the derivative in every direction (calculate kx, ky, kz)
    rhs(t, x, kx, n);
    rhs(t, y, ky, n);
    rhs(t, z, kz, n);

    //do the iteration step (update the positions)
    for(int i=0; i<n; i++) x[i] += dt * kx[i];
    for(int i=0; i<n; i++) y[i] += dt * ky[i];
    for(int i=0; i<n; i++) z[i] += dt * kz[i];
    //for(int i=0; i<n; i++) vx[i] += dt * kvx;
    //for(int i=0; i<n; i++) vy[i] += dt * kvy;
    //for(int i=0; i<n; i++) vz[i] += dt * kvz;
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

/*
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
            for(int i=0; i<n; i++) file << vz[i] << "; ";
        file << endl;

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
            // for(int i=0; i<n; i++) file << vz[i] << "; ";  // kein Semicolon am Ende der Zeile, sonst kann Python die Daten nicht einlesen.
            for(int i=0; i<n; i++) file << vz[i];
        file << endl;

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
            for(int i=0; i<n; i++) file << vz[i] << "; ";
        file << endl;

        //Calculate next timestep
        lf_step(t, dt, x, y, z, vx, vy, vz, testsymplectic, n);

        //update time - so the loop will have a chance to end
        t += dt;
    }

    //close the output file after the iterations are done
    file.close();
}
*/

void driver(double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, int n, Step_function step){

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
        step(t, dt, x, y, z, vx, vy, vz, testfunction, n);

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
    double t_end = 10*M_PI;
    //double t_end = 10*PI;
    double dt = pow(10,-3);
    double t = 0.;

    if (command == commands[0]){
        initialize(n, x, y, z, vx, vy, vz);
        driver(t, t_end, dt, x, y, z, vx, vy, vz, n, fwd_step);
        // driver_fwd(t, t_end, dt, x, y, z, vx, vy, vz, n);
    }

    if (command == commands[1]){
        initialize(n, x, y, z, vx, vy, vz);
        driver(t, t_end, dt, x, y, z, vx, vy, vz, n, rk4_step);
        // driver_rk4(t, t_end, dt, x, y, z, vx, vy, vz, n);
    }

    if (command == commands[2]){
        initialize_symplectic(n, x, y, z, vx, vy, vz);
        driver(t, t_end, dt, x, y, z, vx, vy, vz, n, lf_step);
        // driver_lf(t, t_end, dt, x, y, z, vx, vy, vz, n);
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
