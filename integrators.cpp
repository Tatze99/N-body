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

//#include<boost/python.hpp>


//for global "short-hand" notation - need not to write 'std::' in front of most things
using namespace std;

typedef void (*Step_function)(double, double, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>, void* (*)(double, vector<double>, vector<double>, vector<double>&, int), int);
//const double PI = 4.*atan(1.);

//Functions for read ----------------------------------------------------------------------------

inline bool fileexists (const string& name) {
    //fstream file;
    //file.open(name.c_str(), ios::in)
    if (FILE* file = fopen(name.c_str(), "r")) {
        fclose(file);
        free(file);
        return true;
    } else {
        return false;
    }
}

//void read_file(vector<string> &v, string &filename, string &str){
void read_file(vector<string> &v, string &filename){
    //read complete file and store to auxiliary vector v
    //create file s
    fstream s;
    char cstring[256];

    //open file and read data
    s.open(filename, ios::in);

    while (!s.eof()){
        s.getline(cstring, sizeof(cstring));
        v.push_back(cstring);
    }

    s.close();
}

//void seperate_to_files(vector<string> &file, vector<string> &sngspecific, vector<string> &data, string &s){
void seperate_to_files(int n, vector<string> help, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> &m){
   //create from vector (file) with whole file two seperate files containing songbeamber specific information (sngspecific) or song data (data)
   int counter;
   string s, tmp_s;
   //s = "012;3456789";
   string semi = ";"; 
 //int i = 0;
   //iterate over auxiliary vector
    for(int i=0; i<n; i++){
        s = help[i];

        for(int l=0; l<7; l++){
            //iterate over each line --- mabye replace by find_first_of (see cppreference)
            for(int j=0; j<s.length(); j++){
                //cout << s[j];
                if(s[j] == semi[0]){
                    counter = j;
                    break;
                }
            }
            cout << "l = " << l << "; counter = " << counter << endl;
            //cut off the first part of s 
            tmp_s = s.substr(0, counter);

            //set s to the remaining string
            if ((counter-2) < (s.length()-1)) s = s.substr(counter+2, s.length()-1);

            //set values
            switch (l)
            {
                case 0: x[i] = stod(tmp_s);
                        break;
                case 1: y[i] = stod(tmp_s);
                        break;
                case 2: z[i] = stod(tmp_s);
                        break;
                case 3: vx[i] = stod(tmp_s);
                        break;
                case 4: vy[i] = stod(tmp_s);
                        break;
                case 5: vz[i] = stod(tmp_s);
                        break;
                case 6: m[i] = stod(tmp_s);
            }  
            cout << "Test " << l << endl;         
        }
    }
}

//Functions for read ----------------------------------------------------------------------------

void initialize(int n, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> &m){
    //Change size of vectors
    x.resize(n);
    y.resize(n);
    z.resize(n);
    vx.resize(n);
    vy.resize(n);
    vz.resize(n);
    m.resize(n);

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
    string name = "Input.csv";
    vector<string> help = {};

    read_file(help, name);

    seperate_to_files(n, help, x, y, z, vx, vy, vz, m);

    for(double i : x) cout << i << "; ";
    cout << endl;
    for(double i : y) cout << i << "; ";
    cout << endl;

    //cout << x[0] << "; " << y[0] << endl;
    //mabye initialize natural constants
}

void initialize_symplectic(int n, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> &m){
    //Change size of vectors
    x.resize(n);
    y.resize(n);
    z.resize(n);
    vx.resize(n);
    vy.resize(n);
    vz.resize(n);
    m.resize(n);

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

void *testfunction(double t, vector<double> x, vector<double> m, vector<double> &u_rhs, int n){
    int i;
    u_rhs.resize(n);
    for(int i=0; i<n; i++) u_rhs[i] = cos(t);
}

void *testsymplectic(double t, vector<double> x, vector<double> m, vector<double> &u_rhs, int n){
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

void fwd_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> m, void* rhs(double t, vector<double> x, vector<double> m, vector<double> &u_rhs, int n), int n){
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
    rhs(t, x, m, kx, n);
    rhs(t, y, m, ky, n);
    rhs(t, z, m, kz, n);

    //do the iteration step (update the positions)
    for(int i=0; i<n; i++) x[i] += dt * kx[i];
    for(int i=0; i<n; i++) y[i] += dt * ky[i];
    for(int i=0; i<n; i++) z[i] += dt * kz[i];
    //for(int i=0; i<n; i++) vx[i] += dt * kvx;
    //for(int i=0; i<n; i++) vy[i] += dt * kvy;
    //for(int i=0; i<n; i++) vz[i] += dt * kvz;
}

void rk4_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> m, void* rhs(double t, vector<double> x, vector<double> m, vector<double> &u_rhs, int n), int n){
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
    rhs(t, x, m, kx1, n);
    rhs(t, y, m, ky1, n);
    rhs(t, z, m, kz1, n);
    for(int i=0; i<n; i++) tmpx[i] = x[i] + (dt/2.) * kx1[i];
    for(int i=0; i<n; i++) tmpy[i] = y[i] + (dt/2.) * ky1[i];
    for(int i=0; i<n; i++) tmpz[i] = z[i] + (dt/2.) * kz1[i];

    //second rk4 step
    rhs(t+dt/2., tmpx, m, kx2, n);
    rhs(t+dt/2., tmpy, m, ky2, n);
    rhs(t+dt/2., tmpz, m, kz2, n);
    for(int i=0; i<n; i++) tmpx[i] = x[i] + (dt/2.) * kx2[i];
    for(int i=0; i<n; i++) tmpy[i] = y[i] + (dt/2.) * ky2[i];
    for(int i=0; i<n; i++) tmpz[i] = z[i] + (dt/2.) * kz2[i];

    //third rk4 step
    rhs(t+dt/2., tmpx, m, kx3, n);
    rhs(t+dt/2., tmpy, m, ky3, n);
    rhs(t+dt/2., tmpz, m, kz3, n);
    for(int i=0; i<n; i++) tmpx[i] = x[i] + dt * kx3[i];
    for(int i=0; i<n; i++) tmpy[i] = y[i] + dt * ky3[i];
    for(int i=0; i<n; i++) tmpz[i] = z[i] + dt * kz3[i];

    //fourth rk4 step
    rhs(t+dt, tmpx, m, kx4, n);
    rhs(t+dt, tmpy, m, ky4, n);
    rhs(t+dt, tmpz, m, kz4, n);

    //do the iteration step (update the positions)
    for(int i=0; i<n; i++) x[i] += (dt/6.) * (kx1[i] + 2.*kx2[i] + 2.*kx3[i] + kx4[i]);
    for(int i=0; i<n; i++) y[i] += (dt/6.) * (ky1[i] + 2.*ky2[i] + 2.*ky3[i] + ky4[i]);
    for(int i=0; i<n; i++) z[i] += (dt/6.) * (kz1[i] + 2.*kz2[i] + 2.*kz3[i] + kz4[i]);
}

void lf_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> m, void* rhs(double t, vector<double> x, vector<double> m, vector<double> &u_rhs, int n), int n){
    //Initialize vectors for the steps - only one step here!
    vector<double> kvx1, kvy1, kvz1;

    for(int i=0; i<n; i++) kvx1.push_back(0.);
    for(int i=0; i<n; i++) kvy1.push_back(0.);
    for(int i=0; i<n; i++) kvz1.push_back(0.);

    //lf step - calculate derivative of v
    rhs(t, x, m, kvx1, n);
    rhs(t, y, m, kvy1, n);
    rhs(t, z, m, kvz1, n);

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

void driver(double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, int n, vector<double> m, Step_function step, string command){

    //Create and open output file
    fstream file;
    file.open(command+"-solution.csv", ios::out);
    file.precision(10);

    double e_kin, e_pot, e_tot;

    //For the LF integrator we will possibly end up with another
    //function to integrate than for the fwd and rk4 method
    //Wrong function was a problem for the integrator scheme
    //Problem is resolved by the following if-construction
    if(command=="lf"){
        //loop that iterates up to a certain chosen time (end)
        while((t_end - t) > DBL_EPSILON){
            e_kin = 0.;
            e_pot = 0.;
            for(int i; i<n; i++) e_kin += 0.5 * m[i] * (pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2));
            //for(int i; i<n; i++) e_pot -= m[i] * acceleration[i] * sqrt(pow(x[i],2) + pow(y[i],2) + pow(z[i],2));
            e_tot = e_kin + e_pot;
            
            //Output current values to file - "; " is needed as delimiter for cells
            //Iterations are needed to generally output for n objects without adjusting anything
            file << t << "; ";
                for(int i=0; i<n; i++) file << x[i] << "; ";
                for(int i=0; i<n; i++) file << y[i] << "; ";
                for(int i=0; i<n; i++) file << z[i] << "; ";
                for(int i=0; i<n; i++) file << vx[i] << "; ";
                for(int i=0; i<n; i++) file << vy[i] << "; ";
                for(int i=0; i<n; i++) file << vz[i] << "; ";
            file << e_kin << "; " << e_pot << "; " << e_tot << endl;

            //Calculate next timestep
            step(t, dt, x, y, z, vx, vy, vz, m, testsymplectic, n);

            //update time - so the loop will have a chance to end
            t += dt;
        }
    }
    else{
        //loop that iterates up to a certain chosen time (end)
        while((t_end - t) > DBL_EPSILON){
            e_kin = 0.;
            e_pot = 0.;
            for(int i; i<n; i++) e_kin += 0.5 * m[i] * (pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2));
            //for(int i; i<n; i++) e_pot -= m[i] * acceleration[i] * sqrt(pow(x[i],2) + pow(y[i],2) + pow(z[i],2));
            e_tot = e_kin + e_pot;
            
            //Output current values to file - "; " is needed as delimiter for cells
            //Iterations are needed to generally output for n objects without adjusting anything
            file << t << "; ";
                for(int i=0; i<n; i++) file << x[i] << "; ";
                for(int i=0; i<n; i++) file << y[i] << "; ";
                for(int i=0; i<n; i++) file << z[i] << "; ";
                for(int i=0; i<n; i++) file << vx[i] << "; ";
                for(int i=0; i<n; i++) file << vy[i] << "; ";
                for(int i=0; i<n; i++) file << vz[i] << "; ";
            file << e_kin << "; " << e_pot << "; " << e_tot << endl;

            //Calculate next timestep
            step(t, dt, x, y, z, vx, vy, vz, m, testfunction, n);

            //update time - so the loop will have a chance to end
            t += dt;
        }
    }
 
    //close the output file after the iterations are done
    file.close();
}

void programmteil(string command, vector<string> &commands, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> &m){
    commands.resize(3);
    commands[0] = "fwd";    //forward euler
    commands[1] = "rk4";    //runge kutta 4
    commands[2] = "lf";     //leapfrog

    int n = 2;              //Number of particles
    double t_end = 2*M_PI;
    double dt = pow(10,-2);
    double t = 0.;

    if (command == commands[0]){
        initialize(n, x, y, z, vx, vy, vz, m);
        //driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, fwd_step, command);
        // driver_fwd(t, t_end, dt, x, y, z, vx, vy, vz, n);
    }

    if (command == commands[1]){
        initialize(n, x, y, z, vx, vy, vz, m);
        driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, rk4_step, command);
        // driver_rk4(t, t_end, dt, x, y, z, vx, vy, vz, n);
    }

    if (command == commands[2]){
        initialize_symplectic(n, x, y, z, vx, vy, vz, m);
        driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, lf_step, command);
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

        vector<double> m = {};

        programmteil(command, commands, x, y, z, vx, vy, vz, m);
        auto t2 = chrono::high_resolution_clock::now();
        auto time = chrono::duration<float>(t2-t1).count();
        cout << "Reached end of main." << endl;
        cout << "Die Berechnung hat " << time << " Sekunden gedauert." << endl;

        return 0;
    }
}
