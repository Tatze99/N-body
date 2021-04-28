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

typedef void (* Step_function)(double, double, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>(double, vector<double>, vector<double>, vector<double>, int, vector<double>), int, vector<double>);

typedef vector<double> (DGL)(double t, vector<double> x, vector<double> y, vector<double> z, int n, vector<double> m);

// Testfunctions ----------------------------------------------------------------------------
vector<double> acceleration(double t, vector<double> x, vector<double> y, vector<double> z, int n, vector<double> m){
  double Matrix[n][n];

  vector<double> a(n,0.0);
  for(int i=0; i<n; i++) a.push_back(0.);

  for(int i=0; i<n; i++){
    for(int j=0; j<i; j++) {
      Matrix[i][j] = (x[i]-x[j])/pow((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]),1.5);
    }
  }

  for(int i=0; i<n; i++){
    Matrix[i][i] = 0.;
    for(int j=i+1; j<n; j++) Matrix[i][j] = -Matrix[j][i];
    for(int j=0; j<n; j++) {
      a[i] += Matrix[i][j]*m[j];
    }
    a[i] *= -4*M_PI*M_PI;
  }
  return a;
}

vector<double> cosine(double t, vector<double> x, vector<double> y, vector<double> z, int n, vector<double> m){
  vector<double> u_rhs(n,0.0);
  for(int i=0; i<n; i++) u_rhs[i] = cos(t);
}

vector<double> testsymplectic(double t, vector<double> x, vector<double> y, vector<double> z, int n, vector<double> m){
  vector<double> u_rhs(n,0.0);
  for(int i=0; i<n; i++) u_rhs[i] = -sin(t);
}

// Initialize variables ----------------------------------------------------------------------------
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
    x[1] = 1.;
    y[0] = 0.;
    y[1] = 0.;
    z[0] = 0.;
    z[1] = 0.;

    //fwd and rk4 method need initialisation of v at the same timestep as x
    vx[0] = 0.;
    vx[1] = 0.;
    vy[0] = 0.;
    vy[1] = 2*M_PI;
    vz[0] = 0.;
    vz[1] = 0.;

    //mass [in units of the mass of the sun]
    m[0] = 1.;
    m[1] = 1/333000.;
}

void initialize_symplectic(int n, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, double dt, vector<double> &m){
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
    x[1] = 1.;
    y[0] = 0.;
    y[1] = 0.;
    z[0] = 0.;
    z[1] = 0.;
    
    //contrary the lf method needs initialisation of v a half timestep before the first timestep of x
    //fwd and rk4 method need initialisation of v at the same timestep as x
    vx[0] = 0.;
    vx[1] = 0.;
    vy[0] = 0.;
    vy[1] = 2*M_PI*(1-dt);
    vz[0] = 0.;
    vz[1] = 0.;

    //mass [in units of the mass of the sun]
    m[0] = 1.;
    m[1] = 1/333000.;
}

//Functions for read ----------------------------------------------------------------------------
inline bool fileexists (const std::string& name) {
    if (FILE* file = fopen(name.c_str(), "r")) {
        fclose(file);
        free(file);
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

   /*
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
    */

   //Set

    //set breaking command for printing and line counting functions
    sngspecific.push_back(s);
    data.push_back(s);
}

//Step functions ----------------------------------------------------------------------------
void fwd_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, DGL rhs, int n, vector<double> m){
    //Initialize vectors for the steps - only one step here!
    vector<double> kx(n), ky(n), kz(n);

    //determine the derivative in every direction (calculate kx, ky, kz)
    kx = rhs(t, x, y, z, n, m);
    ky = rhs(t, y, z, x, n, m);
    kz = rhs(t, z, x, y, n, m);

    //do the iteration step (update the positions)
    for(int i=0; i<n; i++) x[i] += dt * vx[i];
    for(int i=0; i<n; i++) y[i] += dt * vy[i];
    for(int i=0; i<n; i++) z[i] += dt * vz[i];
    for(int i=0; i<n; i++) vx[i] += dt * kx[i];
    for(int i=0; i<n; i++) vy[i] += dt * ky[i];
    for(int i=0; i<n; i++) vz[i] += dt * kz[i];
}

void rk4_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, DGL rhs, int n, vector<double> m){

    //Initialize vectors for the steps - only one step here!
    vector<double> kx1(n), kx2(n), kx3(n), kx4(n), tmpx(n);
    vector<double> ky1(n), ky2(n), ky3(n), ky4(n), tmpy(n);
    vector<double> kz1(n), kz2(n), kz3(n), kz4(n), tmpz(n);
    vector<double> vx1(n), vx2(n), vx3(n), vx4(n);
    vector<double> vy1(n), vy2(n), vy3(n), vy4(n);
    vector<double> vz1(n), vz2(n), vz3(n), vz4(n);

    //first rk4 step
    for(int i=0; i<n; i++) vx1[i] = vx[i];
    for(int i=0; i<n; i++) vy1[i] = vy[i];
    for(int i=0; i<n; i++) vz1[i] = vz[i];
    kx1 = rhs(t, x, y, z, n, m);
    ky1 = rhs(t, y, z, x, n, m);
    kz1 = rhs(t, z, x, y, n, m);
    //second rk4 step
    for(int i=0; i<n; i++) vx2[i] = vx[i] + (dt/2.) * kx1[i];
    for(int i=0; i<n; i++) vy2[i] = vy[i] + (dt/2.) * ky1[i];
    for(int i=0; i<n; i++) vz2[i] = vz[i] + (dt/2.) * kz1[i];

    for(int i=0; i<n; i++) tmpx[i] = x[i] + (dt/2.) * vx1[i];
    for(int i=0; i<n; i++) tmpy[i] = y[i] + (dt/2.) * vy1[i];
    for(int i=0; i<n; i++) tmpz[i] = z[i] + (dt/2.) * vz1[i];
    kx2 = rhs(t+dt/2., tmpx, tmpy, tmpz, n, m);
    ky2 = rhs(t+dt/2., tmpy, tmpz, tmpx, n, m);
    kz2 = rhs(t+dt/2., tmpz, tmpx, tmpy, n, m);
    //third rk4 step
    for(int i=0; i<n; i++) vx3[i] = vx[i] + (dt/2.) * kx2[i];
    for(int i=0; i<n; i++) vy3[i] = vy[i] + (dt/2.) * ky2[i];
    for(int i=0; i<n; i++) vz3[i] = vz[i] + (dt/2.) * kz2[i];

    for(int i=0; i<n; i++) tmpx[i] = x[i] + (dt/2.) * vx2[i];
    for(int i=0; i<n; i++) tmpy[i] = y[i] + (dt/2.) * vy2[i];
    for(int i=0; i<n; i++) tmpz[i] = z[i] + (dt/2.) * vz2[i];
    kx3 = rhs(t+dt/2., tmpx, tmpy, tmpz, n, m);
    ky3 = rhs(t+dt/2., tmpy, tmpz, tmpx, n, m);
    kz3 = rhs(t+dt/2., tmpz, tmpx, tmpy, n, m);
    //fourth rk4 step
    for(int i=0; i<n; i++) vx4[i] = vx[i] + dt * kx3[i];
    for(int i=0; i<n; i++) vy4[i] = vy[i] + dt * ky3[i];
    for(int i=0; i<n; i++) vz4[i] = vz[i] + dt * kz3[i];

    for(int i=0; i<n; i++) tmpx[i] = x[i] + dt * vx3[i];
    for(int i=0; i<n; i++) tmpy[i] = y[i] + dt * vy3[i];
    for(int i=0; i<n; i++) tmpz[i] = z[i] + dt * vz3[i];
    kx4 = rhs(t+dt, tmpx, tmpy, tmpz, n, m);
    ky4 = rhs(t+dt, tmpy, tmpz, tmpx, n, m);
    kz4 = rhs(t+dt, tmpz, tmpx, tmpy, n, m);

    //do the iteration step (update the positions)
    for(int i=0; i<n; i++) vx[i] += (dt/6.) * (kx1[i] + 2.*kx2[i] + 2.*kx3[i] + kx4[i]);
    for(int i=0; i<n; i++) vy[i] += (dt/6.) * (ky1[i] + 2.*ky2[i] + 2.*ky3[i] + ky4[i]);
    for(int i=0; i<n; i++) vz[i] += (dt/6.) * (kz1[i] + 2.*kz2[i] + 2.*kz3[i] + kz4[i]);
    for(int i=0; i<n; i++)  x[i] += (dt/6.) * (vx1[i] + 2.*vx2[i] + 2.*vx3[i] + vx4[i]);
    for(int i=0; i<n; i++)  y[i] += (dt/6.) * (vy1[i] + 2.*vy2[i] + 2.*vy3[i] + vy4[i]);
    for(int i=0; i<n; i++)  z[i] += (dt/6.) * (vz1[i] + 2.*vz2[i] + 2.*vz3[i] + vz4[i]);

}

void lf_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, DGL rhs, int n, vector<double> m){
    //Initialize vectors for the steps - only one step here!
    vector<double> kx(n), ky(n), kz(n);

    //lf step - calculate derivative of v
    kx = rhs(t, x, y, z, n, m);
    ky = rhs(t, y, z, x, n, m);
    kz = rhs(t, z, x, y, n, m);

    //calculate n+1/2 value of v
    for(int i=0; i<n; i++) vx[i] += dt * kx[i];
    for(int i=0; i<n; i++) vy[i] += dt * ky[i];
    for(int i=0; i<n; i++) vz[i] += dt * kz[i];
    //do the iteration step (update the positions (n+1))
    for(int i=0; i<n; i++) x[i] += dt * vx[i];
    for(int i=0; i<n; i++) y[i] += dt * vy[i];
    for(int i=0; i<n; i++) z[i] += dt * vz[i];
}

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
            step(t, dt, x, y, z, vx, vy, vz, acceleration, n, m);

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
            step(t, dt, x, y, z, vx, vy, vz, acceleration, n, m);

            //update time - so the loop will have a chance to end
            t += dt;
        }
    }

    //close the output file after the iterations are done
    file.close();
}

void programmteil(string command, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> &m){
    int n = 2;              //Number of particles
    double t_end = 1;       //final time
    double dt = pow(10,-4); //time steps
    double t = 0.;

    if (command == "fwd"){  // forward euler
        initialize(n, x, y, z, vx, vy, vz, m);
        driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, fwd_step, command);
    }
    else if (command == "rk4"){ // Runge Kutta 4
        initialize(n, x, y, z, vx, vy, vz, m);
        driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, rk4_step, command);
    }
    else if (command == "lf"){ // leap frog
        initialize_symplectic(n, x, y, z, vx, vy, vz, dt, m);
        driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, fwd_step, command);
    }
    else cout << "Wrong parameter!" << endl;
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

        //vectors chosen such that n-particles can be realized
        //Only current values are stored and after the output overridden with the new ones
        vector<double> x = {};
        vector<double> y = {};
        vector<double> z = {};

        vector<double> vx = {};
        vector<double> vy = {};
        vector<double> vz = {};

        vector<double> m = {};

        programmteil(command, x, y, z, vx, vy, vz, m);
        auto t2 = chrono::high_resolution_clock::now();
        auto time = chrono::duration<float>(t2-t1).count();
        cout << "Reached end of main." << endl;
        cout << "Die Berechnung hat " << time << " Sekunden gedauert." << endl;

        return 0;
    }
}
