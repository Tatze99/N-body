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
#include <tuple>
#include <omp.h>

#include "header.cpp"
//for global "short-hand" notation - need not to write 'std::' in front of most things
using namespace std;

vector<double> all_from_target(vector<string> &tmp, string input, double objdist, double time){

    int s = tmp.size();
    vector<double> t_max(s), r_max(s), v_max(s), v_0(s);
    vector<double> all = {0., 0., 0., 0., 0.};
    int numbers = s-1;

    vector<vector<double>> values = {t_max,r_max,v_max,v_0};
    values = set_values(values, numbers, input);

    t_max = values[0];
    r_max = values[1];
    v_max = values[2];
    v_0 = values[3];

    int counter = 1;

    for(int i=1; i<s; i++){
        if( (objdist - r_max[i]) < DBL_EPSILON){
            counter = i;
            break;
        }
    }

    all[0] = (t_max[counter] + t_max[counter-1]) / 2.;  //Time needed to reach maxdist
    all[1] = (r_max[counter] + r_max[counter-1]) / 2.;  //Maxdist reached
    all[2] = (v_max[counter] + v_max[counter-1]) / 2.;  //Satellite velocity at maxdist
    all[3] = (v_0[counter] + v_0[counter-1]) / 2.;      //Initial velocity for satellite
    all[4] = time - all[0];                             //Starting time for satellite from earth

    return all;
}

int initialize_satellites(bool final, int counter, double v_min, double v_max, vector<double> x, vector<double> y, vector<double> z, vector<double> vx, vector<double> vy, vector<double> vz, vector<double> r, vector<double> &xs, vector<double> &ys, vector<double> &zs, vector<double> &vxs, vector<double> &vys, vector<double> &vzs, vector<double> &ms, vector<double> &rs, int sat, int startobject){
    //here sat satellites are initiallized given lowest and/or highest velocity and initial position
    double vso = sqrt(pow(vx[startobject],2) + pow(vy[startobject],2) + pow(vz[startobject],2));
    double rso = r[startobject];
    int prefactor = 0;

    //erase all previous values of satellites to only append current ones
    xs.erase(xs.begin(), xs.end());
    ys.erase(ys.begin(), ys.end());
    zs.erase(zs.begin(), zs.end());
    vxs.erase(vxs.begin(), vxs.end());
    vys.erase(vys.begin(), vys.end());
    vzs.erase(vzs.begin(), vzs.end());
    ms.erase(ms.begin(), ms.end());

    if (!final){
        if ( (sat != 100) && (sat != 1000) ){
            //start that only if sat!=100 AND sat != 1000
            if(v_min > DBL_EPSILON){
                //initialize satellites for v_min
                //satellite initial positions
                for(int i=0; i<sat; i++) xs.push_back(x[startobject] + vx[startobject] / vso * rso);
                for(int i=0; i<sat; i++) ys.push_back(y[startobject] + vy[startobject] / vso * rso);
                for(int i=0; i<sat; i++) zs.push_back(z[startobject] + vz[startobject] / vso * rso);

                //since velocity will not be an int
                for(int i=0; i<sat; i++) vxs.push_back((v_min - i*pow(10.,-counter)) * vx[startobject] / vso);
                for(int i=0; i<sat; i++) vys.push_back((v_min - i*pow(10.,-counter)) * vy[startobject] / vso);
                for(int i=0; i<sat; i++) vzs.push_back((v_min - i*pow(10.,-counter)) * vz[startobject] / vso);

                //satellite masses
                for(int i=0; i<sat; i++) ms.push_back(0.);
                prefactor++;
            }

            if(v_max > DBL_EPSILON){
                //initialize satellites for v_max
                //satellite initial positions
                for(int i=0; i<sat; i++) xs.push_back(x[startobject] + vx[startobject] / vso * rso);
                for(int i=0; i<sat; i++) ys.push_back(y[startobject] + vy[startobject] / vso * rso);
                for(int i=0; i<sat; i++) zs.push_back(z[startobject] + vz[startobject] / vso * rso);

                //satellite initial velocities
                for(int i=0; i<sat; i++) vxs.push_back((v_max + i*pow(10.,-counter)) / vso * vx[startobject]);
                for(int i=0; i<sat; i++) vys.push_back((v_max + i*pow(10.,-counter)) / vso * vy[startobject]);
                for(int i=0; i<sat; i++) vzs.push_back((v_max + i*pow(10.,-counter)) / vso * vz[startobject]);

                //satellite masses
                for(int i=0; i<sat; i++) ms.push_back(0.);
                prefactor++;
            }
        }
        else{
            //no velocity worked before, we know only v0 = v_min initially, which is assumed too high -> only test lower velocities
            for(int i=0; i<sat; i++) xs.push_back(x[startobject] + vx[startobject] / vso * rso);
            for(int i=0; i<sat; i++) ys.push_back(y[startobject] + vy[startobject] / vso * rso);
            for(int i=0; i<sat; i++) zs.push_back(z[startobject] + vz[startobject] / vso * rso);

            //satellite initial velocities
            for(double i=0; i<sat; i++) vxs.push_back((v_min - 10.*i/sat) * vx[startobject] / vso);
            for(double i=0; i<sat; i++) vys.push_back((v_min - 10.*i/sat) * vy[startobject] / vso);
            for(double i=0; i<sat; i++) vzs.push_back((v_min - 10.*i/sat) * vz[startobject] / vso);

            //satellite masses
            for(int i=0; i<sat; i++) ms.push_back(0.);
            prefactor = 1;
        }
    }
    else{
        double di = pow(10,-counter);
        for(double i=v_min; i<v_max+di;){
            //satellite initial positions
            xs.push_back(x[startobject] + vx[startobject] / vso * rso);
            ys.push_back(y[startobject] + vy[startobject] / vso * rso);
            zs.push_back(z[startobject] + vz[startobject] / vso * rso);

            //satellite initial velocities
            vxs.push_back(i * vx[startobject] / vso);
            vys.push_back(i * vy[startobject] / vso);
            vzs.push_back(i * vz[startobject] / vso);

            //satellite masses
            ms.push_back(0.);
            i += di;
            prefactor++;
        }
    }

    rs.resize(sat*prefactor);
    for(int i=0; i<sat*prefactor; i++) rs[i] = 0.;

    return prefactor;
}

void driver(double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, int n, vector<double> m, vector<double> &r, Step_function step, string command, double i){
    //Create and open output file
    fstream file;
    file.open(command+"-solution.csv", ios::out);
    file.precision(16);
    //double check = sqrt(pow(vx[10],2) + pow(vy[10],2) + pow(vz[10],2));
    //vx[10] = vx[10] * i / check;
    //vy[10] = vy[10] * i / check;
    //vz[10] = vz[10] * i / check;

    vector<double> x_old, y_old, z_old, vx_old, vy_old, vz_old;

    double Delta = 0.;
    double Delta_aim = 1e-16;
    double stepsize = 0;
    int count  = 0;
    int timestep = 100;
    //loop that iterates up to a certain chosen time (end)
    while((t_end - t) > DBL_EPSILON){
        if(count % timestep == 0){
          file << t << "; ";
              for(int i=0; i<n; i++) file << x[i] << "; ";
              for(int i=0; i<n; i++) file << y[i] << "; ";
              for(int i=0; i<n; i++) file << z[i] << "; ";
              for(int i=0; i<n; i++) file << vx[i] << "; ";
              for(int i=0; i<n; i++) file << vy[i] << "; ";
              for(int i=0; i<n-1; i++) file << vz[i] << "; ";
          file << vz[n-1]<< endl;
          count = 0;
        }

        x_old = x;
        y_old = y;
        z_old = z;
        vx_old = vx;
        vy_old = vy;
        vz_old = vz;

        //Calculate next timestep

        rk5_step(t, dt, x_old, y_old, z_old, vx_old, vy_old, vz_old, acceleration, n, m);

        Delta = 0;
        for(int i=0; i<n; i++) {
            Delta += x_old[i]+y_old[i]+z_old[i];
        }

        dt *= pow(Delta_aim/Delta, 1./5.);
        step(t, dt, x, y, z, vx, vy, vz, acceleration, n, m);

        //if the satellite collided with an object remove it from further calculations
        if ((n == 11) && (crash_check(t, x, y, z, r, x, y, z, n-1, n-1))){
            x.erase(x.begin()+(n-1));
            y.erase(y.begin()+(n-1));
            z.erase(z.begin()+(n-1));
            vx.erase(vx.begin()+(n-1));
            vy.erase(vy.begin()+(n-1));
            vz.erase(vz.begin()+(n-1));
            m.erase(m.begin()+(n-1));
            r.erase(r.begin()+(n-1));

            n -= 1;
            break;
        }

        t += dt;
        count ++; //Count  the number of time steps for saving only every 100th value.
    }
    file.close();
}

vector<double> check_for_boundaries(int precision, int n, double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> &m, vector<double> &r, vector<double> &xs, vector<double> &ys, vector<double> &zs, vector<double> &vxs, vector<double> &vys, vector<double> &vzs, vector<double> &ms, vector<double> &rs, double upper, double lower, string name, int startobject){
    //Know approximately (integer of) initial velocity (v0 ~ 8 < 10)
    //Idea send out some satellites and determine if they return inside a fixed interval of distance to the origin
    //Test 10 sattelites (one for each decimal point) wheather one returns in the interval
        //If no satellite returns test 100 satellites (one for each of the first two digits after the decimal point)
        //It at least one returns then check if we could go lower/higher with the initial velocity v0 by
        //sending out 10 satellites below the lowest digit that worked and ten above the highest digit that worked
    //Repeat sending out lower and higher satellites in even further digits behind the decimal point
    //If some given accuracy of v0 is reached, return highest and lowest velocity that let the satellite turn in the fixed intervall (v_max and v_min)
    //Those values are the "boundary values" for the velocity
    cout << "entered cfb" << endl;
    string input = "Orbits.csv";
    vector<string> tmp = {};
    vector<double> boundary = {0., 0.};

    int counter = 0; //keeps track of the digits after the decimal point
    int sat = 10;
    int prefactor;
    //double v0 = 10.;
    double v0 = 16.; //----just for testing
    double v_min = v0; //initial velocity and calculate below
    double v_max = 0.;
    bool out = false;
    bool final = false;

    vector<double> t_maxdist = {};
    vector<double> v_maxdist = {};
    vector<double> maxdist = {};
    vector<double> v0_sat = {};

    while(counter < (precision+1)){ //determines the precision of v0 in number of digits after the comma
        cout << "while loop " << counter << endl;
        //Old objects and satellites are destroyed and new ones created
        initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
        prefactor = initialize_satellites(final, counter, v_min, v_max, x, y, z, vx, vy, vz, r, xs, ys, zs, vxs, vys, vzs, ms, rs, sat, startobject);
        cout << "initialized satellites" << endl;

        //run programm to the end and get t_maxdist, maxdist, v_maxdist and v0_sat of satellites that returned in suitable interval
        sat_driver(counter, t, t_end, dt, x, y, z, vx, vy, vz, m, r, xs, ys, zs, vxs, vys, vzs, ms, rs, n, rk4_step, lower, upper, prefactor*sat, t_maxdist, maxdist, v_maxdist, v0_sat, out);

        //know from satellite driver how many satellites made it and know their initial velocities
        v_min = findmin(v0_sat);
        v_max = findmax(v0_sat);

        cout << "vmin = " << v_min << "; vmax = " << v_max << endl;

        if ((v_min == infinity()) || (v_max == -infinity())){
            sat *= 10;
            v_min = v0;
            v_max = 0.;
            counter++;
        }
        else{
            sat = 10;
            counter++;
        }
    }

    boundary[0] = v_min;
    boundary[1] = v_max;
    cout << "exit cfb" << endl;
    return boundary;
}

void calc_sat(vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> m, vector<double> &r, int n, Step_function step, string name){
    //Create and open output file
    string input = "Orbits.csv";
    bool out = true;
    bool final = true;
    double v_min, v_max;

    double t = 0.;
    double t_end = 20.;
    double dt = pow(2,-19);
    int startobject = 3;
    int endobject = 7; //no. planet -1; (Pluto = 8)
    int precision = 4;
    int satdummy = 10;
    int sat; //aka prefactor at another point

    vector<double> a(n), e(n), b(n), l(n), u(n);
    vector<vector<double>> values = {a,e,b,l,u};
    values = set_values(values, n-1, input);

    l = values[3];
    u = values[4];

    double lower = 99 * l[endobject]; //read 1% difference from max min orbit
    double upper = 101 * u[endobject];
    cout << "calculated upper lower" << endl;

    vector<double> t_maxdist = {};
    vector<double> maxdist = {};
    vector<double> v_maxdist = {};
    vector<double> v0_sat = {};

    vector<double> xs = {};
    vector<double> ys = {};
    vector<double> zs = {};
    vector<double> vxs = {};
    vector<double> vys = {};
    vector<double> vzs = {};
    vector<double> ms = {};
    vector<double> rs = {};

    vector<double> boundaries = check_for_boundaries(precision, n, t, t_end, dt, x, y, z, vx, vy, vz, m, r, xs, ys, zs, vxs, vys, vzs, ms, rs, upper, lower, name, startobject);
    v_min = boundaries[0];
    v_max = boundaries[1];

    cout << "boundaries[0] = " << boundaries[0] << "; boundaries[1] = " << boundaries[1] << endl;

    //initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
    //sat = initialize_satellites(final, precision, v_min, v_max, x, y, z, vx, vy, vz, r, xs, ys, zs, vxs, vys, vzs, ms, rs, satdummy, startobject);
    //sat_driver(t, t_end, dt, x, y, z, vx, vy, vz, m, r, xs, ys, zs, vxs, vys, vzs, ms, rs, n, rk4_step, lower, upper, sat, t_maxdist, maxdist, v_maxdist, v0_sat, out);
}

void programmteil(string command){
    //vectors chosen such that n-particles can be realized
    //Only current values are stored and after the output overridden with the new ones
    vector<double> x = {};
    vector<double> y = {};
    vector<double> z = {};
    vector<double> vx = {};
    vector<double> vy = {};
    vector<double> vz = {};
    vector<double> m = {};
    vector<double> r = {};

    int n = 10;                  //Number of objects
    double t_end = 10.;           //final time
    double dt = pow(2.,-24);     //time steps
    double t = 0.;

    string name = "Input.csv";

    if(fileexists(name)){
        if (command == "fwd"){  // forward euler
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
            driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, fwd_step, command ,t);
        }
        else if (command == "rk4"){ // Runge Kutta 4
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
            driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step, command, t);
        }
        else if (command == "lf"){ // leap frog
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
            driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, lf_step, command, t);
        }
        else if (command == "sat"){ // satellites
            calc_sat(x, y, z, vx, vy, vz, m, r, n, rk4_step, name);
        }
        else if (command == "test"){ //testing suite for sat velocity boundary applied with rk4 scheme
            n = 11;
            for(double i=9.18; i<9.19;){
                cout << "i = " << i << endl;
                t_end = 10.;
                dt = pow(2,-27); //~1E-6; 2-19 ~ 1e-7 ~ 3.76s -> 1e-6 ~ 37s
                initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
                driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step, command, i);
                i += 0.01;
            }
        }
        else cout << "Wrong parameter!" << endl;
    }else{
        cout << "Input file not found!" << endl;
    }
}

int main(int argc, char** argv){
    auto t1 = chrono::high_resolution_clock::now();
    //arguments of main used to hand over parameters at the programm start
    //1st parameter is programm name
    if(argc != 2) {
        cout << "usage:\n" << argv[0] << " <parameter> || possible parameters: 'fwd', 'rk4' or 'lf' " << endl;
        return -1;
    }else{
        stringstream input{argv[1]};
        string command;
        input >> command;

        programmteil(command);

        auto t2 = chrono::high_resolution_clock::now();
        auto time = chrono::duration<float>(t2-t1).count();
        cout << "Reached end of main after " << time << " seconds." << endl;

        return 0;
    }
}
