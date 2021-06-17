/*  Dies ist ein Programm zur Implementierung des forward-euler,
    rk4 und des leapfrog-Verfahrens.
    C://Users/Martin/Documents/GitHub/N-body
    g++ integrators.cpp -o integrators
    ./integrators rk4
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

int initialize_satellites(int trackinit, bool final, int counter, double v_min, double v_max, vector<double> x, vector<double> y, vector<double> z, vector<double> vx, vector<double> vy, vector<double> vz, vector<double> r, vector<double> &xs, vector<double> &ys, vector<double> &zs, vector<double> &vxs, vector<double> &vys, vector<double> &vzs, vector<double> &ms, vector<double> &rs, int sat, int startobject){
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

        if (sat == 100){
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

        if (sat == 1000){
            sat = 100;
            //no velocity worked before, we know only v0 = v_min initially, which is assumed too high -> only test lower velocities
            for(int i=0; i<sat; i++) xs.push_back(x[startobject] + vx[startobject] / vso * rso);
            for(int i=0; i<sat; i++) ys.push_back(y[startobject] + vy[startobject] / vso * rso);
            for(int i=0; i<sat; i++) zs.push_back(z[startobject] + vz[startobject] / vso * rso);

            //satellite initial velocities
            for(double i=0; i<sat; i++) vxs.push_back((vso + trackinit + i/sat) * vx[startobject] / vso);
            for(double i=0; i<sat; i++) vys.push_back((vso + trackinit + i/sat) * vy[startobject] / vso);
            for(double i=0; i<sat; i++) vzs.push_back((vso + trackinit + i/sat) * vz[startobject] / vso);

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

    if (sat != 1000){
        rs.resize(sat*prefactor);
        for(int i=0; i<sat*prefactor; i++) rs[i] = 0.;
    }
    else{
        rs.resize(100);
        for(int i=0; i<100; i++) rs[i] = 0.;
    }
    return prefactor;

}

void integrator(double t, double t_end, double &dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, int n, vector<double> m, vector<double> &r, Step_function step){
  vector<double> x_old, y_old, z_old, vx_old, vy_old, vz_old;

  double Delta = 0.;
  double Delta_aim = 1e-16;
  //loop that iterates up to a certain chosen time (end)
  while((t_end - t) > DBL_EPSILON){
      x_old = x;
      y_old = y;
      z_old = z;
      vx_old = vx;
      vy_old = vy;
      vz_old = vz;

      //Calculate next timestep
      rk5_step(t, dt, x_old, y_old, z_old, vx_old, vy_old, vz_old, acceleration, n, m);

      Delta = 0;
      for(int i=0; i<n; i++) Delta += x_old[i]+y_old[i]+z_old[i];

      dt *= pow(Delta_aim/Delta, 1./5.);
      step(t, dt, x, y, z, vx, vy, vz, acceleration, n, m);
      t += dt;
  }
}

void driver(double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, int n, vector<double> m, vector<double> &r, Step_function step, string command, double i){
    //Create and open output file
    fstream file;
    file.open(command+"-solution.csv", ios::out);
    file.precision(16);
    // double check = sqrt(pow(vx[10],2) + pow(vy[10],2) + pow(vz[10],2));
    // vx[10] = vx[10] * i / check;
    // vy[10] = vy[10] * i / check;
    // vz[10] = vz[10] * i / check;

    vector<double> x_old, y_old, z_old, vx_old, vy_old, vz_old;

    double Delta = 0.;
    double Delta_aim = 1e-16;
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

            n --;
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
    //int trackinit = -1; // -1 to start at planet velocity
    int trackinit = 2; //startint at planet velocity+2 -- just speeding things up we start at 9 (2)
    int prefactor;
    double v0 = 10.;
    //double v0 = 16.; //----just for testing
    double v_min = v0; //initial velocity and calculate below
    double v_max = 0.;
    bool out = false;
    bool final = false;

    vector<double> t_maxdist, v_maxdist, maxdist, v0_sat;

    while(counter < (precision+1)){ //determines the precision of v0 in number of digits after the comma
        cout << "while loop " << counter << endl;
        //Old objects and satellites are destroyed and new ones created
        initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
        prefactor = initialize_satellites(trackinit, final, counter, v_min, v_max, x, y, z, vx, vy, vz, r, xs, ys, zs, vxs, vys, vzs, ms, rs, sat, startobject);
        if (sat == 1000) sat = 100;
        cout << "initialized satellites" << endl;

        //run programm to the end and get t_maxdist, maxdist, v_maxdist and v0_sat of satellites that returned in suitable interval
        sat_driver(t, t_end, dt, x, y, z, vx, vy, vz, m, r, xs, ys, zs, vxs, vys, vzs, ms, rs, n, rk4_step, lower, upper, prefactor*sat, t_maxdist, maxdist, v_maxdist, v0_sat, out);

        //know from satellite driver how many satellites made it and know their initial velocities
        v_min = findmin(v0_sat);
        v_max = findmax(v0_sat);

        cout << "vmin = " << v_min << "; vmax = " << v_max << endl;

        if (((DBL_MAX-v_min) < DBL_EPSILON) || (v_max == -M_PI)){
            sat *= 10;
            v_min = v0;
            v_max = 0.;
            if (sat == 100) counter++;
            if (sat == 1000) trackinit++;
            //if (trackinit == 0) counter++;
            if (trackinit == 3) counter++; //just for speedingup purposes and starting at 9
            if (trackinit > 5){
                cout << "Broke while" << endl;
                break; //break in order not to go to infinity
            }
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
    double t_end = 5.;
    double dt = pow(2,-22);
    int startobject = 3;
    int endobject = 4; //no. planet -1; (Pluto = 8)
    int precision = 4;
    int sat; //aka prefactor at another point

    vector<double> a(n), e(n), b(n), l(n), u(n);
    vector<vector<double>> values = {a,e,b,l,u};
    cout << "Test" << endl;
    values = set_values(values, n-1, input);
    l = values[3];
    u = values[4];

    double lower = 99 * l[endobject]; //read 1% difference from max min orbit
    double upper = 101 * u[endobject];
    cout << "calculated upper lower" << endl;

    vector<double> t_maxdist, maxdist, v_maxdist, v0_sat;
    vector<double> xs, ys, zs, vxs, vys, vzs, ms, rs;

    vector<double> boundaries = check_for_boundaries(precision, n, t, t_end, dt, x, y, z, vx, vy, vz, m, r, xs, ys, zs, vxs, vys, vzs, ms, rs, upper, lower, name, startobject);
    v_min = boundaries[0];
    v_max = boundaries[1];

    cout << "boundaries[0] = " << boundaries[0] << "; boundaries[1] = " << boundaries[1] << endl;

    //initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
    //sat = initialize_satellites(0, final, precision, v_min, v_max, x, y, z, vx, vy, vz, r, xs, ys, zs, vxs, vys, vzs, ms, rs, 12345, startobject);
    //sat_driver(t, t_end, dt, x, y, z, vx, vy, vz, m, r, xs, ys, zs, vxs, vys, vzs, ms, rs, n, rk4_step, lower, upper, sat, t_maxdist, maxdist, v_maxdist, v0_sat, out);
}

void calc_angle(double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, int n, vector<double> m, int Planet, double time_orbit, double Umlauf, vector<double> &r, double vsat){

  double Delta0 = 0, Delta1 = 0.1;
  vector<double> x_old, y_old, z_old, vx_old, vy_old, vz_old; // save current position
  vector<double> x2, y2, z2, vx2, vy2, vz2;

  x_old = x; y_old = y; z_old = z; vx_old = vx; vy_old = vy; vz_old = vz;

  double omega = 2*M_PI/Umlauf; // angular velocity of the planet
  // cout << "Delta aim = " << M_PI - omega*time_orbit << endl;

  while((t_end - t) > DBL_EPSILON){

    // calculate angle between the Planet and the earth
    Delta0 = angle(x[3], y[3], z[3], x[Planet], y[Planet], z[Planet]);
    integrator(t, t+dt/2, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step); // ensure only one step
    Delta1 = angle(x[3], y[3], z[3], x[Planet], y[Planet], z[Planet]);

    // angle should be Delta = phi(t=torbit, Planet)
    if (Delta0 > Delta1){
      if ((Delta1 < (M_PI - omega*time_orbit)) && (Delta1 > (M_PI - omega*time_orbit-0.1))){
        cout << "Planet: " << x[Planet] <<"; "<< y[Planet] <<"; "<< z[Planet] << endl;
        cout << "Erde: " << x[3] <<"; "<< y[3] <<"; "<< z[3] << endl;
        cout << "end time = " << t << ";  Angle = " << Delta1 << endl;
        break;
      }
    }
    t += dt;
  }

  double time_before = 0.1;
  double time_est = t-time_before;
  // integrate to a time slightly before the estimate
  integrator(0, t-time_before, dt, x_old, y_old, z_old, vx_old, vy_old, vz_old, n, m, r, rk4_step);

  t = 0; // set t back to zero
  // Set the position of the probe to earths position + radius and set the velocity
  set_satellite(x_old,y_old,z_old,vx_old,vy_old,vz_old,vsat,r[3]);

  // stepsize for integration time
  double stepsize = 0.005;
  double step = 0.;
  // distance between planet and earth
  double dist, dist2;
  double dist_old = 100, vsat_old = vsat; // groß wählen

  while (vsat < 10){
    for(int i=0; i<100; i++){ // arbitrary upper boundary, probably never reached

      x = x_old; y = y_old; z = z_old; vx = vx_old; vy = vy_old; vz = vz_old;
      x2 = x_old; y2 = y_old; z2 = z_old; vx2 = vx_old; vy2 = vy_old; vz2 = vz_old;

      // integrate towards two starting times
      integrator(t, t+step, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step);
      integrator(t, t+step+stepsize, dt, x2, y2, z2, vx2, vy2, vz2, n, m, r, rk4_step);

      // set the probe positions for the two times
      set_satellite(x,y,z,vx,vy,vz,vsat,r[3]);
      set_satellite(x2,y2,z2,vx2,vy2,vz2,vsat,r[3]);
      // integrate for the time of Hohmann transfer orbit
      integrator(t+step, t+step+time_orbit, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step);
      integrator(t+step+stepsize, t+step+stepsize+time_orbit, dt, x2, y2, z2, vx2, vy2, vz2, n, m, r, rk4_step);

      // calculate the distances
      dist = distance(x[10], y[10], z[10], x[Planet], y[Planet], z[Planet]);
      dist2 = distance(x2[10], y2[10], z2[10], x2[Planet], y2[Planet], z2[Planet]);

      cout << dist << " ; "<< dist2 << " ; " << step<< endl;
      if (dist > dist2) {
        step += stepsize;

        continue;
        }
      else {
        cout << "time = " << (time_est+step) << endl;
        cout << "vsat = " << vsat << endl;
        break;
      }
    }

    if (dist < 0.01) break; // break if distance is small enough 2*10^6 km
    if (dist_old < dist) {
      vsat -= 0.005; // break if distance gets larger again
      break;
    }

    if (step > stepsize ) step -= stepsize;
    vsat_old = vsat;
    if (dist > 0.2){
      vsat += 0.01;  // larger stepsizes
    } else {
      vsat += 0.001; // if distance is small, reduce the increase of starting velocity
    }

    dist_old = dist;
  }

  cout << "vsat = " << vsat << endl;
  x = x_old; y = y_old; z = z_old; vx = vx_old; vy = vy_old; vz = vz_old;
  integrator(t, t+step, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step);

  set_satellite(x,y,z,vx,vy,vz,vsat,r[3]);
  fstream file;
  file.open("Input_tend.csv", ios::out);
  file.precision(20);
  file.setf(ios_base::fixed);

  for(int i=0; i<n; i++) file << x[i] << ";" << y[i] << ";" << z[i] << ";" << vx[i] << ";" << vy[i] << ";" << vz[i] << ";" << m[i] << ";" << r[i] << endl;
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

    int n = 11;                  //Number of objects
    double t_end = 3.54289877;           //final time  9.09745;
    // double t_end = 5.81159;
    double dt = pow(2.,-24);     //time steps
    double t = 0.;

    string name = "Input_tend.csv";

    if(fileexists(name)){
        if (command == "fwd"){  // forward euler
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
            driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, fwd_step, command, t);
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
            calc_sat(x, y, z, vx, vy, vz, m, r, 10, rk4_step, name);
        }
        else if (command == "test"){ //testing suite for sat velocity boundary applied with rk4 scheme
            n = 11;
            for(double i=9.1706; i<9.18;){
                cout << "i = " << i << endl;
                t_end = 35.;
                dt = pow(2,-27); //~1E-6; 2-19 ~ 1e-7 ~ 3.76s -> 1e-6 ~ 37s
                initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
                driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step, command, i);
                i += 0.01;
            }
        }
        else if (command == "angle"){
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, "Input.csv");
            double time_orbit = 2.54289877;
            double vsat = 9.1706;
            double Umlauf = 11.862;
            calc_angle(t, 12., dt, x, y, z, vx, vy, vz, 11, m, 5, time_orbit, Umlauf, r, vsat);
        }
        else if (command == "calc_t"){
            double t_end = 2.54289877;
            // double t_end = 5;
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, "Input_tend.csv");
            driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step, command, t);
        }

        /*else if (command == "time"){ //testing suite for sat velocity boundary applied with rk4 scheme
            n = 11;

            //hier objdist berechnen
            //hier time berechnen
            vector<double> all = all_from_target(name, );
            initialize_objects(n-1, x, y, z, vx, vy, vz, m, r, name);
            driver(t, all[4], dt, x, y, z, vx, vy, vz, n, m, r, rk4_step, command, i);


            n = 11;
            for(double i=9.55625; i<9.557;){
                cout << "i = " << i << endl;
                t_end = 35.;
                dt = pow(2,-27); //~1E-6; 2-19 ~ 1e-7 ~ 3.76s -> 1e-6 ~ 37s
                initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
                driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step, command, i);
                i += 0.01;
            }
        }*/
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
