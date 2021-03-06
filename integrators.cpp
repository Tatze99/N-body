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
  double Delta_aim = 1e-14;
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

void driver(double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, int n, vector<double> m, vector<double> &r, Step_function step, string command){
    //Create and open output file
    fstream file;
    file.open(command+"-solution.csv", ios::out);
    file.precision(16);

    vector<double> x_old, y_old, z_old, vx_old, vy_old, vz_old;

    double Delta = 0.;
    double Delta_aim = 1e-14;
    int count  = 0;
    int timestep = 100;
    double dist_min = 100.;
    double dist;

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
          // cout << "dt = " << dt << endl;
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
        // if (int(t)%(int(t_end)/1000)==0){ cout << t/t_end << endl;}
        dist = distance(x[10], y[10], z[10], x[5], y[5], z[5]);

        if(dist < dist_min) {dist_min = dist;}

        t += dt;
        count ++; //Count  the number of time steps for saving only every 100th value.
    }
    cout << "Distance to Planet: " << dist_min << endl;

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
  double step_old = 0.;
  // distance between planet and earth
  double dist, dist2;
  double dist_old = 100, vsat_old = vsat; // gro?? w??hlen
  int mult = 1;

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

      step_old = step;
      cout << dist << " ; "<< dist2 << " ; " << step<< endl;
      if (dist > dist2) {
        step += stepsize;
        if (dist > 1) {
          step += dist*0.0001;
        }
        continue;
        }
      else {
        cout << "time = " << (time_est+step) << endl;
        cout << "vsat = " << vsat << endl;
        break;
      }
    }


    if (dist < 0.01) break; // break if distance is small enough 2*10^6 km
    if (dist_old < dist & dist < 3) break;

    if (step > stepsize ) step = step_old-stepsize;
    vsat_old = vsat;

    if (distance(x[10], y[10], z[10], 0,0,0)>distance(x[Planet], y[Planet], z[Planet], 0,0,0)){
      mult=-1;
    } else {mult = 1;}

    if (dist > 1){
      vsat += mult*0.005;  // larger stepsizes
    } else {
      vsat += mult*0.001; // if distance is small, reduce the increase of starting velocity
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

double optimize_vsat(double t, double t_end, double dt, vector<double> x, vector<double> y, vector<double> z, vector<double> vx, vector<double> vy, vector<double> vz, int n, vector<double> m, vector<double> r, int Planet, double &vsat, double Factorz){
  vector <double> x_old, y_old, z_old, vx_old, vy_old, vz_old;
  vector<double> x2, y2, z2, vx2, vy2, vz2;
  vector<double> x3, y3, z3, vx3, vy3, vz3;
  vector<double> dist_min, distxy_min, distz_min, velocity;
  double dist1, dist2, dist3, dist1_min, dist2_min, dist3_min, r1, r2, r3;

  dist1_min = DBL_MAX;
  dist2_min = DBL_MAX;
  dist3_min = DBL_MAX;
  x_old = x; y_old = y; z_old = z; vx_old = vx; vy_old = vy; vz_old = vz;

  x2 = x; y2 = y; z2 = z; vx2 = vx; vy2 = vy; vz2 = vz;
  x3 = x; y3 = y; z3 = z; vx3 = vx; vy3 = vy; vz3 = vz;

  double Factor1 = 1.;
  double Factor2 = 0.99;
  double Factor3 = 1.01;
  set_satellite_old(x,y,z,vx,vy,vz,vsat,r[3],Factorz,Factor1,n);
  set_satellite_old(x2,y2,z2,vx2,vy2,vz2,vsat,r[3],Factorz,Factor2,n);
  set_satellite_old(x3,y3,z3,vx3,vy3,vz3,vsat,r[3],Factorz,Factor3,n);

  while (t<t_end){
    // integrate towards two starting times
    integrator(t, t+dt/10, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step);
    integrator(t, t+dt/10, dt, x2, y2, z2, vx2, vy2, vz2, n, m, r, rk4_step);
    integrator(t, t+dt/10, dt, x3, y3, z3, vx3, vy3, vz3, n, m, r, rk4_step);

    dist1 = distance(x[10], y[10], z[10], x[Planet], y[Planet], z[Planet]);
    dist2 = distance(x2[10], y2[10], z2[10], x2[Planet], y2[Planet], z2[Planet]);
    dist3 = distance(x3[10], y3[10], z3[10], x3[Planet], y3[Planet], z3[Planet]);

    if (dist1 < dist1_min) {
      dist1_min=dist1;
      r1 = distance(x[10], y[10], 0, x[Planet], y[Planet], 0);
    }
    if (dist2 < dist2_min) {
      dist2_min=dist2;
      r2 = distance(x2[10], y2[10], 0, x2[Planet], y2[Planet], 0);
    }
    if (dist3 < dist3_min) {
      dist3_min=dist3;
      r3 = distance(x3[10], y3[10], 0, x3[Planet], y3[Planet], 0);
    }
    t+=dt;
  }
  double v1 = Factor1*vsat;
  double v2 = Factor2*vsat;
  double v3 = Factor3*vsat;

  //  https://www.arndt-bruenner.de/mathe/10/parabeldurchdreipunkte.htm
  double vsat_old = vsat;
  vsat = (v2*v2*(r3-r1)-v1*v1*(r3-r2)-v3*v3*(r2-r1))/(2*(v2*(r3-r1)-v1*(r3-r2)-v3*(r2-r1)));

  cout << "v1 =" << v1 << "; v2 =" << v2 << "; v3 =" << v3 << endl;
  cout << "r1 =" << r1 << "; r2 =" << r2 << "; r3 =" << r3 << endl;
  cout << "vsat =" << vsat << endl;
  fstream file;
  file.open("velocity_dependence.csv", ios::out);
  file.precision(20);
  file.setf(ios_base::fixed);

  double v = vsat + 0.1;
  while (v>(vsat-0.1)){
    vector <double> dist, distxy, distz;
    x = x_old; y = y_old; z = z_old; vx = vx_old; vy = vy_old; vz = vz_old;
    set_satellite_old(x,y,z,vx,vy,vz,vsat,r[3],1, v/vsat,n);
    t = 0;

    while (t < t_end){
      integrator(t, t+dt/10, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step);
      dist.push_back(distance(x[10], y[10], z[10], x[Planet], y[Planet], z[Planet]));
      distxy.push_back(distance(x[10], y[10], 0, x[Planet], y[Planet], 0));
      distz.push_back(distance(0, 0, z[10], 0, 0, z[Planet]));
      t+=dt;
    }
    dist_min.push_back(findmin(dist));
    distxy_min.push_back(findmin(distxy));
    distz_min.push_back(findmin(distz));
    velocity.push_back(v);
    int arg_dist_min = findargmin(dist);

    file << v << ";" << findmin(dist) << ";" << distxy[arg_dist_min] << ";" << distz[arg_dist_min] << endl;
    cout << "v = " << v << ";  rmin = " << findmin(dist) <<endl;
    if (findmin(dist) < 0.04){
      v -= 0.001;
    } else{
      v -= 0.005;
    }

  }
  int arg_vmin = findargmin(dist_min);
  cout << "velocity min = " << velocity[arg_vmin] << endl;
  return velocity[arg_vmin]/vsat;
}

double optimize_initial_angle(double t, double t_end, double dt, vector<double> x, vector<double> y, vector<double> z, vector<double> vx, vector<double> vy, vector<double> vz, int n, vector<double> m, vector<double> r, int Planet, double vsat, double Factorxy){
  vector<double> x2, y2, z2, vx2, vy2, vz2;
  double dist, dist2, dist_min, dist2_min, distz_min, distz2_min;
  dist_min = 100.;
  dist2_min =100.;

  x2 = x; y2 = y; z2 = z; vx2 = vx; vy2 = vy; vz2 = vz;

  double Factor1 = 1.;
  double Factor2 = 0.95;
  set_satellite_old(x,y,z,vx,vy,vz,vsat,r[3],Factor1,Factorxy,n);
  set_satellite_old(x2,y2,z2,vx2,vy2,vz2,vsat,r[3],Factor2,Factorxy,n);

  while (t<t_end){
    // integrate towards two starting times
    integrator(t, t+dt/10, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step);
    integrator(t, t+dt/10, dt, x2, y2, z2, vx2, vy2, vz2, n, m, r, rk4_step);
    dist = distance(x[10], y[10], z[10], x[Planet], y[Planet], z[Planet]);
    dist2 = distance(x2[10], y2[10], z2[10], x2[Planet], y2[Planet], z2[Planet]);

    if (dist < dist_min) {
      dist_min=dist;
      distz_min = distance(0, 0, z[10], 0, 0, z[Planet]);
    }
    if (dist2 < dist2_min) {
      dist2_min=dist2;
      distz2_min = distance(0, 0, z2[10], 0, 0, z2[Planet]);
    }
    t+=dt;
  }
  double Anstieg, Achsenabschnitt;
  cout << "dist 1 = " << dist_min << endl;
  cout << "dist 2 = " << dist2_min << endl;
  Anstieg = (distz_min-distz2_min)/(Factor1-Factor2);
  Achsenabschnitt = distz_min - Anstieg;
  cout << "Anstieg " << Anstieg << endl;
  cout << "Achsenabschnitt" << Achsenabschnitt << endl;

  return -Achsenabschnitt/Anstieg;
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
    // double t_end = 0.462693;           //final time  9.09745;
    int Planet = 3;  // 0 = Merkur, 2 = Erde
    // double t_end = 36.5573329;   // New Horizon
    // double t_end = 36.5648329;   // New Horizon test
    double t_end = 100;
    double dt = pow(2.,-25);     //time steps
    double t = 0.;

    string name = "Input.csv";

    if(fileexists(name)){
        if (command == "fwd"){  // forward euler
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
            driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, fwd_step, command);
        }
        else if (command == "rk4"){ // Runge Kutta 4
            n=11;
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
            driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step, command);

            // double vsat = 9.53855;
            // double vsat = pow(vx[3]*vx[3]+vy[3]*vy[3]+vz[3]*vz[3], 0.5) + 3.4172; //  New Horizon
            double vsat = pow(vx[3]*vx[3]+vy[3]*vy[3]+vz[3]*vz[3], 0.5) + 3.449; //  New Horizon korrigiert
            set_satellite_old(x,y,z,vx,vy,vz,vsat,r[3],1,1,n);
            fstream file;
            // file.open("Input_tend"+to_string(t_end)+".csv", ios::out);
            file.open("Input_tend.csv", ios::out);
            file.precision(20);
            file.setf(ios_base::fixed);

            for(int i=0; i<n; i++) file << x[i] << ";" << y[i] << ";" << z[i] << ";" << vx[i] << ";" << vy[i] << ";" << vz[i] << ";" << m[i] << ";" << r[i] << endl;
        }
        else if (command == "lf"){ // leap frog
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
            driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, lf_step, command);
        }
        else if (command == "sat"){ // satellites
            calc_sat(x, y, z, vx, vy, vz, m, r, 10, rk4_step, name);
        }
        else if (command == "test"){ //testing suite for sat velocity boundary applied with rk4 scheme
            n = 11;
            for(double i=9.2189; i<9.22;){
                cout << "i = " << i << endl;
                t_end = 10;
                dt = pow(2,-27); //~1E-6; 2-19 ~ 1e-7 ~ 3.76s -> 1e-6 ~ 37s
                double check = sqrt(pow(vx[10],2) + pow(vy[10],2) + pow(vz[10],2));
                vx[10] *= i / check;
                vy[10] *= i / check;
                vz[10] *= i / check;
                initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
                driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step, command);
                i += 0.01;
            }
        }
        else if (command == "angle"){
            int n = 11;
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, "Input.csv");
            vector<double> a(n), e(n), b(n), l(n), u(n), Umlaufzeit(n), Hohmann(n), vmin(n);
            vector<vector<double>> values = {a,e,b,l,u, Umlaufzeit, Hohmann, vmin};
            values = set_values(values, 9, "Orbits.csv");

            Umlaufzeit = values[5];
            Hohmann = values[6];
            vmin = values[7];

            double Umlauf = Umlaufzeit[Planet];
            double time_orbit = Hohmann[Planet];
            double vsat = vmin[Planet];
            cout << Umlauf << ";" << time_orbit << ";" << vsat << endl;

            calc_angle(t, Umlauf, dt, x, y, z, vx, vy, vz, 11, m, Planet+1, time_orbit, Umlauf, r, vsat);
        }
        else if (command == "calc_t"){
            int n = 11;
            vector<double> a(n), e(n), b(n), l(n), u(n), Umlaufzeit(n), Hohmann(n), vmin(n);
            vector<vector<double>> values = {a,e,b,l,u, Umlaufzeit, Hohmann, vmin};
            values = set_values(values, 9, "Orbits.csv");
            Hohmann = values[6];
            double t_end = Hohmann[Planet];
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, "Input_tend.csv");
            driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step, command);
        }
        else if (command == "Swing"){
          n=11;
          t_end = 1.5;
          initialize_objects(n, x, y, z, vx, vy, vz, m, r, "Input_Horizon.csv");
          double vsat = pow(vx[3]*vx[3]+vy[3]*vy[3]+vz[3]*vz[3], 0.5) + 3.449; //  New Horizon
          double Faktorxy = optimize_vsat(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, 5, vsat, 1);
          double Faktorz = optimize_initial_angle(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, 5, vsat, Faktorxy);
          cout << "Faktor z = " << Faktorz <<  endl;
          cout << "Faktor xy = " << Faktorxy <<  endl;
          cout << "vsat = " << vsat << endl;

          set_satellite_old(x,y,z,vx,vy,vz,vsat,r[3],Faktorz,Faktorxy, n);

          double t_end = 9.5;
          driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step, command);
        }

        else if (command == "Swing-by"){
          int Planet = 9;
          n=11;
          vector<double> a(n), e(n), b(n), l(n), u(n), Umlaufzeit(n), Hohmann(n), vmin(n);
          vector<vector<double>> values = {a,e,b,l,u, Umlaufzeit, Hohmann, vmin};
          values = set_values(values, 9, "Orbits.csv");
          Hohmann = values[6];
          double t_end = Hohmann[Planet-1];
          double vsat;

          if (Planet == 4){
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, "Input_tend_Mars.csv");
            vsat = 8.621; //Mars

          } else if (Planet == 5){
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, "Input_tend_Jupiter.csv");
            vsat = 9.35;  // Jupiter

          } else if (Planet == 6){
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, "Input_tend_Saturn.csv");
            vsat = 9.389;  // Saturn

          } else if (Planet == 7){
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, "Input_tend_Uranus.csv");
            vsat = 9.71176;  // Uranus

          } else if (Planet == 8){
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, "Input_tend_Neptun.csv");
            vsat = 9.78886;  // Uranus
          } else if (Planet == 9){
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, "Input_tend_Pluto.csv");
            vsat = 9.78886;  // Uranus
          }
          // const double v0 = 6.179554659505929;
          // double ve = pow(vx[3]*vx[3]+vy[3]*vy[3]+vz[3]*vz[3], 0.5);
          // vsat += ve - v0;

          double Faktorxy = optimize_vsat(t, t_end+0.1, dt, x, y, z, vx, vy, vz, n, m, r, Planet, vsat, 1);
          double Faktorz = 1; //optimize_initial_angle(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, Planet, vsat, Faktorxy);
          cout << "Faktor z = " << Faktorz <<  endl;
          cout << "Faktor xy = " << Faktorxy <<  endl;
          cout << "vsat = " << vsat << endl;

          set_satellite_old(x,y,z,vx,vy,vz,vsat,r[3],Faktorz,Faktorxy, n);
          // driver(t, t_end*6, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step, command);
        }

        else if (command == "NewHorizon"){
          n=11;
          initialize_objects(n, x, y, z, vx, vy, vz, m, r, "Input_Horizon.csv");
          n=13;
          double vsat= pow(vx[3]*vx[3]+vy[3]*vy[3]+vz[3]*vz[3], 0.5) + 3.449;
          x.resize(n); y.resize(n); z.resize(n); vx.resize(n); vy.resize(n); vz.resize(n); m.resize(n); r.resize(n);
          double Faktorz = 0.911673; // first z, then xy
          double Faktorxy = 1.00083;
          set_satellite_old(x,y,z,vx,vy,vz,vsat,r[3],Faktorz,Faktorxy, 12);
          Faktorxy = 1.00123; // first xy, then z
          Faktorz = 0.914429;
          set_satellite_old(x,y,z,vx,vy,vz,vsat,r[3],Faktorz,Faktorxy, 13);
          double t_end = 9.48;
          driver(t, t_end, dt, x, y, z, vx, vy, vz, 13, m, r, rk4_step, command);
        }

        else if (command == "test2"){
          n=11;
          initialize_objects(n, x, y, z, vx, vy, vz, m, r, "Input_Horizon.csv");
          double vsat = pow(vx[3]*vx[3]+vy[3]*vy[3]+vz[3]*vz[3], 0.5) + 3.4172; //  New Horizon
          set_satellite_old(x,y,z,vx,vy,vz,vsat,r[3],1,1,n);
          cout << vsat << endl;
          double t_end = 1.5;
          double Faktor = optimize_vsat(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, 5, vsat, 1);
          cout << Faktor << endl;
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
