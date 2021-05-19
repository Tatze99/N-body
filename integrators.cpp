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

//for global "short-hand" notation - need not to write 'std::' in front of most things
using namespace std;

typedef void (* Step_function)(double, double, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, tuple<vector<double>,vector<double>,vector<double>>(double, vector<double>, vector<double>, vector<double>, int, vector<double>), int, vector<double>);

typedef void (* Step_satellite)(double, double, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, tuple<vector<double>,vector<double>,vector<double>>(double, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, int, int, vector<double>), int, vector<double>, int);

typedef tuple<vector<double>,vector<double>,vector<double>> (DGL)(double t, vector<double> x, vector<double> y, vector<double> z, int n, vector<double> m);

typedef tuple<vector<double>,vector<double>,vector<double>> (DGL_sat)(double t, vector<double> xs, vector<double> ys, vector<double> zs, vector<double> x, vector<double> y, vector<double> z, int n, int sat, vector<double> m);

tuple<vector<double>,vector<double>,vector<double>> acceleration(double t, vector<double> x, vector<double> y, vector<double> z, int n, vector<double> m){
    //cout << "execute acc" << endl;
  double Matrix[n][n];

  vector<double> ax(n);
  vector<double> ay(n);
  vector<double> az(n);
  for(int i=0; i<n; i++) ax.push_back(0.);
  for(int i=0; i<n; i++) ay.push_back(0.);
  for(int i=0; i<n; i++) az.push_back(0.);

  for(int i=0; i<n; i++){
    for(int j=0; j<i; j++) {
      Matrix[i][j] = pow((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]),-1.5);
    }
  }

  for(int i=0; i<n; i++){
    Matrix[i][i] = 0.;
    for(int j=i+1; j<n; j++) Matrix[i][j] = Matrix[j][i];
    for(int j=0; j<n; j++) {
      double c = Matrix[i][j]*m[j];
      ax[i] += (x[i]-x[j])*c;
      ay[i] += (y[i]-y[j])*c;
      az[i] += (z[i]-z[j])*c;
    }
    ax[i] *= -4*M_PI*M_PI;
    ay[i] *= -4*M_PI*M_PI;
    az[i] *= -4*M_PI*M_PI;
  }
  return make_tuple(ax, ay, az);
}

tuple<vector<double>,vector<double>,vector<double>> sat_acceleration(double t, vector<double> xs, vector<double> ys, vector<double> zs, vector<double> x, vector<double> y, vector<double> z, int n, int sat, vector<double> m){
  //cout << "execute sat_acc" << endl;
  double Matrix[sat][n];

  vector<double> ax(sat);
  vector<double> ay(sat);
  vector<double> az(sat);
  for(int i=0; i<sat; i++) ax.push_back(0.);
  for(int i=0; i<sat; i++) ay.push_back(0.);
  for(int i=0; i<sat; i++) az.push_back(0.);

  for(int i=0; i<sat; i++){
    for(int j=0; j<n; j++) {
      Matrix[i][j] = pow((xs[i]-x[j])*(xs[i]-x[j])+(ys[i]-y[j])*(ys[i]-y[j])+(zs[i]-z[j])*(zs[i]-z[j]),-1.5);
    }
  }

  for(int i=0; i<sat; i++){
    for(int j=0; j<n; j++){
        double c = Matrix[i][j]*m[j];
        ax[i] += (x[i]-x[j])*c;
        ay[i] += (y[i]-y[j])*c;
        az[i] += (z[i]-z[j])*c;
    }
    ax[i] *= -4*M_PI*M_PI;
    ay[i] *= -4*M_PI*M_PI;
    az[i] *= -4*M_PI*M_PI;
  }
  //cout << "Matrix[2][0-5] = " << Matrix[2][0] << "; " << Matrix[2][1] << "; " << Matrix[2][2] << "; " << Matrix[2][3] << Matrix[2][4] << "; " << Matrix[2][5] << "; " << endl;
  //value of venus and earth too large

  return make_tuple(ax, ay, az);
}

/*
vector<double> cosine(double t, vector<double> x, vector<double> y, vector<double> z, int n, vector<double> m){
  vector<double> u_rhs(n,0.0);
  for(int i=0; i<n; i++) u_rhs[i] = cos(t);
}

vector<double> testsymplectic(double t, vector<double> x, vector<double> y, vector<double> z, int n, vector<double> m){
  vector<double> u_rhs(n,0.0);
  for(int i=0; i<n; i++) u_rhs[i] = -sin(t);
}
*/

bool crash_check(vector<double> x, vector<double> y, vector<double> z, vector<double> r, vector<double> xs, vector<double> ys, vector<double> zs, int n, int satno){
    cout << "run crach_check" << endl;
    vector<double> dist = {};
    bool counter = false;
    vector<string> objects = {"Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"};

    //if distance of sattelite number satno to the i-th object is smaller than its radius return true
    for(int i=0; i<n; i++) dist.push_back(sqrt(pow(x[i]-xs[satno],2) + pow(y[i]-ys[satno],2) + pow(z[i]-zs[satno],2)));
    for(int i=0; i<n; i++){
        if ((dist[i] - r[i]) < DBL_EPSILON){
            cout << "Satellite No. " << satno << " crashed with " << objects[i] << endl;
            counter = true;
        }
    }
    return counter;
}

double findmin(vector<double> v){
    double min = infinity();
    for(double i : v) if ((min-i) > DBL_EPSILON) min = i;
    return min;
}

double findmax(vector<double> v){
    double max = -infinity();
    for(double i : v) if ((i-max) > DBL_EPSILON) max = i;
    return max;
}

inline bool fileexists(const string& name){
    if (FILE* file = fopen(name.c_str(), "r")){
        fclose(file);
        free(file);
        return true;
    }else{
        return false;
    }
}

void read_file(vector<string> &v, string &filename){
    //read complete file and store to auxiliary vector v
    fstream s;
    char cstring[256];
    string tmp;
    string str = "#";

    //open file and read data
    s.open(filename, ios::in);

    while (!s.eof()){
        s.getline(cstring, sizeof(cstring));
        tmp = cstring;
        //comment in input files with "#" as first character
        if (tmp[0] != str[0]) v.push_back(tmp);
    }

    s.close();
}

vector<double> upperlower(vector<string> &tmp, string input, int endobject){
    string s, tmp_s;
    string semi = ";";
    int counter;
    vector<double> ul = {0., 0.};

    if (fileexists(input)) read_file(tmp, input);
    s = tmp[endobject];
    //iterate over each double in one line
    for(int i=0; i<5; i++){
        //iterate over each line
        for(int j=0; j<s.length(); j++){
            if(s[j] == semi[0]){
                counter = j;
                break;
            }
        }

        //cut off the first part of s
        tmp_s = s.substr(0, counter);

        //set s to the remaining string
        if ((counter-2) < (s.length()-1)) s = s.substr(counter+2, s.length()-1);

        //set values
        switch (i)
        {
            case 3: ul[0] = stod(tmp_s);
                    break;
            case 4: ul[1] = stod(tmp_s);
        }
    }
    return ul;
}

void set_startvalues(int n, vector<string> help, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> &m, vector<double> &r){
   int counter;
   string s, tmp_s;
   string semi = ";";  //define delimiter which shall be searched for

   //iterate over auxiliary vector
    for(int i=0; i<n; i++){
        s = help[i];
        //iterate over each double in one line
        for(int l=0; l<8; l++){
            //iterate over each line --- mabye replace by find_first_of (see cppreference)
            for(int j=0; j<s.length(); j++){
                if(s[j] == semi[0]){
                    counter = j;
                    break;
                }
            }

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
                case 3: vx[i] = 365.245*stod(tmp_s);
                        break;
                case 4: vy[i] = 365.245*stod(tmp_s);
                        break;
                case 5: vz[i] = 365.245*stod(tmp_s);
                        break;
                case 6: m[i] = stod(tmp_s);
                        break;
                case 7: r[i] = stod(tmp_s);
            }
        }
    }
}

int initialize_satellites(bool final, int counter, double v_min, double v_max, vector<double> x, vector<double> y, vector<double> z, vector<double> vx, vector<double> vy, vector<double> vz, vector<double> r, vector<double> &xs, vector<double> &ys, vector<double> &zs, vector<double> &vxs, vector<double> &vys, vector<double> &vzs, vector<double> &ms, int sat, int startobject){
    //here sat satellites are initiallized given lowest and/or highest velocity and initial position
    double vso = sqrt(pow(vx[startobject],2) + pow(vy[startobject],2) + pow(vz[startobject],2));
    double rso = 4.2644 * pow(10,-5);
    //double rso = r[startobject];
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

                //since velocity will nor be an int
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
    return prefactor;
}

void initialize_objects(int n, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> &m, vector<double> &r, string name){
    //auxiliary vector
    vector<string> help = {};

    x.erase(x.begin(), x.end());
    y.erase(y.begin(), y.end());
    z.erase(z.begin(), z.end());
    vx.erase(vx.begin(), vx.end());
    vy.erase(vy.begin(), vy.end());
    vz.erase(vz.begin(), vz.end());
    m.erase(m.begin(), m.end());

    //Change size of vectors
    x.resize(n);
    y.resize(n);
    z.resize(n);
    vx.resize(n);
    vy.resize(n);
    vz.resize(n);
    m.resize(n);
    r.resize(n);

    read_file(help, name);
    set_startvalues(n, help, x, y, z, vx, vy, vz, m, r);
}

void fwd_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, DGL rhs, int n, vector<double> m){
    //Initialize vectors for the steps - only one step here!
    vector<double> ax(n), ay(n), az(n);

    //determine the acceleration
    tie(ax, ay, az) = rhs(t, x, y, z, n, m);
    // ay = rhs(t, y, z, x, n, m);
    // az = rhs(t, z, x, y, n, m);

    //do the iteration step (update the positions)
    for(int i=0; i<n; i++) x[i] += dt * vx[i];
    for(int i=0; i<n; i++) y[i] += dt * vy[i];
    for(int i=0; i<n; i++) z[i] += dt * vz[i];
    for(int i=0; i<n; i++) vx[i] += dt * ax[i];
    for(int i=0; i<n; i++) vy[i] += dt * ay[i];
    for(int i=0; i<n; i++) vz[i] += dt * az[i];
}

void rk4_step_sat(double t, double dt, vector<double> x, vector<double> y, vector<double> z, vector<double> vx, vector<double> vy, vector<double> vz, vector<double> &xs, vector<double> &ys, vector<double> &zs, vector<double> &vxs, vector<double> &vys, vector<double> &vzs, DGL_sat rhs_sat, int n, vector<double> m, int sat){
    //Initialize vectors for the steps - only one step here!
    vector<double> ax1(sat), ax2(sat), ax3(sat), ax4(sat), tmpxs(sat);
    vector<double> ay1(sat), ay2(sat), ay3(sat), ay4(sat), tmpys(sat);
    vector<double> az1(sat), az2(sat), az3(sat), az4(sat), tmpzs(sat);
    vector<double> vx1(sat), vx2(sat), vx3(sat), vx4(sat);
    vector<double> vy1(sat), vy2(sat), vy3(sat), vy4(sat);
    vector<double> vz1(sat), vz2(sat), vz3(sat), vz4(sat);

    //tmp vars refer to the satellite, do not update the planet positions
    //first rk4 step
    for(int i=0; i<sat; i++) vx1[i] = vxs[i];
    for(int i=0; i<sat; i++) vy1[i] = vys[i];
    for(int i=0; i<sat; i++) vz1[i] = vzs[i];

    //cout << "run sat_acceleration here" << endl;
    tie(ax1, ay1, az1) = rhs_sat(t, xs, ys, zs, x, y, z, n, sat, m);
    //cout << "a1[2] = " << sqrt(pow(ax1[2],2) + pow(ay1[2],2) + pow(az1[2],2)) << endl;

    //second rk4 step
    for(int i=0; i<sat; i++) vx2[i] = vxs[i] + (dt/2.) * ax1[i];
    for(int i=0; i<sat; i++) vy2[i] = vys[i] + (dt/2.) * ay1[i];
    for(int i=0; i<sat; i++) vz2[i] = vzs[i] + (dt/2.) * az1[i];
    for(int i=0; i<sat; i++) tmpxs[i] = xs[i] + (dt/2.) * vx1[i];
    for(int i=0; i<sat; i++) tmpys[i] = ys[i] + (dt/2.) * vy1[i];
    for(int i=0; i<sat; i++) tmpzs[i] = zs[i] + (dt/2.) * vz1[i];

    tie(ax2, ay2, az2) = rhs_sat(t+dt/2., tmpxs, tmpys, tmpzs, x, y, z, n, sat, m);

    //third rk4 step
    for(int i=0; i<sat; i++) vx3[i] = vxs[i] + (dt/2.) * ax2[i];
    for(int i=0; i<sat; i++) vy3[i] = vys[i] + (dt/2.) * ay2[i];
    for(int i=0; i<sat; i++) vz3[i] = vzs[i] + (dt/2.) * az2[i];
    for(int i=0; i<sat; i++) tmpxs[i] = xs[i] + (dt/2.) * vx2[i];
    for(int i=0; i<sat; i++) tmpys[i] = ys[i] + (dt/2.) * vy2[i];
    for(int i=0; i<sat; i++) tmpzs[i] = zs[i] + (dt/2.) * vz2[i];

    tie(ax3, ay3, az3) = rhs_sat(t+dt/2., tmpxs, tmpys, tmpzs, x, y, z, n, sat, m);

    //fourth rk4 step
    for(int i=0; i<sat; i++) vx4[i] = vxs[i] + dt * ax3[i];
    for(int i=0; i<sat; i++) vy4[i] = vys[i] + dt * ay3[i];
    for(int i=0; i<sat; i++) vz4[i] = vzs[i] + dt * az3[i];
    for(int i=0; i<sat; i++) tmpxs[i] = xs[i] + dt * vx3[i];
    for(int i=0; i<sat; i++) tmpys[i] = ys[i] + dt * vy3[i];
    for(int i=0; i<sat; i++) tmpzs[i] = zs[i] + dt * vz3[i];

    tie(ax4, ay4, az4) = rhs_sat(t+dt, tmpxs, tmpys, tmpzs, x, y, z, n, sat, m);

    //do the iteration step (update the positions)
    for(int i=0; i<sat; i++) vxs[i] += (dt/6.) * (ax1[i] + 2.*ax2[i] + 2.*ax3[i] + ax4[i]);
    for(int i=0; i<sat; i++) vys[i] += (dt/6.) * (ay1[i] + 2.*ay2[i] + 2.*ay3[i] + ay4[i]);
    for(int i=0; i<sat; i++) vzs[i] += (dt/6.) * (az1[i] + 2.*az2[i] + 2.*az3[i] + az4[i]);
    for(int i=0; i<sat; i++)  xs[i] += (dt/6.) * (vx1[i] + 2.*vx2[i] + 2.*vx3[i] + vx4[i]);
    for(int i=0; i<sat; i++)  ys[i] += (dt/6.) * (vy1[i] + 2.*vy2[i] + 2.*vy3[i] + vy4[i]);
    for(int i=0; i<sat; i++)  zs[i] += (dt/6.) * (vz1[i] + 2.*vz2[i] + 2.*vz3[i] + vz4[i]);
}

void rk4_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, DGL rhs, int n, vector<double> m){
    //Initialize vectors for the steps - only one step here!
    vector<double> ax1(n), ax2(n), ax3(n), ax4(n), tmpx(n);
    vector<double> ay1(n), ay2(n), ay3(n), ay4(n), tmpy(n);
    vector<double> az1(n), az2(n), az3(n), az4(n), tmpz(n);
    vector<double> vx1(n), vx2(n), vx3(n), vx4(n);
    vector<double> vy1(n), vy2(n), vy3(n), vy4(n);
    vector<double> vz1(n), vz2(n), vz3(n), vz4(n);

    //first rk4 step
    for(int i=0; i<n; i++) vx1[i] = vx[i];
    for(int i=0; i<n; i++) vy1[i] = vy[i];
    for(int i=0; i<n; i++) vz1[i] = vz[i];

    tie(ax1, ay1, az1) = rhs(t, x, y, z, n, m);

    //second rk4 step
    for(int i=0; i<n; i++) vx2[i] = vx[i] + (dt/2.) * ax1[i];
    for(int i=0; i<n; i++) vy2[i] = vy[i] + (dt/2.) * ay1[i];
    for(int i=0; i<n; i++) vz2[i] = vz[i] + (dt/2.) * az1[i];
    for(int i=0; i<n; i++) tmpx[i] = x[i] + (dt/2.) * vx1[i];
    for(int i=0; i<n; i++) tmpy[i] = y[i] + (dt/2.) * vy1[i];
    for(int i=0; i<n; i++) tmpz[i] = z[i] + (dt/2.) * vz1[i];

    tie(ax2, ay2, az2) = rhs(t+dt/2., tmpx, tmpy, tmpz, n, m);

    //third rk4 step
    for(int i=0; i<n; i++) vx3[i] = vx[i] + (dt/2.) * ax2[i];
    for(int i=0; i<n; i++) vy3[i] = vy[i] + (dt/2.) * ay2[i];
    for(int i=0; i<n; i++) vz3[i] = vz[i] + (dt/2.) * az2[i];
    for(int i=0; i<n; i++) tmpx[i] = x[i] + (dt/2.) * vx2[i];
    for(int i=0; i<n; i++) tmpy[i] = y[i] + (dt/2.) * vy2[i];
    for(int i=0; i<n; i++) tmpz[i] = z[i] + (dt/2.) * vz2[i];

    tie(ax3, ay3, az3) = rhs(t+dt/2., tmpx, tmpy, tmpz, n, m);

    //fourth rk4 step
    for(int i=0; i<n; i++) vx4[i] = vx[i] + dt * ax3[i];
    for(int i=0; i<n; i++) vy4[i] = vy[i] + dt * ay3[i];
    for(int i=0; i<n; i++) vz4[i] = vz[i] + dt * az3[i];
    for(int i=0; i<n; i++) tmpx[i] = x[i] + dt * vx3[i];
    for(int i=0; i<n; i++) tmpy[i] = y[i] + dt * vy3[i];
    for(int i=0; i<n; i++) tmpz[i] = z[i] + dt * vz3[i];

    tie(ax4, ay4, az4) = rhs(t+dt, tmpx, tmpy, tmpz, n, m);

    //do the iteration step (update the positions)
    for(int i=0; i<n; i++) vx[i] += (dt/6.) * (ax1[i] + 2.*ax2[i] + 2.*ax3[i] + ax4[i]);
    for(int i=0; i<n; i++) vy[i] += (dt/6.) * (ay1[i] + 2.*ay2[i] + 2.*ay3[i] + ay4[i]);
    for(int i=0; i<n; i++) vz[i] += (dt/6.) * (az1[i] + 2.*az2[i] + 2.*az3[i] + az4[i]);
    for(int i=0; i<n; i++)  x[i] += (dt/6.) * (vx1[i] + 2.*vx2[i] + 2.*vx3[i] + vx4[i]);
    for(int i=0; i<n; i++)  y[i] += (dt/6.) * (vy1[i] + 2.*vy2[i] + 2.*vy3[i] + vy4[i]);
    for(int i=0; i<n; i++)  z[i] += (dt/6.) * (vz1[i] + 2.*vz2[i] + 2.*vz3[i] + vz4[i]);
}

void lf_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, DGL rhs, int n, vector<double> m){
    //Initialize vectors for the steps - only one step here!
    vector<double> ax(n), ay(n), az(n);

    //lf step - calculate derivative of v at n-th position
    tie(ax, ay, az) = rhs(t, x, y, z, n, m);

    //calculate n+1/2 value of v
    for(int i=0; i<n; i++) vx[i] += dt * ax[i];
    for(int i=0; i<n; i++) vy[i] += dt * ay[i];
    for(int i=0; i<n; i++) vz[i] += dt * az[i];

    //do the iteration step (update the positions (n+1))
    for(int i=0; i<n; i++) x[i] += dt * vx[i];
    for(int i=0; i<n; i++) y[i] += dt * vy[i];
    for(int i=0; i<n; i++) z[i] += dt * vz[i];
}

void driver(double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, int n, vector<double> m, vector<double> &r, Step_function step, string command){
    //Create and open output file
    fstream file;
    file.open(command+"-solution.csv", ios::out);
    file.precision(16);

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

        //Calculate next timestep
        step(t, dt, x, y, z, vx, vy, vz, acceleration, n, m);

        //update time - so the loop will have a chance to end
        t += dt;
        count ++;
    }

    //close the output file after the iterations are done
    file.close();
}

void sat_driver(double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> m, vector<double> &r, vector<double> &xs, vector<double> &ys, vector<double> &zs, vector<double> &vxs, vector<double> &vys, vector<double> &vzs, vector<double> &ms, int n, Step_function step, double lower, double upper, int sat, vector<double> &t_maxdist, vector<double> &maxdist, vector<double> &v_maxdist, vector<double> &v0_sat, bool out, Step_satellite step_sat){
    //Create and open output file
    fstream file;
    if (out){
        file.open("Satellites.csv", ios::out);
        file.precision(16);
    }

    vector<double> dist = {};
    vector<double> v_dist = {};

    //clear outputvalues in the end since the no of satellites changes
    maxdist.erase(maxdist.begin(), maxdist.end());
    v_maxdist.erase(v_maxdist.begin(), v_maxdist.end());
    t_maxdist.erase(t_maxdist.begin(), t_maxdist.end());
    v0_sat.erase(v0_sat.begin(), v0_sat.end());

    for(int i=0; i<sat; i++) dist.push_back(0.);
    for(int i=0; i<sat; i++) v_dist.push_back(0.);
    for(int i=0; i<sat; i++) maxdist.push_back(0.);
    for(int i=0; i<sat; i++) v_maxdist.push_back(sqrt(pow(vxs[i],2) + pow(vys[i],2) + pow(vzs[i],2)));
    for(int i=0; i<sat; i++) t_maxdist.push_back(0.);
    for(int i=0; i<sat; i++) v0_sat.push_back(sqrt(pow(vxs[i],2) + pow(vys[i],2) + pow(vzs[i],2)));

    int counter = 0; //counts the number of crashed satellites in each timestep
    //cout << "entered sat driver and set startvalues; t = " << t << endl;
    // cout << "vdist0[2] = " << v0_sat[2] << endl;
    //loop that iterates up to a certain chosen time (end)
    while((t_end - t) > DBL_EPSILON){
        //cout << "sat while " << t << endl;
        //since we work in ICRF2 frame the "Schwerpunkt" is an fixed point we can reference the distance to.
        //Calculate next timestep for satellites
        //cout << "run sat_acceleration after" << endl;
        step_sat(t, dt, x, y, z, vx, vy, vz, xs, ys, zs, vxs, vys, vzs, sat_acceleration, n, m, sat);
        //cout << "run sat_acceleration before" << endl;

        //calculate the velocities for satellites at their distances
        for(int i=0; i<sat; i++) v_dist[i] = sqrt(pow(vxs[i],2) + pow(vys[i],2) + pow(vzs[i],2));
        cout << "vdist[2] = " << v_dist[2] << endl;
        for(int i=0; i<sat; i++) dist[i] = sqrt(pow(xs[i],2) + pow(ys[i],2) + pow(zs[i],2));

        //if v at dist lower than v_maxdist: update v (since should be min at maximum distance, that we know for sure)
        for(int i=0; i<sat; i++){
            if ((v_maxdist[i]-v_dist[i]) > DBL_EPSILON){
                v_maxdist[i] = v_dist[i];
                maxdist[i] = dist[i];
                t_maxdist[i] = t;
                // cout << "I never get here----------------------v_dist increases due to very high accelerations" << endl; //but acceleration function seems to be correct
                //mabye sat_initialisation function that is not correct?
            }
        }

        for(int i=(sat-1); i=0; i--){
            //if a satellite collided with an object remove it from further calculations
            if (crash_check(x, y, z, r, xs, ys, zs, n, i)){
                xs.erase(xs.begin()+i);
                ys.erase(ys.begin()+i);
                zs.erase(zs.begin()+i);
                vxs.erase(vxs.begin()+i);
                vys.erase(vys.begin()+i);
                vzs.erase(vzs.begin()+i);
                ms.erase(ms.begin()+i);

                dist.erase(dist.begin()+i);
                maxdist.erase(maxdist.begin()+i);
                v_dist.erase(v_dist.begin()+i);
                v_maxdist.erase(v_maxdist.begin()+i);
                t_maxdist.erase(t_maxdist.begin()+i);
                v0_sat.erase(v0_sat.begin()+i);
                counter++;
            }
        }

        //update number of remaining satellites
        sat -= counter;

        //Calculate next timestep for objects
        step(t, dt, x, y, z, vx, vy, vz, acceleration, n, m);

        //update time
        t += dt;
        counter = 0;
    }
    cout << "exit sat while" << endl;

    //delete satellites if they did not return in a suitable interval
    for(int i=(sat-1); i=0; i--){
        if (((upper-maxdist[i]) < DBL_EPSILON) || ((maxdist[i] - lower) < DBL_EPSILON)){
            t_maxdist.erase(t_maxdist.begin()+i);
            maxdist.erase(maxdist.begin()+i);
            v_maxdist.erase(v_maxdist.begin()+i);
            v0_sat.erase(v0_sat.begin()+i);
        }
    }

    if (out){
        for(int i=0; i<sat; i++){
            //if a satellite returned in the given interval output that satellite to satellite file
            file << t_maxdist[i] << "; " << maxdist[i] << "; " << v_maxdist[i] << "; " << v0_sat[i] << endl; //need not to check lower upper here
        }
    }
    //close the output file after the iterations are done
    file.close();
    cout << "exit sat driver" << endl;
}

vector<double> check_for_boundaries(int precision, int n, double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> &m, vector<double> &r, vector<double> &xs, vector<double> &ys, vector<double> &zs, vector<double> &vxs, vector<double> &vys, vector<double> &vzs, vector<double> &ms, double upper, double lower, string name, int startobject){
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
    double v0 = 10.;
    double v_min = v0; //initial velocity and calculate below
    double v_max = 0.;
    bool out = false;
    bool final = false;

    //inserted for testing
    precision = 0;

    vector<double> t_maxdist = {};
    vector<double> v_maxdist = {};
    vector<double> maxdist = {};
    vector<double> v0_sat = {};

    while(counter < (precision+1)){ //determines the precision of v0 in number of digits after the comma
        cout << "while loop " << counter << endl;
        //Old objects and satellites are destroyed and new ones created
        initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
        cout << "initialized objects" << endl;
        prefactor = initialize_satellites(final, counter, v_min, v_max, x, y, z, vx, vy, vz, r, xs, ys, zs, vxs, vys, vzs, ms, sat, startobject);
        cout << "initialized satellites; t = " << t << "; dt = " << dt << endl;

        //run programm to the end and get t_maxdist, maxdist, v_maxdist and v0_sat of satellites that returned in suitable interval
        sat_driver(t, t_end, dt, x, y, z, vx, vy, vz, m, r, xs, ys, zs, vxs, vys, vzs, ms, n, rk4_step, lower, upper, prefactor*sat, t_maxdist, maxdist, v_maxdist, v0_sat, out, rk4_step_sat);

        //know from satellite driver how many satellites made it and know their initial velocities
        v_min = findmin(v0_sat);
        v_max = findmax(v0_sat);

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

void calc_sat(vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> m, vector<double> &r, int n, Step_function step, Step_satellite step_sat, string name){
    //Create and open output file
    fstream file;
    string input = "Orbits.csv";
    file.open("Satellites.csv", ios::out);
    file.precision(16);
    bool out = true;
    bool final = true;
    double v_min, v_max;
    //double counter = 4;

    double t = 0.;
    double t_end = 20.;
    //double dt = pow(2,-13);
    double dt = pow(2,-13); //--just for testing purposes
    int startobject = 2; //no. planet -1; (Earth = 2)
    int endobject = 8;
    int precision = 4;
    int satdummy = 10;
    int sat; //aka prefactor at another point

    vector<string> tmp = {};
    double upper = 101 * upperlower(tmp, input, endobject)[0]; //read 1% difference from max min orbit
    double lower = 101 * upperlower(tmp, input, endobject)[1];
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

    vector<double> boundaries = check_for_boundaries(precision, n, t, t_end, dt, x, y, z, vx, vy, vz, m, r, xs, ys, zs, vxs, vys, vzs, ms, upper, lower, name, startobject);
    v_min = boundaries[0];
    v_max = boundaries[1];

    cout << "boundaries[0] = " << boundaries[0] << "; boundaries[1] = " << boundaries[1] << endl;

    //initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
    //sat = initialize_satellites(final, precision, v_min, v_max, x, y, z, vx, vy, vz, r, xs, ys, zs, vxs, vys, vzs, ms, satdummy, startobject);
    //sat_driver(t, t_end, dt, x, y, z, vx, vy, vz, m, r, xs, ys, zs, vxs, vys, vzs, ms, n, rk4_step, lower, upper, sat, t_maxdist, maxdist, v_maxdist, v0_sat, out, rk4_step_sat);
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
    double t_end = 20.;           //final time
    double dt = pow(2.,-13);     //time steps
    double t = 0.;

    string name = "Input2.csv";

    if(fileexists(name)){
        if (command == "fwd"){  // forward euler
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
            driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, fwd_step, command);
        }
        else if (command == "rk4"){ // Runge Kutta 4
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
            driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, rk4_step, command);
        }
        else if (command == "lf"){ // leap frog
            initialize_objects(n, x, y, z, vx, vy, vz, m, r, name);
            driver(t, t_end, dt, x, y, z, vx, vy, vz, n, m, r, lf_step, command);
        }
        else if (command == "sat"){ // satellites
            calc_sat(x, y, z, vx, vy, vz, m, r, n, rk4_step, rk4_step_sat, name);
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
