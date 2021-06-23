#pragma once

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>
#include <cfloat>
#include <tuple>
#include <omp.h>

using namespace std;

typedef void (* Step_function)(double, double, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, tuple<vector<double>,vector<double>,vector<double>>(double, vector<double>, vector<double>, vector<double>, int, vector<double>), int, vector<double>);
typedef tuple<vector<double>,vector<double>,vector<double>> (DGL)(double t, vector<double> x, vector<double> y, vector<double> z, int n, vector<double> m);

tuple<vector<double>,vector<double>,vector<double>> acceleration(double t, vector<double> x, vector<double> y, vector<double> z, int n, vector<double> m){
  double Matrix[n][n];

  vector<double> ax(n, 0.);
  vector<double> ay(n, 0.);
  vector<double> az(n, 0.);

  #pragma omp parallel for
  for(int i=0; i<n; i++){
    for(int j=0; j<i; j++) {
      if ( (i > 9) && (j > 9) ) Matrix[i][j] = 0.;
      else{
        Matrix[i][j] = pow((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]),-1.5);
      }
    }
  }

  for(int i=0; i<n; i++){
    Matrix[i][i] = 0.;
    for(int j=i+1; j<n; j++) Matrix[i][j] = Matrix[j][i];
    for(int j=0; j<n; j++) {
      //if (isinf(Matrix[i][j]) == true) Matrix[i][j] = 0;
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

void set_satellite(vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, double vsat, double r){
  const double v0 = 6.179554659505929;
  const double h0 = 1.016334532486522;
  const double g0 = -38.21962069010524;

  double hs = pow(x[3]*x[3]+y[3]*y[3]+z[3]*z[3], 0.5);
  double ve = pow(vx[3]*vx[3]+vy[3]*vy[3]+vz[3]*vz[3], 0.5);
  double gs = -4*M_PI*M_PI/(hs*hs);

  vsat += ve - v0;

  // double Deltav = abs(v0-pow(v0*v0-2*(g0*h0-gs*hs),0.5));
  // if (h0 > hs){vsat += Deltav;}
  // else{vsat -= Deltav;}

  x[10] = x[3] + vx[3]*r/ve;
  y[10] = y[3] + vy[3]*r/ve;
  z[10] = z[3] + vz[3]*r/ve;

  vx[10] = vx[3]*vsat/ve;
  vy[10] = vy[3]*vsat/ve;
  vz[10] = vz[3]*vsat/ve;
}

double angle(double x1, double y1, double z1, double x2, double y2, double z2){
  return acos((x1*x2+y1*y2+z1*z2)/(pow(x1*x1+y1*y1+z1*z1,0.5)*pow(x2*x2+y2*y2+z2*z2,0.5)));
}

double distance(double x1, double y1, double z1, double x2, double y2, double z2){
  return pow((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2),0.5);
}

bool crash_check(double t, vector<double> x, vector<double> y, vector<double> z, vector<double> r, vector<double> xs, vector<double> ys, vector<double> zs, int n, int satno){
    vector<double> dist = {};
    bool counter = false;
    vector<string> objects = {"Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"};

    //if distance of sattelite number satno to the i-th object is smaller than its radius return true
    for(int i=0; i<n; i++) dist.push_back(sqrt(pow(x[i]-xs[satno],2) + pow(y[i]-ys[satno],2) + pow(z[i]-zs[satno],2)));
    for(int i=0; i<n; i++){
        if ((dist[i] - r[i]) < DBL_EPSILON){
            cout <<"t = " << t << "; Satellite No. " << satno << " crashed with " << objects[i] << endl;
            counter = true;
        }
    }
    return counter;
}

double findmin(vector<double> v){
    double min = DBL_MAX;
    for(double i : v) if ((min-i) > DBL_EPSILON) min = i;
    return min;
}

double findmax(vector<double> v){
    double max = -M_PI;
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

vector<vector<double>> set_values(vector<vector<double>> values, int n, string name){
  vector<string> help = {};
  read_file(help, name);
  int counter;
  string s, tmp_s;
  string semi = ";";  //define delimiter which shall be searched for

  int k = values.size();
  //iterate over auxiliary vector
  for(int i=0; i<n; i++){
    s = help[i];
    //iterate over each double in one line
    for(int l=0; l<k; l++){
      //iterate over each line --- mabye replace by find_first_of (see cppreference)
      for(int j=0; j<s.length(); j++){
        if(s[j] == semi[0]){
          counter = j;
          break;
        }
        // if the end of line has no ";"
        if(j == s.length()-1){
            counter = j-1;
            break;
        }
      }
      //cut off the first part of s
      tmp_s = s.substr(0, counter);

      //set s to the remaining string
      if ((counter-2) < s.length()) s = s.substr(counter+1, s.length()-1);
      //set values
      values[l][i] = stod(tmp_s);
    }
  }
  return values;
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
    r.erase(r.begin(), r.end());

    //Change size of vectors
    x.resize(n);
    y.resize(n);
    z.resize(n);
    vx.resize(n);
    vy.resize(n);
    vz.resize(n);
    m.resize(n);
    r.resize(n);

    vector<vector<double>> values = {x, y, z, vx, vy, vz, m, r};
    values = set_values(values, n, name);
    x = values[0];
    y = values[1];
    z = values[2];
    vx = values[3];
    vy = values[4];
    vz = values[5];
    m = values[6];
    r = values[7];
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
    for(int i=0; i<n; i++) {
      vx2[i] = vx[i] + (dt/2.) * ax1[i];
      vy2[i] = vy[i] + (dt/2.) * ay1[i];
      vz2[i] = vz[i] + (dt/2.) * az1[i];
      tmpx[i] = x[i] + (dt/2.) * vx1[i];
      tmpy[i] = y[i] + (dt/2.) * vy1[i];
      tmpz[i] = z[i] + (dt/2.) * vz1[i];
    }
    tie(ax2, ay2, az2) = rhs(t+dt/2., tmpx, tmpy, tmpz, n, m);

    //third rk4 step
    for(int i=0; i<n; i++) {
      vx3[i] = vx[i] + (dt/2.) * ax2[i];
      vy3[i] = vy[i] + (dt/2.) * ay2[i];
      vz3[i] = vz[i] + (dt/2.) * az2[i];
      tmpx[i] = x[i] + (dt/2.) * vx2[i];
      tmpy[i] = y[i] + (dt/2.) * vy2[i];
      tmpz[i] = z[i] + (dt/2.) * vz2[i];
    }
    tie(ax3, ay3, az3) = rhs(t+dt/2., tmpx, tmpy, tmpz, n, m);

    //fourth rk4 step
    for(int i=0; i<n; i++) {
      vx4[i] = vx[i] + dt * ax3[i];
      vy4[i] = vy[i] + dt * ay3[i];
      vz4[i] = vz[i] + dt * az3[i];
      tmpx[i] = x[i] + dt * vx3[i];
      tmpy[i] = y[i] + dt * vy3[i];
      tmpz[i] = z[i] + dt * vz3[i];
    }
    tie(ax4, ay4, az4) = rhs(t+dt, tmpx, tmpy, tmpz, n, m);

    //do the iteration step (update the positions)
    for(int i=0; i<n; i++) {
      vx[i] += (dt/6.) * (ax1[i] + 2.*ax2[i] + 2.*ax3[i] + ax4[i]);
      vy[i] += (dt/6.) * (ay1[i] + 2.*ay2[i] + 2.*ay3[i] + ay4[i]);
      vz[i] += (dt/6.) * (az1[i] + 2.*az2[i] + 2.*az3[i] + az4[i]);
      x[i] += (dt/6.) * (vx1[i] + 2.*vx2[i] + 2.*vx3[i] + vx4[i]);
      y[i] += (dt/6.) * (vy1[i] + 2.*vy2[i] + 2.*vy3[i] + vy4[i]);
      z[i] += (dt/6.) * (vz1[i] + 2.*vz2[i] + 2.*vz3[i] + vz4[i]);
    }
}

void rk5_step(double t, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, DGL rhs, int n, vector<double> m){
    //Initialize vectors for the steps - only one step here!
    vector<double> ax1(n), ax2(n), ax3(n), ax4(n), ax5(n), ax6(n), tmpx(n);
    vector<double> ay1(n), ay2(n), ay3(n), ay4(n), ay5(n), ay6(n), tmpy(n);
    vector<double> az1(n), az2(n), az3(n), az4(n), az5(n), az6(n), tmpz(n);
    vector<double> vx1(n), vx2(n), vx3(n), vx4(n), vx5(n), vx6(n);
    vector<double> vy1(n), vy2(n), vy3(n), vy4(n), vy5(n), vy6(n);
    vector<double> vz1(n), vz2(n), vz3(n), vz4(n), vz5(n), vz6(n);
    vector<double> x_4(n), y_4(n), z_4(n), vx_4(n), vy_4(n), vz_4(n);
    vector<double> x_5(n), y_5(n), z_5(n), vx_5(n), vy_5(n), vz_5(n);

    //first rk4 step
    for(int i=0; i<n; i++) vx1[i] = vx[i];
    for(int i=0; i<n; i++) vy1[i] = vy[i];
    for(int i=0; i<n; i++) vz1[i] = vz[i];

    tie(ax1, ay1, az1) = rhs(t, x, y, z, n, m);

    //second rk4 step
    for(int i=0; i<n; i++) {
      vx2[i] = vx[i] + (dt/5.) * ax1[i];
      vy2[i] = vy[i] + (dt/5.) * ay1[i];
      vz2[i] = vz[i] + (dt/5.) * az1[i];
      tmpx[i] = x[i] + (dt/5.) * vx[i];
      tmpy[i] = y[i] + (dt/5.) * vy[i];
      tmpz[i] = z[i] + (dt/5.) * vz[i];
    }

    tie(ax2, ay2, az2) = rhs(t+dt/5., tmpx, tmpy, tmpz, n, m);

    //third rk4 step
    for(int i=0; i<n; i++) {
      vx3[i] = vx[i] + ((3./40.) * ax1[i] + (9./40.) * ax2[i])*dt;
      vy3[i] = vy[i] + ((3./40.) * ay1[i] + (9./40.) * ay2[i])*dt;
      vz3[i] = vz[i] + ((3./40.) * az1[i] + (9./40.) * az2[i])*dt;
      tmpx[i] = x[i] + ((3./40.) * vx1[i] + (9./40.) * vx2[i])*dt;
      tmpy[i] = y[i] + ((3./40.) * vy1[i] + (9./40.) * vy2[i])*dt;
      tmpz[i] = z[i] + ((3./40.) * vz1[i] + (9./40.) * vz2[i])*dt;
    }

    tie(ax3, ay3, az3) = rhs(t+dt*(3./10.), tmpx, tmpy, tmpz, n, m);

    //fourth rk4 step
    for(int i=0; i<n; i++) {
      vx4[i] = vx[i] + ((3./10.) * ax1[i] - (9./10.) * ax2[i] + (6./5.) * ax3[i])*dt;
      vy4[i] = vy[i] + ((3./10.) * ay1[i] - (9./10.) * ay2[i] + (6./5.) * ay3[i])*dt;
      vz4[i] = vz[i] + ((3./10.) * az1[i] - (9./10.) * az2[i] + (6./5.) * az3[i])*dt;
      tmpx[i] = x[i] + ((3./10.) * vx1[i] - (9./10.) * vx2[i] + (6./5.) * vx3[i])*dt;
      tmpy[i] = y[i] + ((3./10.) * vy1[i] - (9./10.) * vy2[i] + (6./5.) * vy3[i])*dt;
      tmpz[i] = z[i] + ((3./10.) * vz1[i] - (9./10.) * vz2[i] + (6./5.) * vz3[i])*dt;
    }

    tie(ax4, ay4, az4) = rhs(t+dt*(3./5.), tmpx, tmpy, tmpz, n, m);

    //fifth rk4 step
    for(int i=0; i<n; i++) {
      vx5[i] = vx[i] + (-(11./54.) * ax1[i] + (5./2.) * ax2[i] - (70./27.) * ax3[i] + (35./27.) * ax4[i])*dt;
      vy5[i] = vy[i] + (-(11./54.) * ay1[i] + (5./2.) * ay2[i] - (70./27.) * ay3[i] + (35./27.) * ay4[i])*dt;
      vz5[i] = vz[i] + (-(11./54.) * az1[i] + (5./2.) * az2[i] - (70./27.) * az3[i] + (35./27.) * az4[i])*dt;
      tmpx[i] = x[i] + (-(11./54.) * vx1[i] + (5./2.) * vx2[i] - (70./27.) * vx3[i] + (35./27.) * vx4[i])*dt;
      tmpy[i] = y[i] + (-(11./54.) * vy1[i] + (5./2.) * vy2[i] - (70./27.) * vy3[i] + (35./27.) * vy4[i])*dt;
      tmpz[i] = z[i] + (-(11./54.) * vz1[i] + (5./2.) * vz2[i] - (70./27.) * vz3[i] + (35./27.) * vz4[i])*dt;
    }

    tie(ax5, ay5, az5) = rhs(t+dt, tmpx, tmpy, tmpz, n, m);

    //sixth rk4 step
    for(int i=0; i<n; i++) {
      vx6[i] = vx[i] + ((1631./55296.) * ax1[i] + (175./512.) * ax2[i] + (575./13824.) * ax3[i] + (44275./110592.) * ax4[i] + (253./4096.) * ax5[i])*dt;
      vy6[i] = vy[i] + ((1631./55296.) * ay1[i] + (175./512.) * ay2[i] + (575./13824.) * ay3[i] + (44275./110592.) * ay4[i] + (253./4096.) * ay5[i])*dt;
      vz6[i] = vz[i] + ((1631./55296.) * az1[i] + (175./512.) * az2[i] + (575./13824.) * az3[i] + (44275./110592.) * az4[i] + (253./4096.) * az5[i])*dt;
      tmpx[i] = x[i] + ((1631./55296.) * vx1[i] + (175./512.) * vx2[i] + (575./13824.) * vx3[i] + (44275./110592.) * vx4[i] + (253./4096.) * vx5[i])*dt;
      tmpy[i] = y[i] + ((1631./55296.) * vy1[i] + (175./512.) * vy2[i] + (575./13824.) * vy3[i] + (44275./110592.) * vy4[i] + (253./4096.) * vy5[i])*dt;
      tmpz[i] = z[i] + ((1631./55296.) * vz1[i] + (175./512.) * vz2[i] + (575./13824.) * vz3[i] + (44275./110592.) * vz4[i] + (253./4096.) * vz5[i])*dt;
    }

    tie(ax6, ay6, az6) = rhs(t+dt*(7./8.), tmpx, tmpy, tmpz, n, m);

    //do the iteration step (update the positions)
    for(int i=0; i<n; i++) {
      vx_5[i] = ((37./378.) * ax1[i] + (250./621.) * ax3[i] + (125./594.) * ax4[i] + (512./1771.) * ax6[i])*dt;
      vy_5[i] = ((37./378.) * ay1[i] + (250./621.) * ay3[i] + (125./594.) * ay4[i] + (512./1771.) * ay6[i])*dt;
      vz_5[i] = ((37./378.) * az1[i] + (250./621.) * az3[i] + (125./594.) * az4[i] + (512./1771.) * az6[i])*dt;
      x_5[i]  = ((37./378.) * vx1[i] + (250./621.) * vx3[i] + (125./594.) * vx4[i] + (512./1771.) * vx6[i])*dt;
      y_5[i]  = ((37./378.) * vy1[i] + (250./621.) * vy3[i] + (125./594.) * vy4[i] + (512./1771.) * vy6[i])*dt;
      z_5[i]  = ((37./378.) * vz1[i] + (250./621.) * vz3[i] + (125./594.) * vz4[i] + (512./1771.) * vz6[i])*dt;
    }

    // rk4 scheme
    for(int i=0; i<n; i++) {
      vx_4[i] = ((2825./27648.) * ax1[i] + (18575./48384.) * ax3[i] + (13525./55296.) * ax4[i] + (277./14336.) * ax5[i] + (1./4.) * ax6[i])*dt;
      vy_4[i] = ((2825./27648.) * ay1[i] + (18575./48384.) * ay3[i] + (13525./55296.) * ay4[i] + (277./14336.) * ay5[i] + (1./4.) * ay6[i])*dt;
      vz_4[i] = ((2825./27648.) * az1[i] + (18575./48384.) * az3[i] + (13525./55296.) * az4[i] + (277./14336.) * az5[i] + (1./4.) * az6[i])*dt;
      x_4[i]  = ((2825./27648.) * vx1[i] + (18575./48384.) * vx3[i] + (13525./55296.) * vx4[i] + (277./14336.) * vx5[i] + (1./4.) * vx6[i])*dt;
      y_4[i]  = ((2825./27648.) * vy1[i] + (18575./48384.) * vy3[i] + (13525./55296.) * vy4[i] + (277./14336.) * vy5[i] + (1./4.) * vy6[i])*dt;
      z_4[i]  = ((2825./27648.) * vz1[i] + (18575./48384.) * vz3[i] + (13525./55296.) * vz4[i] + (277./14336.) * vz5[i] + (1./4.) * vz6[i])*dt;
    }

    for(int i=0; i<n; i++) {
      vx[i] = abs(vx_5[i] - vx_4[i]);
      vy[i] = abs(vy_5[i] - vy_4[i]);
      vz[i] = abs(vz_5[i] - vz_4[i]);
      x[i]  = abs(x_5[i]  - x_4[i]);
      y[i]  = abs(y_5[i]  - y_4[i]);
      z[i]  = abs(z_5[i]  - z_4[i]);
    }
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

void sat_driver(double t, double t_end, double dt, vector<double> &x, vector<double> &y, vector<double> &z, vector<double> &vx, vector<double> &vy, vector<double> &vz, vector<double> m, vector<double> &r, vector<double> &xs, vector<double> &ys, vector<double> &zs, vector<double> &vxs, vector<double> &vys, vector<double> &vzs, vector<double> &ms, vector<double> &rs, int n, Step_function step, double lower, double upper, int sat, vector<double> &t_maxdist, vector<double> &maxdist, vector<double> &v_maxdist, vector<double> &v0_sat, bool out){
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

    vector<double> x_old, y_old, z_old, vx_old, vy_old, vz_old;
    double Delta;
    double Delta_aim = 1e-16;

    int counter = 0; //counts the number of crashed satellites in each timestep

    x.insert(x.end(), xs.begin(), xs.end());
    y.insert(y.end(), ys.begin(), ys.end());
    z.insert(z.end(), zs.begin(), zs.end());
    vx.insert(vx.end(), vxs.begin(), vxs.end());
    vy.insert(vy.end(), vys.begin(), vys.end());
    vz.insert(vz.end(), vzs.begin(), vzs.end());
    r.insert(r.end(), rs.begin(), rs.end());
    m.insert(m.end(), ms.begin(), ms.end());
    cout << "sat = " << sat << endl;

    //loop that iterates up to a certain chosen time (end)
    while((t_end - t) > DBL_EPSILON){
        x_old = x;
        y_old = y;
        z_old = z;
        vx_old = vx;
        vy_old = vy;
        vz_old = vz;

        //Calculate next timestep
        rk5_step(t, dt, x_old, y_old, z_old, vx_old, vy_old, vz_old, acceleration, n+sat, m);

        Delta = 0;
        for(int i=0; i<n+sat; i++) {
            Delta += x_old[i]+y_old[i]+z_old[i];
        }

        dt *= pow(Delta_aim/Delta, 1./5.);

        //Calculate next timestep for objects and sats
        step(t, dt, x, y, z, vx, vy, vz, acceleration, n+sat, m);

        for(int i=n; i<n+sat; i++){
            xs[i-n] = x[i];
            ys[i-n] = y[i];
            zs[i-n] = z[i];
            vxs[i-n] = vx[i];
            vys[i-n] = vy[i];
            vzs[i-n] = vz[i];
        }

        //calculate the velocities for satellites at their distances
        for(int i=0; i<sat; i++) v_dist[i] = sqrt(pow(vxs[i],2) + pow(vys[i],2) + pow(vzs[i],2));
        for(int i=0; i<sat; i++) dist[i] = sqrt(pow(xs[i],2) + pow(ys[i],2) + pow(zs[i],2));

        //if v at dist lower than v_maxdist: update v (since should be min at maximum distance, that we know for sure)
        for(int i=0; i<sat; i++){
            if ((v_maxdist[i]-v_dist[i]) > DBL_EPSILON){
                v_maxdist[i] = v_dist[i];
                maxdist[i] = dist[i];
                t_maxdist[i] = t;
            }
        }

        for(int i=sat; i>0; i--){
            //if a satellite collided with an object remove it from further calculations
            if (crash_check(t, x, y, z, r, xs, ys, zs, n, i-1)){
                x.erase(x.begin()+(n+i-1));
                y.erase(y.begin()+(n+i-1));
                z.erase(z.begin()+(n+i-1));
                vx.erase(vx.begin()+(n+i-1));
                vy.erase(vy.begin()+(n+i-1));
                vz.erase(vz.begin()+(n+i-1));
                m.erase(m.begin()+(n+i-1));
                r.erase(r.begin()+(n+i-1));

                dist.erase(dist.begin()+(i-1));
                maxdist.erase(maxdist.begin()+(i-1));
                v_dist.erase(v_dist.begin()+(i-1));
                v_maxdist.erase(v_maxdist.begin()+(i-1));
                t_maxdist.erase(t_maxdist.begin()+(i-1));
                v0_sat.erase(v0_sat.begin()+(i-1));
                counter++;
            }
        }

        //update number of remaining satellites
        sat -= counter;

        //update time
        t += dt;
        counter = 0;
    }
    cout << "exit sat while; sat = " << sat << endl;

    //delete satellites if they did not return in a suitable interval
    for(int i=sat; i>0; i--){
        if (((upper-maxdist[i-1]) < DBL_EPSILON) || ((maxdist[i-1] - lower) < DBL_EPSILON)){
            if((upper-maxdist[i-1]) < DBL_EPSILON) cout << v0_sat[i-1] << " is too high" << endl;
            if((maxdist[i-1] - lower) < DBL_EPSILON) cout << v0_sat[i-1] << " is too low" << endl;
            t_maxdist.erase(t_maxdist.begin()+(i-1));
            maxdist.erase(maxdist.begin()+(i-1));
            v_maxdist.erase(v_maxdist.begin()+(i-1));
            v0_sat.erase(v0_sat.begin()+(i-1));
            counter++;
        }
    }
    sat -= counter;

    if (out){
        for(int i=0; i<sat; i++){
            //if a satellite returned in the given interval output that satellite to satellite file
            file << t_maxdist[i] << "; " << maxdist[i] << "; " << v_maxdist[i] << "; " << v0_sat[i] << endl; //need not to check lower upper here
        }
    }
    //close the output file after the iterations are done
    file.close();
    cout << "exit sat driver; sat = " << sat << endl;
}
