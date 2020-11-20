#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>

using namespace std;

double Temp(double x, double t){
  return 200+8000*t*(0.01-x);
}


double e(double x, double t, double dt){
  return 3000*exp((-6000)/(Temp(x,t)*1.));
}


int main(){
  int B=3000;
  int Ta=6000;
  int n=100;
  int p=100;
  double rhov=1500;
  double arho=1500;
  double brhon=0;
  double dt(0), dx(0);
  vector<double> rho(p),rhon(p);
  dt=10/(n*1.);
  dx=0.01/(p*1.);
  ofstream mon_flux;
  string name_file="Result"+to_string(0);
  mon_flux.open(name_file, ios::out); // Ouvre un fichier appelé name_file
  cout.precision(15);
  
  for(int i=0;i<p;i++){
    rho[i]=rhov;
    mon_flux << i*dx << " " << rho[i] << endl;
  }
  mon_flux.close();
 
  
  for (int k=0;k<p;k++){
    ofstream mon_flux;
    //name_file="Result"+to_string(j);
    name_file="Result"+to_string(k);
    mon_flux.open(name_file, ios::out);
    arho=1500;
    for(int j=1;j<n;j++){
      //rhon[k]=(1-e(k*dx,j*dt,dt))*rho[k]+1000*e(k*dx,j*dt,dt);
      brhon=(1-e(k*dx,j*dt,dt))*arho+1000*e(k*dx,j*dt,dt);
      arho=brhon;
      //mon_flux << k*dx << " " << rhon[k] << endl;
      mon_flux << j*dt << " " << brhon << endl;
    }
    //rho=rhon;
    mon_flux.close();

  }
  //cout<<"la valeur et"<<rho.size()<<endl;
  


  return 0;
}
