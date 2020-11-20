#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include "LU.h"

using namespace std;

double Temp(double x, double t){
  return 200+8000*t*(0.01-x);
}


double e(double x, double t, double dt){
  return 3000*exp((-6000)/(Temp(x,t)*1.));
}

double u_seconde(vector<double> x,vector<double> rho, int j){
  if(j==0){
    return 2*((rho[j+1]-rho[j])/pow(x[j+1]-x[j],2.));
  }
  if(j==x.size()-1){
    return 2*((rho[j-1]-rho[j])/pow(x[j]-x[j-1],2.));
  }

  double a=(rho[j+1]-rho[j])/(x[j+1]-x[j]);
  double b=(rho[j-1]-rho[j])/(x[j]-x[j-1]);
  return 2*(a+b)/(x[j+1]-x[j-1]);
}  

vector<vector<double>> remplissage(int n, vector<double> x, double c, double v, double min){
  vector<double> kr(n+1);
  vector<double> metr(n+1);
  vector<vector<double>> kk(n+1,vector<double>(n+1,0.));
  
  for(int j=0;j<n+1;j++){
    metr[j]=max(abs(u_seconde(x[j],c,v)),min);
    // cout<<"vl"<<abs(u_seconde(x_j,c,v))<<endl;
  }

  //cout<<"----------------------------------------------"<<endl;

  kr[0]=0.;
  for(int k=1;k<n+1;k++){
    kr[k]=(metr[k-1]+metr[k])/2.;
  }

  kk[0][0]=pow(10,12);
  kk[n][n]=pow(10,12);
  for(int i=1;i<n;i++){
    kk[i][i]=(kr[i]+kr[i+1]);  
  }

  for(int i=1;i<n;i++){
    kk[i][i+1]=-kr[i+1];
    kk[i][i-1]=-kr[i];
  }
  return kk;
}



int main(){
  int B=3000;
  int Ta=6000;
  int p=100;
  int n=100;
  int zmax=30;
  double rhov=1500;
  double rhop=1000;
  double dt(0), dx(0);
  vector<double> x(n+1,0.),b(n+1,0.),xn(n+1,0.);
  vector<double> rho(n+1),rhon(n+1);
  vector<vector<double>> kk(n+1,vector<double>(n+1,0.));

  
  dt=10/(p*1.);
  dx=0.01/(n*1.);
  ofstream mon_flux;
  string name_file="Result"+to_string(0);
  mon_flux.open(name_file, ios::out); // Ouvre un fichier appelé name_file
  cout.precision(15);
  
  for(int i=0;i<n+1;i++){
    rho[i]=rhov;
    x[i]=i*dx;
    mon_flux << i*dx << " " << rho[i] << endl;
  }
  b[n]=pow(10,12)*1.;
  mon_flux.close();
 
  
  for (int j=1;j<p;j++){
    ofstream mon_flux;
    name_file="Result"+to_string(j);
    name_file="Result"+to_string(k);
    mon_flux.open(name_file, ios::out);
    for(int z=0; z<zmax+1;z++ ){
      for(int k=0;k<n+1;k++){
	rhon[k]=(1-e(x[k],j*dt,dt))*rho[k]+1000*e(x[k],j*dt,dt);
      }
      kk=remplissage(n,x,min);
      xn=resol(kk,b,n);
      cout<<"error_1= "<<max_error_vector(x,xn,n)<<endl;
      if(max_error_vector(x,xn,n)<epsilon){
	break;
      }
      x=xn;
    }
    for(int k=0;k<n+1;k++){
      mon_flux << xn[k] << " " << rhon[k] << endl;
    }
    rho=rhon;
    mon_flux.close();

  }
  //cout<<"la valeur et"<<rho.size()<<endl;
  


  return 0;
}
