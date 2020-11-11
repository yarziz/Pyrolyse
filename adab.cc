#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include "LU.h"

using namespace std;

double u_seconde(double x, double c, double v){
  double a=0;
  a=1/(c*(1-exp(c/v)));
  return a*pow(c/v,2.)*exp((c*x)/v);
}  

double u(double x, double c, double v){
  double a=0;
  a=1/(c*(1-exp(c/v)));
  return a*(exp((c*x)/v)-1)+x/c;}

/*double u_seconde(double x, double c, double v){
  return 1;
  }*/


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


void affichage_matrice(vector<vector<double>> kk,int n){
  string a="";
  for(int j=0;j<n+1;j++){
    a="";
    for(int i=0;i<n+1;i++){
      a=a+" "+to_string(kk[j][i]);
    }
    cout<<a<<endl;
  }
}


void affichage_vector(vector<double> x,int n){
  for(int j=0;j<n+1;j++){
       
    cout<<x[j]<<endl;
  }

}


int main(){
  double c=0.5;
  double v=0.01;
  double min=30;
  double dx=0;
  double x_i=0;
  double epsilon=0.01;
  int n=20;
  int zmax=30;
  vector<double> x(n+1,0.),b(n+1,0.),xn(n+1,0.);
  vector<double> metr(n+1);
  vector<vector<double>> kk(n+1,vector<double>(n+1,0.)),kkchol(n+1,vector<double>(n+1,0.));
  
 
  dx=1/(n*1.);
  for(int j=0;j<n+1;j++){
    x[j]=j*dx;
  }
  
  b[n]=pow(10,12)*1.;

  for(int z=0; z<zmax+1;z++ ){
    kk=remplissage(n,x,c,v,min);
    xn=resol(kk,b,n);
    cout<<"error_1= "<<max_error_vector(x,xn,n)<<endl;
    if(max_error_vector(x,xn,n)<epsilon){
      break;
    }
    x=xn;
  }

  ofstream mon_flux;
  ofstream mon_flux_1;
  string name_file="Result";
  string name_file_1="Result_adap";
  mon_flux.open(name_file, ios::out);
  mon_flux_1.open(name_file_1, ios::out);

  for(int i=0;i<n+1;i++){
    mon_flux<<i*dx<<" "<<u(i*dx,c,v)<<endl;
    mon_flux_1<<xn[i]<<" "<<u(xn[i],c,v)<<endl;
  }

  
  
  
  return 0;
}
