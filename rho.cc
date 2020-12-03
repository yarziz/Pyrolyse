#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include "LU.h"
#include "mod_fonction.h"


using namespace std;

double Temp(double x, double t){
  return 200+8000*t*(0.01-x);
}


double e(double x, double t, double dt){
  return 3000*exp((-6000)/(Temp(x,t)*1.));
}


double ue(double x,double t, double c, double v){
  double a=0;
  a=1/(c*(1-exp(c/v)));
  return a*(exp((c*x)/v)-1)+x/c+t*t;}

double u_seconde(vector<double> x,vector<double> rho, int j,int n){
  if(j==0){
    return 2*((rho[j+1]-rho[j])/pow(x[j+1]-x[j],2.));
  }
  if(j==n){
    return 2*((rho[j-1]-rho[j])/pow(x[j]-x[j-1],2.));
  }

  double a=(rho[j+1]-rho[j])/(x[j+1]-x[j]);
  double b=(rho[j-1]-rho[j])/(x[j]-x[j-1]);
  return 2*(a+b)/(x[j+1]-x[j-1]);
}  

vector<vector<double>> remplissage(int n, vector<double> x,vector<double> rhon, double min){
  vector<double> kr(n+1);
  vector<double> u(n+1),uu(n+1,0);
  vector<double> metr(n+1);
  vector<vector<double>> kk(n+1,vector<double>(n+1,0.));
   double max_u=0;
  
  for(int j=0;j<n+1;j++){
    u[j]=abs(u_seconde(x,rhon,j,n));
    // cout<<"vl"<<abs(u_seconde(x_j,c,v))<<endl;
  }
  max_u=max_error_vector(u,uu,n);


  for(int j=0;j<n+1;j++){
    metr[j]=max(abs(u_seconde(x,rhon,j,n))/max_u,min);
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



vector<double> new_rho(vector<double> x,double rhov,int j,int n,double dt,double c,double v){
  vector<double> rho(n+1,rhov);
  vector<double> rhon(n+1);
  vector<double> un(n+1);
  if(j==0){
    return rho ;
  }
  for(int p=0;p<j;p++){
    for(int k=0;k<n+1;k++){
      //un[k]=u[k]+dt*a*(exp((c*x[k])/v)-1)+x[k]/c+2*j*dt;
      rhon[k]=(1-e(x[k],(p)*dt,dt))*rho[k]+1000*e(x[k],(p)*dt,dt);
    }
    rho=rhon;
  }
  return rhon;
}



int main(){
  int B=3000;
  double c=0.5;
  double v=0.01;
  int Ta=6000;
  int p=100;
  int n=30;
  int zmax=200;
  double rhov=1500;
  double rhop=1000;
  double min=1;
  double epsilon=0.00001;
  double dt(0), dx(0);
  vector<double> x(n+1,0.),xx(n+1,0.),b(n+1,0.),xn(n+1,0.);
  vector<double> rho(n+1,rhov),rhon(n+1),u(n+1,0),un(n+1,0),ui(n+1);
  vector<vector<double>> kk(n+1,vector<double>(n+1,0.));

  
  dt=10/(p*1.);
  dx=0.01/(n*1.);
  ofstream mon_flux;
  string name_file="Result"+to_string(0);
  mon_flux.open(name_file, ios::out); // Ouvre un fichier appelé name_file
  cout.precision(15);
  
  for(int i=0;i<n+1;i++){
    x[i]=i*dx;
    // u(i)=ue(x[i],0,c,v);
    //mon_flux << i*dx << " " << u[i] << endl;
    mon_flux << i*dx << " " << rho[i] << endl;
  }
  b[n]=pow(10,10)*1.;
  mon_flux.close();
  //ui=u
  
  for (int j=0;j<p+1;j++){
    ofstream mon_flux;
    name_file="Result"+to_string(j+1);
    mon_flux.open(name_file, ios::out);
    for(int z=0; z<zmax+1;z++ ){
      for(int k=0;k<n+1;k++){
	//un[k]=u[k]+dt*a*(exp((c*x[k])/v)-1)+x[k]/c+2*j*dt;
	rhon[k]=(1-e(x[k],(j)*dt,dt))*rho[k]+1000*e(x[k],(j)*dt,dt);
      }

      /*
      if(j>39){
	 min=1000000;
      }


      if(j>59){
	 min=100000000;
	 }*/
      kk=remplissage(n,x,rhon,min);
     
      //kk=remplissage(n,x,un,min);
      //affichage_matrice(kk,n);
      xn=resol(kk,b,n);
      //affichage_vector(xn,n);

      if(j>0){
	rho=new_rho(xn,rhov,j,n,dt,c,v);
         }
      //affichage_vector(rho,n);
      
      xx=x;
      if(j==44){
	//	affichage_matrice(kk,n);
      }
      if(max_error_vector(x,xn,n)<epsilon){
	 cout<<"error_1= "<<max_error_vector(x,xn,n)<<" instant"<<j+1<<endl;
	cout<<"*************************"<<endl;
	break;
      }
     
      x=xn;
    }
    
    for(int o=0;o<n+1;o++){
      mon_flux << xx[o] << " " << rhon[o] << endl;
    }
    rho=rhon;
    //u=un;
    mon_flux.close();

  }
  //cout<<"la valeur et"<<rho.size()<<endl;
  


  return 0;
}
