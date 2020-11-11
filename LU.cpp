#include<vector>
#include<cmath>
#include<iostream>
#include "LU.h"

using namespace std;

vector<double> resol(vector<vector<double>> a, vector<double> b, int n){
  vector<double> x(n+1,0.),y(n+1,0.),z(n+1,0.);
  y[0]=a[0][0];
  z[0]=b[0]/y[0];

  for(int i=1;i<n+1;i++){
    y[i]=a[i][i]-(a[i][i-1]*a[i-1][i])/(y[i-1]*1.);
    z[i]=(b[i]-a[i-1][i]*z[i-1])/(y[i]*1.);
  }

  x[n]=z[n];

  for(int j=1;j<n+1;j++){
    x[n-j]=z[n-j]-(a[n-j][n-j+1]*x[n-j+1])/(y[n-j]*1.);
  }
  
  return x;
}

double max_error_vector(vector<double> a,vector<double> b,int n){
  double max_error=0;
  for(int i=0;i<n+1;i++){
    if(abs(a[i]-b[i])>max_error){
      max_error=abs(a[i]-b[i]);
    }
  }
  return max_error;
}
