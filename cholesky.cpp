#include<vector>
#include<cmath>
#include "mod_fonction.h"

using namespace std;

vector<vector<double>> chol(vector<vector<double>> a,int n)
{

  vector<vector<double>> l(n+1,vector<double>(n+1,0.));
  double c(0.),cc(0.);
  int k(0);

  for(int j=0;j<n+1;j++){
    k=0;
    if(k<=j-1){
      c=0.;
      for(k=0;k<=j-1;k++){
	c=c+l[j][k]*l[j][k];
      }
      l[j][j]=sqrt(a[j][j]-c);
    }
    else{
      l[j][j]=sqrt(a[j][j]);
    }

    for(int i=j+1;i<n+1;i++)
      {
	k=0;
	if(k<=j-1){
	  cc=0;
	  for(k=0;k<=j-1;k++){
	    cc=cc+l[i][k]*l[j][k];
	  }
	  l[i][j]=(a[i][j]-cc)/l[j][j];
	}else{
	  l[i][j]=a[i][j]/l[j][j];
	}
      }

  }
  return l;
}



vector<vector<double>> transpose(vector<vector<double>> a,int n)
{
  vector<vector<double>> at(n+1,vector<double>(n+1,0.));
  for(int i=0;i<n+1;i++){
    for(int j=0;j<n+1;j++){
      at[i][j]=a[j][i];
    }
  }
  return at;
}




vector<double> rechol(vector<vector<double>> l, vector<double> b, int n)
{
  vector<double> x(n+1,0.),y(n+1,0.);
  vector<vector<double>> lt(n+1,vector<double>(n+1,0.));
  double c=0.;
  y[0]=b[0]/l[0][0];
  for(int i=1;i<n+1;i++){
    c=0.;
    for(int k=0;k<=i-1;k++){
      c=c+l[i][k]*y[k];
    }
    y[i]=(b[i]-c)/l[i][i];
  }
  lt=transpose(l,n);

  x[n]=y[n]/lt[n][n];

  for(int j=n-1;j>=0;j--){
    c=0.;
    for(int z=j+1;z<n+1;z++){
      c=c+x[z]*lt[j][z];
    }
    x[j]=(y[j]-c)/lt[j][j];
  }
  return x;
} 
  

			







