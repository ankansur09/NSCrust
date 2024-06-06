#include <iostream>
#include <math.h>
#include <gsl/gsl_sf.h>
using namespace std;

int rNum = 40;
int thNum = 40;
double factor = 0.1;
int tNum = 10000000;
int plotSteps = 100;
double thtd=1.0;
double rmin=0.75;
int l=1;
double dt=1;
double const Pi=4*atan(1);



double init_A(double r, double th, double rmin)
{
double z0,z;
z0 = 0.041;
z = 1.0-r;
if (z<z0){
  return (1.0-(z*z/z0*z0))*sin(th)*sin(th) ;
  }
else{ return 0.0; }

}

double init_B(double r, double th, double rmin)
{
 return 0.0; 
}

double n(double r, double th){   return 1.0;}

// Function chi (1/(r^2 sin(th)*sin(th) * n(r,th))) at a given radius and theta

double chi(double r, double th){
    return 1.0/(pow(r*sin(th),2)*n(r,th));}

// Function eta returns the value of resistivity at a given radius and theta

double eta(double r, double th){return 1.0;}

/*void* solve_new_A()
{
for(int i=1;i<rNum+1;i++){
	for(int j=1;j<thNum-1;j++){
	r=rmin+(i-1)*dr;
	th=j*dth;
	if(i!=1){
		Aaux[i][j]=A[i][j]+dt*res_term_A[i][j]*gsA[i][j];
		Aaux[i][j]+= dt*((B[i][j+1]-B[i][j-1])*(A[i+1][j]-A[i-1][j]) -(B[i+1][j]-B[i-1][j])*(A[i][j+1]-A[i][j-1]))*hall_term_A[i][j];

			}
		}
	}
	return A;

}
void* solve_A_boundary(double legendre_comb, double legendre_comb_mid, double boundary_factors1, double boundary_factors2)
{
	
	for(int n=0;n<l;n++){
		a[n]=0;
		for(int j=0;j<thNum-1;j++){
			a[n]+=(A[rNum][j]+A[rNum][j+1])*legendre_comb_mid[n][j]*boundary_factors1[n];
		}
	}

        for(int j=1;j<thNum-1;j++){
		A[1][j]=0;
		A[0][j]=-A[2][j];
		A[rNum+1][j]=A[rNum-1][j];
		for(int n=0;n<l;n++){
			A[rNum+1][j]+=boundary_factors2[n]*a[n]*legendre_comb[n][j];
		}
	}
	return A;
}*/	



int main()
{    

double **B, **dBr, **dBth;
double **A, **Aaux, **gsA;
double a[l];

double dr,dth;
double r,th;
double res_term_A[rNum+2][thNum];
double hall_term_A[rNum+2][thNum];
//Arrays that contain the precalculated quantities neccesary to solve fluxes for B
double hall_rflux[rNum+2][thNum];
double hall_thflux[rNum+2][thNum];
double res_rflux[rNum+2][thNum];
double res_thflux[rNum+2][thNum];
double boundary_factors1[l], boundary_factors2[l];
double legendre_comb[l][thNum];
double legendre_comb_mid[l][thNum];

double sines[thNum], cotans[thNum];
dr=(1.0-rmin)/(rNum-1);
dth=Pi/(thNum-1);



A=new double*[rNum+2];
Aaux=new double*[rNum+2];
gsA=new double*[rNum+2];
//a=new double[l];

B=new double*[rNum+2];
dBr=new double*[rNum+2];
dBth=new double*[rNum+2];
for(int i=0;i<rNum+2;i++){
	B[i]=new double[thNum];
	dBr[i]=new double[thNum];
	dBth[i]=new double[thNum];
	A[i]=new double[thNum];
	Aaux[i]=new double[thNum];
	gsA[i]=new double[thNum];
	}

	
	
for(int j=1;j<thNum-1;j++){
		th=j*dth;
		sines[j]=sin(th);
		cotans[j]=cos(th)/sin(th);
		}

for(int i=0;i<rNum+2;i++){
	r=rmin+(i-1)*dr;
	for(int j=0;j<thNum-1;j++){
		th=j*dth;
		hall_term_A[i][j]= sines[j]*chi(r,th)/4/dr/dth;
		hall_rflux[i][j] = chi(r+dr/2,th)/8.0/dr/dth;
		hall_thflux[i][j]=-chi(r,th+dth/2)/8.0/dr/dth;
		res_term_A[i][j] = thtd*eta(r,th);
		res_rflux[i][j]  = thtd*eta(r+dr/2,th)/sin(th)/dr/dr;
		res_thflux[i][j] = thtd*eta(r,th+dth/2)/r/r/sin(th+dth/2.0)/dth/dth;
	}
}
	
	

for(int n=0;n<l;n++){
	boundary_factors1[n]=(n+1.0)/(n+2.0)*sqrt(Pi*(2*n+3))/2*dth;
	boundary_factors2[n]=-2*dr*(n+1)*sqrt((2*n+3)/(4*Pi));
	for(int j=0;j<thNum-1;j++){
		th=j*dth;
		legendre_comb[n][j]=(cos(th)*gsl_sf_legendre_Pl(n+1,cos(th))-gsl_sf_legendre_Pl(n,cos(th)));
		legendre_comb_mid[n][j]=(cos(th+dth*0.5)*gsl_sf_legendre_Pl(n+1,cos(th+dth*0.5))-gsl_sf_legendre_Pl(n,cos(th+dth*0.5)))/sin(th+dth*0.5);
		}
	}

	
	
for(int i=1;i<rNum+1;i++){
	r=rmin+(i-1)*dr;
	for(int j=1;j<thNum-1;j++){
		th=j*dth;
		B[i][j]=init_B(r,th,rmin);
		A[i][j]=init_A(r,th,rmin);
		}
	}
	

// fixes A boundary
for(int n=0;n<l;n++){
	a[n]=0;
	for(int j=0;j<thNum-1;j++){
		a[n]+=(A[rNum][j]+A[rNum][j+1])*legendre_comb_mid[n][j]*boundary_factors1[n];
	}
	}
for(int j=1;j<thNum-1;j++){
	A[1][j]=0;
	A[0][j]=-A[2][j];
	A[rNum+1][j]=A[rNum-1][j];
	for(int n=0;n<l;n++){
		A[rNum+1][j]+=boundary_factors2[n]*a[n]*legendre_comb[n][j];
	}
}
//end of A boundary

for(int j=1;j<thNum-1;j++){
	B[1][j]=0;
	B[2][j]=B[3][j]/2;
	B[rNum][j]=0;
	//B[rNum-1][j]=B[rNum-2][j]/2;
	B[rNum+1][j]=-B[rNum-1][j];
}

double t=0;

for(int i=1;i<rNum+1;i++){
		for(int j=1;j<thNum-1;j++){
			r=rmin+(i-1)*dr;
			gsA[i][j]=(A[i+1][j]+A[i-1][j]-2*A[i][j])/dr/dr+1/r/r*(A[i][j+1]+A[i][j-1]-2*A[i][j])/dth/dth-1/r/r*cotans[j]*(A[i][j+1]-A[i][j-1])/2/dth;
		}
	}


/*
for(int k=0;k<=tNum;k++){


	double newdt=100000000;
	double localnewdt=100000000;
	double temp=0;
	for(int i=1;i<rNum+1;i++){
		for(int j=1;j<thNum-1;j++){
			r=rmin+(i-1)*dr;
			gsA[i][j]=(A[i+1][j]+A[i-1][j]-2*A[i][j])/dr/dr+1/r/r*(A[i][j+1]+A[i][j-1]-2*A[i][j])/dth/dth-1/r/r*cotans[j]*(A[i][j+1]-A[i][j-1])/2/dth;
		}
	}
	for(int i=2;i<rNum+1;i++){
		for(int j=1;j<thNum-1;j++){
			r=rmin+(i-1)*dr;
			th=j*dth;
			temp=dr/sqrt(pow(gsA[i][j]/(r*sines[j]),2)+pow((B[i+1][j]-B[i-1][j])/(2*dr*r*sines[j]),2)+pow((B[i][j+1]-B[i][j-1])/(2*dth*r*r*sines[j]),2));
			if(temp<localnewdt){
				localnewdt=temp;
			}
		}
		if(localnewdt<newdt){
			newdt=localnewdt;
		}
	}
	//If Ohmic critical step is smaller, use that one
	if(newdt>dr*dr/thtd){
		newdt=dr*dr/thtd;
	}
	dt=factor*newdt;
	if(dt<0.0000000001){
		dt=0.0000000001;
	}
	
	solve_new_A();
	solve_new_B();
	
		
	for(int i=1;i<rNum+1;i++){
		for(int j=1;j<thNum-1;j++){
			if(isinf(dBr[i][j])||isinf(dBth[i][j])||isnan(dBr[i][j])||isnan(dBth[i][j])||isinf(Aaux[i][j])||isnan(Aaux[i][j]))
				{
				break;
				}
			}
			B[i][j]+=(dBr[i][j]-dBr[i-1][j]+dBth[i][j]-dBth[i][j-1])*sines[j];
			
			A[i][j]=Aaux[i][j];
		}
	
	solve_B_boundary();
	solve_A_boundary();
	t+=dt;

	}
	
*/
//cout << A[41][10]<<endl;
//cout << B[41][10]<<endl;
//cout << a[1]<<endl;
//cout<< gsA[36][10];

cout << gsl_sf_legendre_Pl(2,1);
cout << gsl_sf_legendre_Pl(2,0.5);
cout << gsl_sf_legendre_Pl(2,0.8);

return 0;

}


