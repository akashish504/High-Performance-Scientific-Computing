#include <iostream>	
#include <time.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <omp.h>

using namespace std;
#define NUM_THREADS 4

int main(  int argc, char *argv[])
{	

	ofstream outdata;
	omp_set_num_threads(NUM_THREADS);
	clock_t tStart = clock();
	int my_PE_num , n1;
	//ofstream outdata;
	int n = 100;
	double Dt=0.01;
	int no_intervals=1000;
	//stiffnesses of individual springs
	double Kx[n+1],m[n];
	
	//initializing Kx and m matrix
	#pragma omp parallel for
	for(int i=0;i<n+1;i++){
		m[i]=1;
		Kx[i]=1000;
	} 
	Kx[n]=1000;


	outdata.open("u.csv",ios::app);
	//K and M matrices
	double M[n][n],K[n][n];
	
	//Populating M and K matrix
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		#pragma omp parallel for
		for(int j=0;j<n;j++){
			M[i][j]=0;
			K[i][j]=0;
			if(i==j){
				M[j][j]=m[i];
				K[j][j]=Kx[i]+Kx[i+1];
			}
			if(j==i-1){
				K[i][j]=-Kx[i];
			}
			if(j==i+1){
				K[i][j]=-Kx[i+1];
			}
		}
	}

	//Initializing MPI 
	double u[n] = {0},u_old[n] = {0};
	double v[n] = {0},v_old[n] = {0};
	double a[n] = {0},a_old[n] = {0};
	double F[n] = {0};
	v_old[49] =2;

	for(int k=0;k<no_intervals;k++){
		//Predictors for u and v.
		#pragma omp parallel for
		for(int i=0;i<n;i++){
			u[i]=u_old[i]+Dt*v_old[i]+Dt*Dt*a_old[i]/2;
			v[i]=v_old[i]+Dt*a_old[i]/2;
			u_old[i]=u[i];
		}
    
		#pragma omp parallel for
		for(int i=0;i<n;i++){
			if(i==0){
				a[i]=(F[i]-K[i][i]*u[i]-K[i][i+1]*u[i+1])/m[i];
			}
			else if(i==n-1){
				a[i]=(F[i]-K[i][i]*u[i]-K[i][i-1]*u[i-1])/m[i];
			}			
			else{
				a[i]=(F[i]-K[i][i]*u[i]-K[i][i-1]*u[i-1]-K[i][i+1]*u[i+1])/m[i];
			}
			a_old[i]=a[i];
		}


		//Corrector for u and v
		#pragma omp parallel for
		for(int i=0;i<n;i++){
			v[i]=v[i]+a[i]*Dt/2;
			v_old[i]=v[i];
		}
		outdata<<u[48]<<','<<u[49]<<','<<u[50]<<endl;
	}
	cout<<("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC)<<endl;
}


