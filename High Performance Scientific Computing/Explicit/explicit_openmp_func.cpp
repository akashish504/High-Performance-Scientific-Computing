#include <iostream>	
#include <time.h>
#include <fstream>
#include <iomanip>
#include <omp.h>

using namespace std;
#define NUM_THREADS 4
int start_s=clock();
const int n =100;
const int sqn = n*n;

void display(float array[n]);

void display(float array[n]){
	for(int i = 0;i<n;i++){
			cout<<array[i]<<endl;
	}
}

void scalarArray(float array[n],float scalar,float final[n]);

void scalarArray(float array[n],float scalar,float final[n]){
	//static float new1[n];
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		final[i] = array[i]*scalar;
	}
}

float arrayArray(float arr1[n], float arr2[n]);

float arrayArray(float arr1[n],float arr2[n]){
	float sum = 0; 
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel for reduction(+:sum)
	for(int i =0;i<n;i++){
		sum = sum + arr1[i]*arr2[i];
	}
	return sum;
}



void matrixFlatArray(float arr1[],float arr2[],float final1[]);

void matrixFlatArray(float arr1[],float arr2[],float final1[]){
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		float temp[n] = {0};
	    #pragma omp parallel for
		for(int j=0;j<n;j++){
			temp[j] = arr1[i*n + j];
		}
		final1[i] = arrayArray(temp,arr2);
	}
}

void addArr3(float a1[],float a2[],float a3[],float sum1[]);

void addArr3(float a1[],float a2[],float a3[],float sum1[]){
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		sum1[i] = a1[i] + a2[i] + a3[i];
	}
}

void addArr2(float a1[],float a2[],float sum1[]);

void addArr2(float a1[],float a2[],float sum1[]){
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		sum1[i] = a1[i] + a2[i];
	}
}

void subArr2(float a1[],float a2[],float sub1[]);

void subArr2(float a1[],float a2[],float sub1[]){
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		sub1[i] = a1[i] - a2[i] ;
	}
}





int main()
{	
	int start_s=clock();
	ofstream outdata;
	//int n = 3;
	float Dt=0.01;
	int no_intervals=1000;
	//stiffnesses of individual springs
	float Kx[n+1],m[n];
	

	//intializing u,v,a,u_old,v_old,a_old
	float u[n],v[n],a[n],F[n]; 
	float u_old[n],v_old[n],a_old[n];
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		u[i]=0;
		if(i==49){
			v[i]=2;
		}
		else{
			v[i]=0;
		}
		a[i]=0;
		F[i]=0;
		

		
		u_old[i]=0;
		v_old[i]=v[i];
		a_old[i]=0;
		m[i]=1;
		Kx[i]=1000;
	}
	Kx[n]=1000;

	outdata.open("u.csv",ios::app);
	//K and M matrices
	float M[n][n],K[n][n];

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

	
	float flatK[sqn] = {0};
	float flatM[sqn] = {0};
	#pragma omp parallel for
	for(int j=0;j<n;j++){
		#pragma omp parallel for
		for(int k=0;k<n;k++){
			flatM[j*n+k]=M[j][k];
			flatK[j*n+k] = K[j][k];
		}
	}

	for(int i=0;i<no_intervals;i++){
		//block getting u_til_n+1;
		float u_til[n] = {0};
		float vt[n] = {0};
		float at[n] = {0};
		scalarArray(v_old,Dt,vt);
		scalarArray(a_old,(0.5*Dt*Dt),at);
		addArr3(u_old,vt,at,u_til);
		//getting v_til_n+1
		float at2[n] = {0};
		float v_til[n] = {0};
		scalarArray(a_old,0.5*Dt,at2);
		addArr2(v_old,at2,v_til);
		//getting a_n+1 in variable an
		float t1[n] = {0};
		matrixFlatArray(flatK,u_til,t1);
		float t2[n]  = {0};
		subArr2(F,t1,t2);
		float an[n] = {0};
		matrixFlatArray(flatM,t2,an);
		//getting u_n+1
		float newt[n] = {0};
		scalarArray(an,(Dt*Dt*0.5),newt);
		addArr2(v_til,newt,v);

		#pragma omp parallel for
		for(int j=0;j<n;j++){
			u[j] = u_til[j];
			a[j] = an[j];
			u_old[j] = u[j];
			a_old[j] = a[j];
			v_old[j] = v[j];
		}

		outdata<<u[48]<<','<<u[49]<<','<<u[50]<<endl;

	

}
int stop_s=clock();
cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << endl;
}

//Print function is too tedious to be made.
	