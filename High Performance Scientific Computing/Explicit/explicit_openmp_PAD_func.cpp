#include <iostream>	
#include <time.h>
#include <fstream>
#include <iomanip>
#include <omp.h>

using namespace std;
#define NUM_THREADS 4
// run following command on your terminal
//getconf LEVEL1_DCACHE_LINESIZE
//you will get cache line size for your system 
//if your pc has 64 bit architecture then divide cache line size by 8 to get PAD number
//If your pc has 32 bit architecture then divide cache line size by 4 to get PAD number
#define PAD 8
const int n =100;
const int sqn = n*n;
//in this code I have increased all array size such that false sharing would not happen between 2 threads
//to increase array size and to ensure that no 2 data elements lies on same cache line I have increase dimension of each array
void display(float array[n][PAD]);

void display(float array[n][PAD]){
	for(int i = 0;i<n;i++){
			cout<<array[i][0]<<endl;
	}
}

void scalarArray(float array[n][PAD],float scalar,float final[n][PAD]);

void scalarArray(float array[n][PAD],float scalar,float final[n][PAD]){
	//static float new1[n];
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		final[i][0] = array[i][0]*scalar;
	}
}

float arrayArray(float arr1[n][PAD], float arr2[n][PAD]);

float arrayArray(float arr1[n][PAD],float arr2[n][PAD]){
	float sum = 0; 
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel for reduction(+:sum)
	for(int i =0;i<n;i++){
		sum = sum + arr1[i][0]*arr2[i][0];
	}
	return sum;
}



void matrixFlatArray(float arr1[][PAD],float arr2[][PAD],float final1[][PAD]);

void matrixFlatArray(float arr1[][PAD],float arr2[][PAD],float final1[][PAD]){
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		float temp[n][PAD] = {0};
	    #pragma omp parallel for
		for(int j=0;j<n;j++){
			temp[j][0] = arr1[i*n + j][0];
		}
		final1[i][0] = arrayArray(temp,arr2);
	}
}

void addArr3(float a1[][PAD],float a2[][PAD],float a3[][PAD],float sum1[][PAD]);

void addArr3(float a1[][PAD],float a2[][PAD],float a3[][PAD],float sum1[][PAD]){
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		sum1[i][0] = a1[i][0] + a2[i][0] + a3[i][0];
	}
}

void addArr2(float a1[][PAD],float a2[][PAD],float sum1[][PAD]);

void addArr2(float a1[][PAD],float a2[][PAD],float sum1[][PAD]){
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		sum1[i][0] = a1[i][0] + a2[i][0] ;
	}
}

void subArr2(float a1[][PAD],float a2[][PAD],float sub1[][PAD]);

void subArr2(float a1[][PAD],float a2[][PAD],float sub1[][PAD]){
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		sub1[i][0] = a1[i][0] - a2[i][0] ;
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
	float Kx[n+1][PAD],m[n][PAD];
	

	//intializing u,v,a,u_old,v_old,a_old
	float u[n][PAD],v[n][PAD],a[n][PAD],F[n][PAD]; 
	float u_old[n][PAD],v_old[n][PAD],a_old[n][PAD];
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		u[i][0]=0;
		if(i==49){
			v[i][0]=2;
		}
		else{
			v[i][0]=0;
		}
		
		a[i][0]=0;
		F[i][0]=0;
				
		u_old[i][0]=0;
		v_old[i][0]=v[i][0];
		a_old[i][0]=0;
		m[i][0]=1;
		Kx[i][0]=1000;
	}

	Kx[n][0]=1000;

	outdata.open("new4.csv",ios::app);
	//K and M matrices
	float M[n][n][PAD],K[n][n][PAD];

	//Populating M and K matrix
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		#pragma omp parallel for
		for(int j=0;j<n;j++){
			M[i][j][0]=0;
			K[i][j][0]=0;
			if(i==j){
				M[j][j][0]=m[i][0];
				K[j][j][0]=Kx[i][0]+Kx[i+1][0];
			}
			if(j==i-1){
				K[i][j][0]=-Kx[i][0];
			}
			if(j==i+1){
				K[i][j][0]=-Kx[i+1][0];
			}
		}
	}

	
	float flatK[sqn][PAD] = {0};
	float flatM[sqn][PAD] = {0};
	#pragma omp parallel for
	for(int j=0;j<n;j++){
		#pragma omp parallel for
		for(int k=0;k<n;k++){
			flatM[j*n+k][0]=M[j][k][0];
			flatK[j*n+k][0] = K[j][k][0];
		}
	}

	for(int i=0;i<no_intervals;i++){
		//block getting u_til_n+1;
		float u_til[n][PAD] = {0};
		float vt[n][PAD] = {0};
		float at[n][PAD] = {0};
		scalarArray(v_old,Dt,vt);
		scalarArray(a_old,(0.5*Dt*Dt),at);
		addArr3(u_old,vt,at,u_til);
		//getting v_til_n+1
		float at2[n][PAD] = {0};
		float v_til[n][PAD] = {0};
		scalarArray(a_old,0.5*Dt,at2);
		addArr2(v_old,at2,v_til);
		//getting a_n+1 in variable an
		float t1[n][PAD] = {0};
		matrixFlatArray(flatK,u_til,t1);
		float t2[n][PAD]  = {0};
		subArr2(F,t1,t2);
		float an[n][PAD] = {0};
		matrixFlatArray(flatM,t2,an);
		//getting u_n+1
		float newt[n][PAD] = {0};
		scalarArray(an,(Dt*Dt*0.5),newt);
		addArr2(v_til,newt,v);
		#pragma omp parallel for
		for(int j=0;j<n;j++){
			u[j][0] = u_til[j][0];
			a[j][0] = an[j][0];
		}
		#pragma omp parallel for
		for(int j=0;j<n;j++){
			u_old[j][0] = u[j][0];
			a_old[j][0] = a[j][0];
			v_old[j][0] = v[j][0];
		}
		outdata<<u[0][0]<<endl;

	

}
int stop_s=clock();
cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << endl;
}

//Print function is too tedious to be made.
	