#include <iostream>	
#include <time.h>
#include <fstream>
#include <iomanip>
#include <time.h>

using namespace std;

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
	for(int i=0;i<n;i++){
		final[i] = array[i]*scalar;
	}
}

float arrayArray(float arr1[n], float arr2[n]);

float arrayArray(float arr1[n],float arr2[n]){
	float sum = 0;
	for(int i =0;i<n;i++){
		sum = sum + arr1[i]*arr2[i];
	}
	return sum;
}



void matrixFlatArray(float arr1[],float arr2[],float final1[]);

void matrixFlatArray(float arr1[],float arr2[],float final1[]){
	for(int i=0;i<n;i++){
		float temp[n] = {0};
		for(int j=0;j<n;j++){
			temp[j] = arr1[i*n + j];
		}
		final1[i] = arrayArray(temp,arr2);
	}
}

void addArr3(float a1[],float a2[],float a3[],float sum1[]);

void addArr3(float a1[],float a2[],float a3[],float sum1[]){
	for(int i=0;i<n;i++){
		sum1[i] = a1[i] + a2[i] + a3[i];
	}
}

void addArr2(float a1[],float a2[],float sum1[]);

void addArr2(float a1[],float a2[],float sum1[]){
	for(int i=0;i<n;i++){
		sum1[i] = a1[i] + a2[i] ;
	}
}

void subArr2(float a1[],float a2[],float sub1[]);

void subArr2(float a1[],float a2[],float sub1[]){
	for(int i=0;i<n;i++){
		sub1[i] = a1[i] - a2[i] ;
	}
}





int main()
{	

	clock_t tStart = clock();
	ofstream outdata;
	float Dt=0.01;
	int no_intervals=1000;
	//stiffnesses of individual springs
	float Kx[n+1],m[n];
	
	//Initailizing Kx matrix
	for(int i=0;i<n+1;i++){
		Kx[i]=1;
	} 
	//initializing m matrix
	for(int i=0;i<n+1;i++){
		m[i]=1;
	} 

	//intializing u,v,a,u_old,v_old,a_old
	float u[n],v[n],a[n],F[n]; 
	float u_old[n],v_old[n],a_old[n];
	for(int i=0;i<n;i++){
		u[i]=0;
		if(i==0){
			v[i]=1;
		}
		else{
			v[i]=0;
		}
		a[i]=0;
		F[i]=0;
		u_old[i]=0;
		v_old[i]=v[i];
		a_old[i]=0;
	}

//cout<<u[0]<<"___"<<u[1]<<"___"<<u[2]<<endl;

	outdata.open("new3.csv",ios::app);
	//K and M matrices
	float M[n][n],K[n][n];

	//Populating M and K matrix
	for(int i=0;i<n;i++){
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
cout<<u[0]<<"-----"<<u[1]<<"___"<<u[2]<<endl;

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
		float flatK[sqn] = {0};
		float flatM[sqn] = {0};
		for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				flatM[j*n+k]=M[j][k];
				flatK[j*n+k] = K[j][k];
			}
		}
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

		cout<<"--------Displaying u at end of each----"<<endl;
		display(u);
		for(int j=0;j<n;j++){
			u[j] = u_til[j];
			a[j] = an[j];
		}
		for(int j=0;j<n;j++){
			u_old[j] = u[j];
			a_old[j] = a[j];
			v_old[j] = v[j];
		}
		outdata<<u[0]<<","<<u[1]<<","<<u[2]<<endl;
		//cout<<u[0]<<","<<u[1]<<","<<u[2]<<endl;

}

	//cout<<("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC)<<endl;
}

//Print function is too tedious to be made.
	