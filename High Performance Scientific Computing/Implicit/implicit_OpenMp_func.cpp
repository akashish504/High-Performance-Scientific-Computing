#include <iostream>	
#include <time.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <omp.h>

using namespace std;
#define NUM_THREADS 4

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
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		final[i] = array[i]*scalar;
	}
}

void scalarmat(float array[sqn],float scalar,float final[sqn]);

void scalarmat(float array[sqn],float scalar,float final[sqn]){
	//static float new1[n];
	#pragma omp parallel for
	for(int i=0;i<n*n;i++){
		final[i] = array[i]*scalar;
	}
}

float arrayArray(float arr1[n], float arr2[n]);

float arrayArray(float arr1[n],float arr2[n]){
	float sum = 0;
	#pragma omp parallel for reduction(+:sum)
	for(int i =0;i<n;i++){
		sum = sum + arr1[i]*arr2[i];
	}
	return sum;
}



void matrixFlatArray(float arr1[],float arr2[],float final1[]);

void matrixFlatArray(float arr1[],float arr2[],float final1[]){
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
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		sum1[i] = a1[i] + a2[i] + a3[i];
	}
}

void addArr2(float a1[],float a2[],float sum1[]);

void addArr2(float a1[],float a2[],float sum1[]){
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		sum1[i] = a1[i] + a2[i] ;
	}
}
void addmat(float a1[],float a2[],float sum1[]);

void addmat(float a1[],float a2[],float sum1[]){
	#pragma omp parallel for
	for(int i=0;i<n*n;i++){
		sum1[i] = a1[i] + a2[i] ;
	}
}
void subArr2(float a1[],float a2[],float sub1[]);

void subArr2(float a1[],float a2[],float sub1[]){
	#pragma omp parallel for
	for(int i=0;i<n;i++){
		sub1[i] = a1[i] - a2[i] ;
	}
}

void inverse(float a[], float b[]);
void inverse(float a[], float b[])
{
	float temp=1.0;
	float inv[n+1][n+1],alpha[n+1], beta[n+2];
	float Tinv[n+1][n+1];
	#pragma omp parallel for
	for (int i = 1; i <= n; i++){
		#pragma omp parallel for
		for (int j = 1; j <= n; j++){
		inv[i][j] = a[(i-1)*n+(j-1)];
		}
	}
	alpha[0] = 1 ;

	for(int i=1;i<=n;i++){

		if(i==1) {alpha[i]=inv[i][i];}
		else {alpha[i]= inv[i][i]*alpha[i-1] - inv[i][i-1]*inv[i-1][i]*alpha[i-2];}
	}

     beta[n+1]=1.0;
	for(int i=n;i>=1;i--){
		if(i==n) beta[i]=inv[i][i];
		else beta[i]=inv[i][i]*beta[i+1] - inv[i+1][i]*inv[i][i+1]*beta[i+2];
	}
	

	Tinv[1][1]= 1/(inv[1][1] - (inv[2][1]*inv[1][2]*beta[3]/beta[2]));
	#pragma omp parallel for
	for(int i=2;i<=n;i++){
		if(i==n) Tinv[i][i] = 1/(inv[i][i] - (inv[i][i-1]*inv[i-1][i]*alpha[i-2]/alpha[i-1]));
		else Tinv[i][i] = 1/(inv[i][i] - (inv[i][i-1]*inv[i-1][i]*alpha[i-2]/alpha[i-1]) - (inv[i+1][i]*inv[i][i+1]*beta[i+2]/beta[i+1]));
	}

	for(int i=1;i<=n;i++){
		for(int j=1;j<=n;j++){
			if(i<j){
				#pragma omp parallel for reduction(*:temp)
				for(int k=1;k<=j-i;k++) {
					temp=temp*inv[j-k][j-k+1];
				}
				Tinv[i][j]=temp*alpha[i-1]*Tinv[j][j]*pow(-1,j-i)/alpha[j-1];
				temp=1.0;
			}
			if(i>j){
				#pragma omp parallel for reduction(*:temp)
				for(int k=1;k<=i-j;k++) {
					temp = temp*inv[j+k][j+k-1];
				}

		       	Tinv[i][j]=temp*beta[i+1]*Tinv[j][j]*pow(-1,i-j)/beta[j+1];	
				temp=1.0;
			}
			else continue;
		}
	}	

	#pragma omp parallel for
	for (int i = 1; i <= n; i++){
		#pragma omp parallel for
		for (int j = 1; j <= n; j++){
		b[(i-1)*n+(j-1)]= Tinv[i][j];
		}
	}

}


int main()
{	
	int start_s=clock();
	omp_set_num_threads(NUM_THREADS);
	ofstream outdata;
	//int n = 3;
	float Dt=0.01;
	int no_intervals=1000;
	//stiffnesses of individual springs
	float Kx[n+1],m[n];
	

	//intializing u,v,a,u_old,v_old,a_old, m and Kx
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

	//Populating M and K matrix#pragma omp parallel for
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


	for(int i=0;i<no_intervals;i++)
	{
		//block getting u_til_n+1;
		float u_til[n] = {0};
		float vt[n] = {0};
		float at[n] = {0};
		scalarArray(v_old,Dt,vt);
		scalarArray(a_old,(0.25*Dt*Dt),at);
		addArr3(u_old,vt,at,u_til);
		//getting v_til_n+1
		float at2[n] = {0};
		float v_til[n] = {0};
		scalarArray(a_old,0.5*Dt,at2);
		addArr2(v_old,at2,v_til);
		//getting a_n+1 in variable an
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
		float flatinv[sqn]= {0};
		float betaK[sqn]={0};
		float Tinv[sqn]={0};
		scalarmat(flatK,(0.25*Dt*Dt),betaK);
		addmat(flatM,betaK,flatinv);
		inverse(flatinv,Tinv);

		float t1[n] = {0};
		matrixFlatArray(flatK,u_til,t1);
		float t2[n]  = {0};
		subArr2(F,t1,t2);
		float an[n] = {0};
		matrixFlatArray(Tinv,t2,an);
		//getting u_n+1
		float newt[n] = {0};
		float newu[n] = {0};
		scalarArray(an,(Dt*0.5),newt);
		scalarArray(an,(0.25*Dt*Dt),newu);
		addArr2(v_til,newt,v);
		addArr2(u_til,newu,u);

		#pragma omp parallel for
		for(int j=0;j<n;j++){
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
	