#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <iostream>
#include <time.h>


using namespace std;

const int n = 3;
const int sqn = n*n;
const float beta = 0.5;





__global__ void add(int a, int b, int *c)
{
 *c = a + b;
 printf("ss\n");
}

__global__ void hello(double *u,double *v, double *a,double *F,double *flatK, double *flatM , double *u_old,double*v_old,double *a_old,int size,double Dt,double *new1,double *test){
	int i = threadIdx.x;
	double sum = 0;
	
	//printf("this thing %f\n",new1[i] );
	//printf("sup\n");
	//printf("%f\n",Dt );
	

	
	for(int j = 0;j<1000;j++){
		//printf("(%f)\n",v[i] );
		u[i] = u_old[i] + Dt*v_old[i] + Dt*Dt*0.5*a_old[i];
		//printf("%f\n",u[i] );
		__syncthreads();
		if(i == 0){
			a[i] = (F[i] - flatK[i*size + i]*u[i] - flatK[i*size + i+1]*u[i+1]);
			
		} 
		else if(i == size-1){
			a[i] = (F[i] - flatK[i*size + i]*u[i] - flatK[i*n + (i-1)]*u[i-1]);
			
		}
		else{
			a[i] = (F[i] - flatK[i*size + i]*u[i] - flatK[i*size + (i-1)]*u[i-1] - flatK[i*size + (i + 1)]*u[i+1]);
			
		}
		__syncthreads();

		for(int k=0;k<size;k++){
			sum = 0;
			for(int p = 0;p<size;p++){
				sum = sum + (new1[k*size + p] * a[i]);
				__syncthreads();
			}
			test[i] = sum;
			__syncthreads();
		}
		a[i] = test[i];
		__syncthreads();


		v[i] = v_old[i] + Dt*(a[i] + a_old[i])*0.5;
		__syncthreads();
		u[i] = u[i] +  Dt*Dt*0.5*(a[i]);
		__syncthreads();
		v_old[i] = v[i];
		__syncthreads();
		u_old[i] = u[i];
		__syncthreads();
		a_old[i] = a[i];
		__syncthreads();
		
	if(i==0)	
	printf("%f\n",u[i]);

	
	//*var1 = u[0];
}
	//printf("%f\n",v[i] );
	
	
}




void inverse(double a[], double b[])
{
	double temp=1.0;
	double inv[n+1][n+1],alpha[n+1], beta[n+2];
	double Tinv[n+1][n+1];
	for (int i = 1; i <= n; i++){
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
	//cout<<endl;

	Tinv[1][1]= 1/(inv[1][1] - (inv[2][1]*inv[1][2]*beta[3]/beta[2]));
	for(int i=2;i<=n;i++){
		if(i==n) Tinv[i][i] = 1/(inv[i][i] - (inv[i][i-1]*inv[i-1][i]*alpha[i-2]/alpha[i-1]));
		else Tinv[i][i] = 1/(inv[i][i] - (inv[i][i-1]*inv[i-1][i]*alpha[i-2]/alpha[i-1]) - (inv[i+1][i]*inv[i][i+1]*beta[i+2]/beta[i+1]));
	}

	for(int i=1;i<=n;i++){
		for(int j=1;j<=n;j++){
			if(i<j){

				for(int k=1;k<=j-i;k++) {
					temp=temp*inv[j-k][j-k+1];
				}
				Tinv[i][j]=temp*alpha[i-1]*Tinv[j][j]*pow(-1,j-i)/alpha[j-1];
				temp=1.0;
			}
			if(i>j){
				for(int k=1;k<=i-j;k++) {
					temp = temp*inv[j+k][j+k-1];
				}

		       	Tinv[i][j]=temp*beta[i+1]*Tinv[j][j]*pow(-1,i-j)/beta[j+1];	
				temp=1.0;
			}
			else continue;
		}
	}	

	for (int i = 1; i <= n; i++){
		for (int j = 1; j <= n; j++){
			b[(i-1)*n+(j-1)]= Tinv[i][j];
		}
	}

}


int main(){

	
	//cudaDeviceReset();

	double Dt = 0.01;
	int no_intervals = 1000;
	//Initialiing Stiffness Array and Mass Array
	double Kx[n+1],m[n];

	for(int i=0;i<n+1;i++){
		Kx[i] = 1000;
		if(i != n){
			m[i] = 1;
		}
	}
	

	

	double *u,*v,*a,*F,*flatK,*flatM,*u_old,*v_old,*a_old,*new1,*test;
	u = (double *)malloc(n*sizeof(double));
	v = (double *)malloc(n*sizeof(double));
	a = (double *)malloc(n*sizeof(double));
	F = (double *)malloc(n*sizeof(double));
	flatM = (double *)malloc(n*n*sizeof(double));
	flatK = (double *)malloc(n*n*sizeof(double));
	u_old = (double *)malloc(n*sizeof(double));
	v_old = (double *)malloc(n*sizeof(double));
	a_old = (double *)malloc(n*sizeof(double));
	new1 = (double *)malloc(n*n*sizeof(double));
	test = (double *)malloc(n*sizeof(double));



	for(int i=0;i<n;i++){
		u[i] = 0;
		if(i ==0){
			v[i] = 2;
		}
		else{
			v[i] = 0;
		}
		a[i] = 0;
		F[i] = 0;
		u_old[i] = u[i];
		v_old[i] = v[i];
		a_old[i] = a[i];
	}

	//Mass and K Matrices
	double M[n][n],K[n][n];
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			M[i][j]=0;
			K[i][j]=0;
			if(i==j){
				M[i][j]=m[i];
				K[i][j]=Kx[i]+Kx[i+1];
			}
			if(j==i-1){
				K[i][j]=-Kx[i];
			}
			if(j==i+1){
				K[i][j]=-Kx[i+1];
			}
		}
	}

	for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				flatM[j*n+k]=M[j][k];
				flatK[j*n+k] = K[j][k];
			}
	}
	
	int size = n;
	double toInv[sqn];

	for(int i =0;i<sqn;i++){
		toInv[i] = beta*Dt*Dt*flatK[i] + flatM[i];
	}


	double new2[sqn];
	inverse(toInv,new2);
	for(int i=0;i<sqn;i++){
		//cout<<toInv[i]<<endl;
	}
	for(int i=0;i<sqn;i++){
		new1[i] = new2[i];
	}


	//Making Device Copies of all matrices (u,v,a,u_old,v_old,a_old,F,flatM,flatK);
	double *dev_u,*dev_v,*dev_a,*dev_F,*dev_flatM,*dev_flatK,*dev_u_old,*dev_v_old,*dev_a_old,*dev_size,*dev_new1,*dev_test;
	cudaMalloc((void **)&dev_u,n*sizeof(double));
	cudaMalloc((void **)&dev_v,n*sizeof(double));
	cudaMalloc((void **)&dev_a,n*sizeof(double));
	cudaMalloc((void **)&dev_F,n*sizeof(double));
	cudaMalloc((void **)&dev_flatK,n*n*sizeof(double));
	cudaMalloc((void **)&dev_flatM,n*n*sizeof(double));
	cudaMalloc((void **)&dev_u_old,n*n*sizeof(double));
	cudaMalloc((void **)&dev_v_old,n*n*sizeof(double));
	cudaMalloc((void **)&dev_a_old,n*n*sizeof(double));
	cudaMalloc((void **)&dev_size,sizeof(int));
	cudaMalloc((void **)&dev_new1,n*n*sizeof(double));
	cudaMalloc((void **)&dev_test,n*sizeof(double));


	//Transfering all arrays to Device

	cudaMemcpy(dev_u,u,n*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v,v,n*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_a,a,n*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_F,F,n*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_flatK,flatK,n*n*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_flatM,flatM,n*n*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_u_old,u_old,n*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v_old,v_old,n*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_a_old,a_old,n*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_new1,new1,n*n*sizeof(double),cudaMemcpyHostToDevice);
	//cudaMemcpy(dev_size,size,sizeof(int),cudaMemcpyHostToDevice);

	//Running Kernal
	int N = 1; //Number of Blocks
	int t = n; // threads per block
	//Constraint -- N*t = n;
	int a1,b,c;
	int *dev_c;
	a1=3;	
	b=4;
	cudaMalloc((void**)&dev_c, sizeof(int));
	//hello<<<1,10>>>();
	//add<<<1,10>>>(a1,b,dev_c);
	hello<<<1,n>>>(dev_u,dev_v,dev_a,dev_F,dev_flatK,dev_flatM,dev_u_old,dev_v_old,dev_a_old,size,Dt,dev_new1,dev_test);
	//cudaMemcpy(&c, dev_c, sizeof(int), cudaMemcpyDeviceToHost);
	//printf("%d + %d is %d\n", a1, b, c);
	cudaFree(dev_c);
	cudaFree(dev_u);
	cudaFree(dev_v);
	cudaFree(dev_a);
	cudaFree(dev_flatK);
	cudaFree(dev_flatM);
	cudaFree(dev_u_old);
	cudaFree(dev_v_old);
	cudaFree(dev_a_old);

	//printf("djdm\n" );
	

	//func<<<1,512>>>(dev_u,dev_v,dev_a,dev_F,dev_flatK,dev_flatM,dev_u_old,dev_v_old,dev_a_old);
	


free(u);
free(v);
free(a);
free(F);
free(flatM);
free(flatK);
free(u_old);
free(v_old);
free(a_old);
free(new1);
free(test);









return 0;
}

