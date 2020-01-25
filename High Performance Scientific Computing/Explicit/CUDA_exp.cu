#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

#include <time.h>


using namespace std;

const int n = 100;    // no of threads and DoF
const int sqn = n*n;



__global__ void add(int a, int b, int *c)
{
 *c = a + b;
 printf("ss\n");
}

__global__ void hello(double *u,double *v, double *a,double *F,double *flatK, double *flatM , double *u_old,double*v_old,double *a_old,int size,double Dt){
	int i = threadIdx.x;
	//printf("sup\n");
	//printf("%f\n",Dt );

	for(int j = 0;j<1000;j++){   // 1000 = no_intervals
		//printf("(%f)\n",v[i] );
		u[i] = u_old[i] + Dt*v_old[i] + Dt*Dt*0.5f*a_old[i];
		//printf("%f\n",u[i] );
		__syncthreads();
		if(i == 0){
			a[i] = (F[i] - flatK[i*size + i]*u[i] - flatK[i*size + i+1]*u[i+1])/flatM[i*size + i];
			
		} 
		else if(i == size-1){
			a[i] = (F[i] - flatK[i*size + i]*u[i] - flatK[i*n + (i-1)]*u[i-1])/flatM[i*size + i];
			
		}
		else{
			a[i] = (F[i] - flatK[i*size + i]*u[i] - flatK[i*size + (i-1)]*u[i-1] - flatK[i*size + (i + 1)]*u[i+1])/flatM[i*size + i];
			
		}
		__syncthreads();

		v[i] = v_old[i] + Dt*(a[i] + a_old[i])*0.5f;
		__syncthreads();
		v_old[i] = v[i];
		__syncthreads();
		u_old[i] = u[i];
		__syncthreads();
		a_old[i] = a[i];
		__syncthreads();
		
	if(i==48)						///////////   49, 50, 51
	printf("%f\n",u[i]);
	
	//*var1 = u[0];
}
	//printf("%f\n",v[i] );
	
}
int main(){
	
//cudaDeviceReset();
	clock_t time_i,time_f;
	time_i = clock();

	double Dt = 0.01f;
	int no_intervals = 1000;       /////////////////////////// dof
	//Initialiing Stiffness Array and Mass Array
	double Kx[n+1],m[n];

	for(int i=0;i<n+1;i++){
		Kx[i] = 1000;     ////////////////////////////////////    SPRING CONST
		if(i != n){
			m[i] = 1;	///////////////////////// MASS
		}
	}
	

	double *u,*v,*a,*F,*flatK,*flatM,*u_old,*v_old,*a_old;
	u = (double *)malloc(n*sizeof(double));
	v = (double *)malloc(n*sizeof(double));
	a = (double *)malloc(n*sizeof(double));
	F = (double *)malloc(n*sizeof(double));
	flatM = (double *)malloc(n*n*sizeof(double));
	flatK = (double *)malloc(n*n*sizeof(double));
	u_old = (double *)malloc(n*sizeof(double));
	v_old = (double *)malloc(n*sizeof(double));
	a_old = (double *)malloc(n*sizeof(double));


	for(int i=0;i<n;i++){
		u[i] = 0;
		if(i ==49){
			v[i] = 2;      /////////////    velocity
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

	//Making Device Copies of all matrices (u,v,a,u_old,v_old,a_old,F,flatM,flatK);
	double *dev_u,*dev_v,*dev_a,*dev_F,*dev_flatM,*dev_flatK,*dev_u_old,*dev_v_old,*dev_a_old,*dev_size;
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
	//for(int ii=0;ii<1000;ii++) 
	hello<<<1,n>>>(dev_u,dev_v,dev_a,dev_F,dev_flatK,dev_flatM,dev_u_old,dev_v_old,dev_a_old,size,Dt);
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

	
	time_f = clock();
	double exTime = time_f - time_i;
	printf("Time Taken: %f\n",exTime/CLOCKS_PER_SEC);
	

	//func<<<1,512>>>(dev_u,dev_v,dev_a,dev_F,dev_flatK,dev_flatM,dev_u_old,dev_v_old,dev_a_old);
















return 0;
}

