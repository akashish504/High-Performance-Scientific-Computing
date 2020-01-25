

#include <iostream>	
#include <time.h>
#include <iomanip>
#include <time.h>
#include "mpi.h"

using namespace std;


int main(  int argc, char *argv[])
{	
	clock_t tStart = clock();
	int my_PE_num , n1;
	//ofstream outdata;
	int n = 3;
	double Dt=0.01;
	int no_intervals=10;
	//stiffnesses of individual springs
	double Kx[n+1],m[n];
	
	//Initailizing Kx matrix
	for(int i=0;i<n+1;i++){
		Kx[i]=1;
	} 
	//initializing m matrix
	for(int i=0;i<n+1;i++){
		m[i]=1;
	} 


	//outdata.open("new1.csv",ios::app);
	//K and M matrices
	double M[n][n],K[n][n];
	
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

	//Initializing MPI 
	MPI_Init( &argc, &argv);
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD,&my_PE_num);
	MPI_Comm_size(MPI_COMM_WORLD,&n1);
	double u[n] = {0},u_old[n] = {0};
	double v[n] = {0},v_old[n] = {0};
	double a[n] = {0},a_old[n] = {0};
	double F[n] = {0};
	v_old[0] =1;
	int flag[n1][4]={0};
	int main_flag[4]={0};
	
	for(int k=0;k<no_intervals;k++){
		//Predictors for u and v.
		
		for(int i=(my_PE_num)*n/n1;i<=(my_PE_num+1)*n/n1-1;i++){
			u[i]=u_old[i]+Dt*v_old[i]+Dt*Dt*a_old[i]/2;
			v[i]=v_old[i]+Dt*a_old[i]/2;
			u_old[i]=u[i];
		}
		
		if(my_PE_num==0){
			MPI_Send(&u[(my_PE_num+1)*n/n1-1],1,MPI_DOUBLE,my_PE_num+1,70,MPI_COMM_WORLD);
			MPI_Send(&u_old[(my_PE_num+1)*n/n1-1],1,MPI_DOUBLE,my_PE_num+1,100,MPI_COMM_WORLD);
			MPI_Recv(&u[(my_PE_num+1)*n/n1],1,MPI_DOUBLE,my_PE_num+1,10,MPI_COMM_WORLD, &status);
			MPI_Recv(&u_old[(my_PE_num+1)*n/n1],1,MPI_DOUBLE,my_PE_num+1,40,MPI_COMM_WORLD, &status);
		}	
		else if(my_PE_num==n1-1){
			MPI_Send(&u[my_PE_num*n/n1],1,MPI_DOUBLE,my_PE_num-1,10,MPI_COMM_WORLD);
			MPI_Send(&u_old[my_PE_num*n/n1],1,MPI_DOUBLE,my_PE_num-1,40,MPI_COMM_WORLD);	
			MPI_Recv(&u[my_PE_num*n/n1-1],1,MPI_DOUBLE,my_PE_num-1,70,MPI_COMM_WORLD, &status);
			MPI_Recv(&u_old[my_PE_num*n/n1-1],1,MPI_DOUBLE,my_PE_num-1,100,MPI_COMM_WORLD, &status);
		}
		else{	
			MPI_Send(&u[my_PE_num*n/n1],1,MPI_DOUBLE,my_PE_num-1,10,MPI_COMM_WORLD);
			MPI_Send(&u_old[my_PE_num*n/n1],1,MPI_DOUBLE,my_PE_num-1,40,MPI_COMM_WORLD);
			MPI_Send(&u[(my_PE_num+1)*n/n1-1],1,MPI_DOUBLE,my_PE_num+1,70,MPI_COMM_WORLD);
			MPI_Send(&u_old[(my_PE_num+1)*n/n1-1],1,MPI_DOUBLE,my_PE_num+1,100,MPI_COMM_WORLD);
			MPI_Recv(&u[my_PE_num*n/n1-1],1,MPI_DOUBLE,my_PE_num-1,70,MPI_COMM_WORLD, &status);
			MPI_Recv(&u_old[my_PE_num*n/n1-1],1,MPI_DOUBLE,my_PE_num-1,100,MPI_COMM_WORLD, &status);
			MPI_Recv(&u[(my_PE_num+1)*n/n1],1,MPI_DOUBLE,my_PE_num+1,10,MPI_COMM_WORLD, &status);
			MPI_Recv(&u_old[(my_PE_num+1)*n/n1],1,MPI_DOUBLE,my_PE_num+1,40,MPI_COMM_WORLD, &status);
		}

		MPI_Barrier(MPI_COMM_WORLD);
    
		for(int i=(my_PE_num)*n/n1;i<=(my_PE_num+1)*n/n1-1;i++){
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
		if(my_PE_num==0){
			MPI_Send(&a[(my_PE_num+1)*n/n1-1],1,MPI_DOUBLE,my_PE_num+1,90,MPI_COMM_WORLD);
			MPI_Send(&a_old[(my_PE_num+1)*n/n1-1],1,MPI_DOUBLE,my_PE_num+1,120,MPI_COMM_WORLD);	
			MPI_Recv(&a[(my_PE_num+1)*n/n1],1,MPI_DOUBLE,my_PE_num+1,30,MPI_COMM_WORLD, &status);
			MPI_Recv(&a_old[(my_PE_num+1)*n/n1],1,MPI_DOUBLE,my_PE_num+1,60,MPI_COMM_WORLD, &status);
		}	
		else if(my_PE_num==n1-1){
			MPI_Send(&a[my_PE_num*n/n1],1,MPI_DOUBLE,my_PE_num-1,30,MPI_COMM_WORLD);
			MPI_Send(&a_old[my_PE_num*n/n1],1,MPI_DOUBLE,my_PE_num-1,60,MPI_COMM_WORLD);	
			MPI_Recv(&a[my_PE_num*n/n1-1],1,MPI_DOUBLE,my_PE_num-1,90,MPI_COMM_WORLD, &status);
			MPI_Recv(&a_old[my_PE_num*n/n1-1],1,MPI_DOUBLE,my_PE_num-1,120,MPI_COMM_WORLD, &status);
			}
		else{	
			MPI_Send(&a[my_PE_num*n/n1],1,MPI_DOUBLE,my_PE_num-1,30,MPI_COMM_WORLD);
			MPI_Send(&a_old[my_PE_num*n/n1],1,MPI_DOUBLE,my_PE_num-1,60,MPI_COMM_WORLD);	
			MPI_Send(&a[(my_PE_num+1)*n/n1-1],1,MPI_DOUBLE,my_PE_num+1,90,MPI_COMM_WORLD);
			MPI_Send(&a_old[(my_PE_num+1)*n/n1-1],1,MPI_DOUBLE,my_PE_num+1,120,MPI_COMM_WORLD);
			MPI_Recv(&a[my_PE_num*n/n1-1],1,MPI_DOUBLE,my_PE_num-1,90,MPI_COMM_WORLD, &status);
			MPI_Recv(&a_old[my_PE_num*n/n1-1],1,MPI_DOUBLE,my_PE_num-1,120,MPI_COMM_WORLD, &status);
			MPI_Recv(&a[(my_PE_num+1)*n/n1],1,MPI_DOUBLE,my_PE_num+1,30,MPI_COMM_WORLD, &status);
			MPI_Recv(&a_old[(my_PE_num+1)*n/n1],1,MPI_DOUBLE,my_PE_num+1,60,MPI_COMM_WORLD, &status);
		}


		MPI_Barrier(MPI_COMM_WORLD);
		

		//Corrector for u and v
		for(int i=(my_PE_num)*n/n1;i<=(my_PE_num+1)*n/n1-1;i++){
			v[i]=v[i]+a[i]*Dt/2;
			v_old[i]=v[i];
		}


		if(my_PE_num==0){
			MPI_Send(&v[(my_PE_num+1)*n/n1-1],1,MPI_DOUBLE,my_PE_num+1,80,MPI_COMM_WORLD);
			MPI_Send(&v_old[(my_PE_num+1)*n/n1-1],1,MPI_DOUBLE,my_PE_num+1,110,MPI_COMM_WORLD);
			MPI_Recv(&v[(my_PE_num+1)*n/n1],1,MPI_DOUBLE,my_PE_num+1,20,MPI_COMM_WORLD, &status);
			MPI_Recv(&v_old[(my_PE_num+1)*n/n1],1,MPI_DOUBLE,my_PE_num+1,50,MPI_COMM_WORLD, &status);
		}	
		else if(my_PE_num==n1-1){
			MPI_Send(&v[my_PE_num*n/n1],1,MPI_DOUBLE,my_PE_num-1,20,MPI_COMM_WORLD);
			MPI_Send(&v_old[my_PE_num*n/n1],1,MPI_DOUBLE,my_PE_num-1,50,MPI_COMM_WORLD);
			MPI_Recv(&v[my_PE_num*n/n1-1],1,MPI_DOUBLE,my_PE_num-1,80,MPI_COMM_WORLD, &status);
			MPI_Recv(&v_old[my_PE_num*n/n1-1],1,MPI_DOUBLE,my_PE_num-1,110,MPI_COMM_WORLD, &status);
			}
		else{	
			MPI_Send(&v[my_PE_num*n/n1],1,MPI_DOUBLE,my_PE_num-1,20,MPI_COMM_WORLD);
			MPI_Send(&v_old[my_PE_num*n/n1],1,MPI_DOUBLE,my_PE_num-1,50,MPI_COMM_WORLD);
			MPI_Send(&v[(my_PE_num+1)*n/n1-1],1,MPI_DOUBLE,my_PE_num+1,80,MPI_COMM_WORLD);
			MPI_Send(&v_old[(my_PE_num+1)*n/n1-1],1,MPI_DOUBLE,my_PE_num+1,110,MPI_COMM_WORLD);
			MPI_Recv(&v[my_PE_num*n/n1-1],1,MPI_DOUBLE,my_PE_num-1,80,MPI_COMM_WORLD, &status);
			MPI_Recv(&v_old[my_PE_num*n/n1-1],1,MPI_DOUBLE,my_PE_num-1,110,MPI_COMM_WORLD, &status);
			MPI_Recv(&v[(my_PE_num+1)*n/n1],1,MPI_DOUBLE,my_PE_num+1,20,MPI_COMM_WORLD, &status);
			MPI_Recv(&v_old[(my_PE_num+1)*n/n1],1,MPI_DOUBLE,my_PE_num+1,50,MPI_COMM_WORLD, &status);
		}
		MPI_Barrier(MPI_COMM_WORLD);

	}
	cout<<u[my_PE_num]<<endl;
	//cout<<("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC)<<endl;
	MPI_Finalize();
}


