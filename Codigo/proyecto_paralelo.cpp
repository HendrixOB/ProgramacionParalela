#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "mpi.h"
#include <stdio.h>

using namespace std;
double productointerno(const vector<double> &V, const vector<double> &W);
vector<double> prodmatriz(const vector<vector<double>> &A,  const vector<double> &V);
vector<double> suma(const vector<double> &V,const double a, const vector<double> &W);


int main(int argc, char *argv[])
{
	MPI_Status status;
    	int rank;                                 // Almacenamos el rank de cada proceso
	int size;                                 // Almacenamos el número de procesos 
    	double tic = 0.0;                         // Almacenamos cuando inicia la ejecución
	double toc = 0.0;                         // Almacenamos cuando termina la ejecución
	double tictoc = 0.0;                      // Tiempo total de ejecucion
    	MPI_Init(&argc, &argv);                   // Inicio del entorno MPI
    	MPI_Comm_size(MPI_COMM_WORLD, &size);     // Almacen del número del proceso
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);     // Almacen del rank el proceso
    	tic = MPI_Wtime();                        // Medimos el tiempo inicial de la ejecución
    	
    	//Definimos la matriz a utilizar
    	vector<vector<double>> A = { {1, 1, 10, 0, 0, 0, 7, 0, 0, -5},
    	{1, 2, 11, 0, 0, 0, 0, 0, 0, 0},
    	{10, 11, 3, -1, 0, 0, 0, 0, 0, 0},
    	{0, 0, -1, 4, 5, 0, 0, 0, 0, 0},
    	{0, 0, 0, 5, 5, 20, 0, 0, 0, 0},
    	{0, 0, 0, 0, 20, 6, -5, 0, 0, 12},
    	{7, 0, 0, 0, 0, -5, 7, 13, 0, 11},
    	{0, 0, 0, 0, 0, 0, 13, 8, -2, 10},
    	{0, 0, 0, 0, 0, 0, 0, -2, 9, -1},
    	{-5, 0, 0, 0, 0, 12, 11, 10, -1, 10}};
    
    	//Definimos el vector resultado 
   	vector<double> B = { 6, 25, -11, 15, 3, 7, 4, 19, 8, -6 };
   	double error = 1.0e-100;                  //Definimos el error
	int n = A.size();                         //Definimos la dimension de la matriz
	int k = 0;                                //K nos dara el numero de iteraciones
	double alpha = 0.0;                       //variable que se usara en el metodo
	double beta = 0.0;                        //variable que se usara en el metodo
	vector<double> r = B;                     //Hacemos el reciduo igual a b
	vector<double> p = r;                     //Hacemos p=r
	vector<double> Ap(n, 0.0);                //Definimos un nuevo vector
	vector<double> X(n, 0.0);                 //Definimos el vector que queremos obtener


    if (rank == 0)
    {
	while(k < n){
			vector<double> r0 = r; 
			Ap = prodmatriz(A,p);
			alpha = productointerno(r,r)/productointerno(p,Ap);
			X = suma(X,alpha,p);
			r = suma(r, -alpha,Ap);
			double norma = sqrt(productointerno(r,r));
			if (norma < error){
				break;
			}
			beta = productointerno(r,r)/productointerno(r0,r0);
			p = suma(r,beta,p);
			k++;
		}
    }

    if (rank == 0)
    {
	cout << "X es\n";
   	for (const auto i: X){
      		cout << i << ' ';
   	}
   	cout << '\n';
	toc = MPI_Wtime();
	tictoc = toc - tic;
	if(rank == 0){
		printf("Tiempo total: %f\n",tictoc);		
	}
    }
    MPI_Finalize();
    //return 0;
}

double productointerno(const vector<double> &V, const vector<double> &W){
	int size, rank;
    	MPI_Status status;
    	double prodint = 0.0;
    	double prod = 0.0;
    	MPI_Comm_size(MPI_COMM_WORLD, &size);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    	int n = V.size();
    	for(int i = 0; i < n; i++){
    		prodint += V[i]*W[i];
    	}
    	if(rank == 0){
    		prod = prodint;
    		for(int j = 1; j < size; j++){
    			MPI_Recv(&prodint, 1, MPI_DOUBLE,j,0,MPI_COMM_WORLD, &status);
    			prod += prodint;
    		}
    	}
    	else{
    		MPI_Send(&prodint,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    	}
    	return prod;
}


vector<double> prodmatriz(const vector<vector<double>> &A,  const vector<double> &V){
 	int n = A.size(); 
 	vector<double> G( n );
 	for (int i = 0; i < n; i++){
 		G[i] = productointerno(A[i], V);
 	}
 	return G;
 }
 
 

vector<double> suma(const vector<double> &V,const double a, const vector<double> &W){
 	int size, rank;
    	MPI_Status status;
    	int n = V.size();
    	vector<double> Z(n);
    	vector<double> Y(n);
    	MPI_Comm_size(MPI_COMM_WORLD, &size);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    	for(int i = 0; i < n; i++){
    		Z[i] = V[i]+a*W[i];
    	}
    	if(rank == 0){
    		Y = Z;
    		for(int j = 1; j < size; j++){
    			MPI_Recv(&Z[j], n, MPI_DOUBLE,j,0,MPI_COMM_WORLD, &status);
    			Y[j] = Z[j];
    		}
    	}
    	else{
    		MPI_Send(&Y,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    	}
    	return Z;
}
