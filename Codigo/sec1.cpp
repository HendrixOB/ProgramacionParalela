#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <bits/stdc++.h>

using namespace std;

//Definimos el producto interno, 
//esto para evitar multiplicar matrices de 1xn por nx1

double productointerno(const vector<double> &V, const vector<double> &W);


//Otra cosa que usa el algoritmo es la multiplicacion de matrices nxn por nx1
vector<double> prodmatriz(const vector<vector<double>> &A,  const vector<double> &V);
 
 
vector<double> suma(const vector<double> &V,const double a, const vector<double> &W);

int main()
{
   	clock_t start, end;
	
   	start = clock();
   	vector<vector<double>> A = {
    	{1, 1, 10, 0, 0, 0, 7, 0, 0, -5},
    	{1, 2, 11, 0, 0, 0, 0, 0, 0, 0},
    	{10, 11, 3, -1, 0, 0, 0, 0, 0, 0},
    	{0, 0, -1, 4, 5, 0, 0, 0, 0, 0},
    	{0, 0, 0, 5, 5, 20, 0, 0, 0, 0},
    	{0, 0, 0, 0, 20, 6, -5, 0, 0, 12},
    	{7, 0, 0, 0, 0, -5, 7, 13, 0, 11},
    	{0, 0, 0, 0, 0, 0, 13, 8, -2, 10},
    	{0, 0, 0, 0, 0, 0, 0, -2, 9, -1},
    	{-5, 0, 0, 0, 0, 12, 11, 10, -1, 10}};
	vector<double> B = { 6, 25, -11, 15, 3, 7, 4, 19, 8, -6 };
   	double error = 1.0e-100;
	int n = A.size();
	int k = 0;
	double alpha = 0.0;
	double beta = 0.0;
	vector<double> r = B;
	vector<double> p = r;
	vector<double> Ap(n, 0.0);
	vector<double> X(n, 0.0);
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

   cout << "X es\n";
   for (const auto i: X){
      cout << i << ' ';
   }
   cout << '\n';
   cout << "Verifiquemos AX nos da b\n";
   
   vector<double> Y = prodmatriz( A, X );
   for (const auto i: Y){
      cout << i << ' ';
   }
   cout << '\n';
   
   end = clock();
  
    	// Calculating total time taken by the program.
    	double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    	cout << "Time taken by program is : " << fixed 
             << time_taken << setprecision(5);
    	cout << " sec " << endl;

 
}

double productointerno(const vector<double> &V, const vector<double> &W){
	int n = V.size(); //definimos el numero de entradas, ya que deben tener las misma
	double prodint = 0.0; //definimos una variable para poder meter el resultado
	for (int i = 0; i < n; i++){
		prodint = prodint + V[i]*W[i];
	}
	return prodint;
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
 	int n = V.size();
 	vector<double> Z(n);
 	for(int i = 0; i < n; i++){
 		Z[i] = V[i] + a*W[i];
 	}
 	return Z;
 }