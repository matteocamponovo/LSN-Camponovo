#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <math.h>
#include "random.h"

using namespace std;

 double error(double * AV2, double * AV, int n){           /* Function for statistical uncertainty estimation */
     if (n==0)
         return 0;
     else
         return sqrt((AV2[n] - AV[n]*AV[n])/n);
 }

int main() {
    Random rnd;
        int seed[4];
        int p1, p2;
        ifstream Primes("Primes");
        if (Primes.is_open()){
            Primes >> p1 >> p2 ;
        } else cerr << "PROBLEM: Unable to open Primes" << endl;
        Primes.close();
        ifstream input("seed.in");
        string property;
        if (input.is_open()){
            while ( !input.eof() ){
            input >> property;
                if( property == "RANDOMSEED" ){
                    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                    rnd.SetRandom(seed,p1,p2);
                }
            }
            input.close();
        } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    /* sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss*/

    double sum = 0;
    int M = 100000;              /* Total number of throws */
    int N = 100;                 /* Number of blocks */
    int L = (M/N);               /* Number of throws in each block, please use for M a multiple of N */
    double *r = new double[M];   /* U[0,1) uniform distribution */
    for (int i=0; i<M; i++) {
        r[i]=rnd.Rannyu();
    }

    double* ave = new double[N];
    double* ave2 = new double[N];
    double* sum_prog = new double[N];
    double* su2_prog = new double[N];
    double* err_prog = new double[N];
    for (int i=0; i<N; i++){
        ave[i]=ave2[i]=sum_prog[i]=su2_prog[i]=err_prog[i]=0;
    }

    for (int i=0; i<N; i++) {
        sum = 0;
        for (int j=0; j<L; j++) {
            int k = j+i*L;
            sum += (r[k]-0.5)*(r[k]-0.5);
        }
        ave[i] = sum/L;
        ave2[i] = (ave[i])*(ave[i]);
    }

    for (int i=0; i<N; i++){
        for (int j=0; j<i+1; j++){
            sum_prog[i] += ave[j];
            su2_prog[i] += ave2[j];
        }
        sum_prog[i]/=(i+1);                              /* Cumulative average */
        su2_prog[i]/=(i+1);                              /* Cumulative square average */
        err_prog[i] = error(su2_prog,sum_prog,i);        /* Statistical uncertainty */
    }
    ofstream fout("sigmamedio.txt");
    for (int i=0; i<N; i++) fout<< i*L << "   "<< sum_prog[i] << "    "<< err_prog[i]<< endl;
    fout.close();
    rnd.SaveSeed();
    return 0;
}
                          
