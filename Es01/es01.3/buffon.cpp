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
    /* sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss*/
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

    double x =0;
    double y =0;
    double X =0;
    double Y =0;
    double phi =0;
    double d =1;
    double l = 0.8;
    int M = 100000;              /* Total number of throws */
    int N = 100;                 /* Number of blocks */
    double L = (M/N);            /* Number of throws in each block, please use for M a multiple of N */
    double Nhit =0;
    
    double* ave = new double[N];
    double* ave2 = new double[N];
    double* sum_prog = new double[N];
    double* su2_prog = new double[N];
    double* err_prog = new double[N];
    for (int i=0; i<N; i++){
        ave[i]=ave2[i]=sum_prog[i]=su2_prog[i]=err_prog[i]=0;
    }

for (int i=0; i<N; i++){
    Nhit = 0;
    for (int j=0; j<L; j++) {

        x = rnd.Rannyu()*(d/2);         /*genero numero casuale tra o e d/2*/
    
        X = -1 + 2*rnd.Rannyu();          /*genero angolo casuale tra o e pi*/
        Y = rnd.Rannyu();
        if ((X*X+Y*Y)<=1) phi = acos(X/sqrt(X*X+Y*Y));
        while ((X*X+Y*Y)>1) {
            X = -1+2*rnd.Rannyu();
            Y = rnd.Rannyu();
            phi = acos(X/sqrt(X*X+Y*Y));
        }
        y = l/2 * sin(phi);            /*calcolo distanza centro ago dalla riga*/
    
        if (x<y) Nhit+=1;
    }
    ave[i] = 2*l/(d*(Nhit/L));
    ave2[i] = (ave[i])*(ave[i]);
    cout << ave[i] << "   " << ave2[i] << endl;
}
    

    for (int i=0; i<N; i++){
        for (int j=0; j<i+1; j++){
            sum_prog[i] += ave[j];
            su2_prog[i] += ave2[j];
        }
        sum_prog[i]/=(i+1);                              /* Cumulative average */
        su2_prog[i]/=(i+1);                              /* Cumulative square average */
        /*cout << sum_prog[i]*sum_prog[i] << " " << su2_prog[i] << endl;*/

        err_prog[i] = error(su2_prog,sum_prog,i);        /* Statistical uncertainty */
    }
    ofstream fout("pigreco.txt");
    for (int i=0; i<N; i++) fout<< i*L << "   "<< sum_prog[i] << "    "<< err_prog[i]<< endl;
    fout.close();
    rnd.SaveSeed();
    return 0;
}
                          
