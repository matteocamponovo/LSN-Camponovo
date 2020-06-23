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
    double dir = 0;
    double x, y, z, theta, phi;
    double a = 1;
    int M = 10000;              /* Total number of throws */
    int N = 100;                 /* Number of blocks */
    int L = (M/N);               /* Number of throws in each block, please use for M a multiple of N */
    int LRW = 100;                /*numero dei passi di ciascun RW*/
    ofstream fout("risultati_continuo.txt");
    double media[LRW][N];
    double media2[LRW][N];
    double r2[LRW][L];
    double* sum_prog = new double[N];
    double* su2_prog = new double[N];
    double* err_prog = new double[N];
    for (int j=0; j<N; j++){
        for (int i=0; i<LRW; i++) media[i][j]=media2[i][j]=0;
    }
    for (int j=0; j<L; j++){
        for (int i=0; i<LRW; i++) r2[i][j]=0;
    }
    
for (int i=0; i<N; i++) {       /*ciclo for per ciclare sui blocchi*/
    
    for (int j=0; j<L; j++) {       /*ciclo for per ciclare sui componenti di un blocco, ossia 100 RW*/
        x=y=z = 0;
        for (int k=0; k<LRW; k++){  /*ciclo per generare le 100 posizioni successive di un RW*/
            dir = rnd.Rannyu();
            theta = acos (1.-2.*rnd.Rannyu());
            phi = rnd.Rannyu()* 2.*M_PI;
            if (dir > 0.5){
                x = x + a*sin(theta)*cos(phi);
                y = x + a*sin(theta)*sin(phi);
                z = z + a*cos(theta);
            }
            else {
                x = x - a*sin(theta)*cos(phi);
                y = x - a*sin(theta)*sin(phi);
                z = z - a*cos(theta);
            }
                r2[k][j] = (x*x + y*y + z*z);
        }
    }
    
    for (int fix=0; fix<LRW; fix++){                /*LRW*/
        sum =0;
        for (int j=0; j<L; j++)sum+=r2[fix][j];
        media[fix][i]=sqrt(sum/L);
        media2[fix][i]=media[fix][i]*media[fix][i];
    }
}
    for (int k=0; k<LRW; k++){
        for (int i=0; i<N; i++){
            for (int j=0; j<i+1; j++){
                sum_prog[i] += media[k][j];
                su2_prog[i] += media2[k][j];
            }
            sum_prog[i]/=(i+1);                              /* Cumulative average */
            su2_prog[i]/=(i+1);                              /* Cumulative square average */
            err_prog[i] = error(su2_prog,sum_prog,i);        /* Statistical uncertainty */
        }
        fout<< k+1 << "   "<< sum_prog[N-1] << "    "<< err_prog[N-1]<< endl;
        for (int i=0; i<N; i++)sum_prog[i]=su2_prog[i]=err_prog[i]=0;
    }
fout.close();
rnd.SaveSeed();
return 0;
}
                          
