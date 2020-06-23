#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <math.h>
#include "random.h"

using namespace std;

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

    int n = 10000;                          /* Numero totale di numeri pseudocasuali prodotti */
    int M = 100;                            /* Numbero di blocchi in cui suddivido [0, 1] */
    double length = 0.01;
    double E = n/M;
    double x;
    double chiquadro=0;
    int *count = new int[M];                /* contatore per ogni intervallo*/
    ofstream fout("chiquadro.txt");
    
    for (int k=0; k<M; k++){
        for (int i=0; i<M; i++) count[i]=0;     /*lo inizializzo a 0*/
        for (int i=0; i<n; i++) {
            x=rnd.Rannyu();
            for (int j=0; j<M; j++){
                if ((x>(j*length)) && (x<((j+1)*length))){
                    count[j]+=1;
                    j=M;
                }
            }
        }
        
        for (int i=0; i<M; i++){                            /*calcolo di un chiquadro*/
            chiquadro += (((count[i]-E)*(count[i]-E))/E);
        }
        fout<< chiquadro << endl;
        chiquadro = 0;
    }
    fout.close();
    rnd.SaveSeed();
    return 0;
}
                          
