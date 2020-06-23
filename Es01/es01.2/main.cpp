#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <math.h>
#include "random.h"
using namespace std;

double exponential(double y, double lambda){
    return -(1/lambda)*log(1-y);
}
double CauchyLorentz(double y, double gamma, double mu){
    return gamma*tan(M_PI*(y-0.5))+mu;
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

    int M = 10000;
    ofstream Nun("N_un.txt");
    ofstream Nexp("N_exp.txt");
    ofstream NCL("N_CL.txt");
    double *S1 = new double[M];
    double *S2 = new double[M];
    double *S10 = new double[M];
    double *S100 = new double[M];
    double *S1exp = new double[M];
    double *S2exp = new double[M];
    double *S10exp = new double[M];
    double *S100exp = new double[M];
    double *S1CL = new double[M];
    double *S2CL = new double[M];
    double *S10CL = new double[M];
    double *S100CL = new double[M];
    double appo, appoCL, appoexp, a, b;
    for (int i=0; i<M; i++)S1[i]=S2[i]=S10[i]=S100[i]=S1exp[i]=S2exp[i]=S10exp[i]=S100exp[i]=S1CL[i]=S2CL[i]=S10CL[i]=S100CL[i]=0;
    
   /* riempio S1 per tutte e tre le distrubuzioni */
    for (int i=0; i<M; i++){
        S1[i]=rnd.Rannyu();
        S1exp[i]=exponential(S1[i],1);
        S1CL[i]=CauchyLorentz(S1[i],1,0);
    }
    
    /* riempio S2*/
    for (int k=0; k<M; k++) {
        a=rnd.Rannyu();
        b=rnd.Rannyu();
        S2[k]=(a+b)/2;
        S2exp[k]=(exponential(a,1)+exponential(b,1))/2;
        S2CL[k]=(CauchyLorentz(a,1,0)+CauchyLorentz(b,1,0))/2;
    }
    
    /* riempio S10*/
    for (int k=0; k<M; k++) {
        for (int j=0; j<10; j++){
            appo=rnd.Rannyu();
            appoexp=exponential(appo,1);
            appoCL=CauchyLorentz(appo,1,0);
            S10[k]+=appo;
            S10exp[k]+=appoexp;
            S10CL[k]+=appoCL;
        }
        S10[k]=S10[k]/10;
        S10exp[k]=S10exp[k]/10;
        S10CL[k]=S10CL[k]/10;
    }
    
    /* riempio S100 sommando a 100 a 100 su S1 */
     for (int k=0; k<M; k++) {
         for (int j=0; j<100; j++){
               appo=rnd.Rannyu();
               appoexp=exponential(appo,1);
               appoCL=CauchyLorentz(appo,1,0);
               S100[k]+=appo;
               S100exp[k]+=appoexp;
               S100CL[k]+=appoCL;
           }
           S100[k]=S100[k]/100;
           S100exp[k]=S100exp[k]/100;
           S100CL[k]=S100CL[k]/100;
       }
    
    /*caricamento file*/
    for (int i=0; i<M; i++) Nun<< S1[i] <<" "<<S2[i]<<" "<<S10[i]<<" "<<S100[i]<<endl;
    for (int i=0; i<M; i++) Nexp<< S1exp[i] << " "<<S2exp[i]<< " "<<S10exp[i]<<" "<<S100exp[i]<<endl;
    for (int i=0; i<M; i++) NCL<< S1CL[i] << " "<<S2CL[i]<< " "<<S10CL[i]<<" "<<S100CL[i]<<endl;

    Nun.close();
    Nexp.close();
    NCL.close();
    rnd.SaveSeed();
    return 0;
}
                          
