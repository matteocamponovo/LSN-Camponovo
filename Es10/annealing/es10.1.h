/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "random.h"
#include "Chrom.h"
#include "Maps.h"
#include <vector>

int seed[4];
Random *rnd;

//parameters, observables
int N_p, rect, N_c, nstep;
const int N_t=100;

double P_m=0.1;
double accepted,attempted, beta;

Chrom new_male;
Chrom male;

Maps Mappa;

//functions
void Input(void); 
void Annealing(double);

