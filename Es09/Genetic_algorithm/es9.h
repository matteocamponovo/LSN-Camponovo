#include "random.h"
#include "Chrom.h"
#include "Maps.h"
#include <vector>

int seed[4];
Random *rnd;

int N_p, N_c, rect, N_gen;

double r=2.; //esponente per la roulette truccata
double P_m=0.1;
double P_c=0.95;
double P_e=0.05;

std::vector<Chrom> pop; //pop è un vettore di Chrom, quindi ogni elemento di pop
                      //è un Chrom, ossia una sequenza di 32 città.

Maps Mappa;

//metodi
void Input(void);
void Generation(void);
bool sort_fitness( Chrom &, Chrom &);
double Half_mean(void);
bool Check();
