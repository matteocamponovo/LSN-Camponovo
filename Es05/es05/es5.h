#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;
Random rnd;
int conta;
int eq;
int distribuzione;
string nome_stato;
double x,y,z, delta;
int M, N, L;
double xtest, ytest, ztest, alpha;
double sum;

void Input ();
void Passo ();
double psi1S (double);
double psi2P (double, double, double);
void Equilibrazione ();
double error(double *, double *, int);
