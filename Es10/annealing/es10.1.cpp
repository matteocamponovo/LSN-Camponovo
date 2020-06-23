#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include<algorithm>
#include "es10.1.h" 

using namespace std;

int main()
{
  Input(); //Inizialization
  ofstream Best;
	Best.open("Best.dat");
	beta=0.1;
	for(int i=0;i<N_t;i++){
  	Annealing(beta);
		cout<<"Temperatura :"<<"\t"<<1./beta<<endl;
		cout<<"Lunghezza percorso :"<<"\t"<<male.Cost_fn_L1(Mappa)<<endl;
		cout<<"------------------------------------"<<endl;
		Best<<setprecision(9)<<i<<"\t"<<male.Cost_fn_L1(Mappa)<<endl;
		beta=beta*1.1;
	}
	Best.close();
	male.Save_conf(Mappa);
  delete rnd;
  return 0;
}


void Input(void){
  
  ifstream ReadInput;

//Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  rnd =new Random;
  
  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3]; 
  rnd->SetRandom(seed,p1,p2);
   
  input.close();   
  //Read input informations
  ReadInput.open("input.dat");

  ReadInput >> N_c;

  ReadInput >> nstep;

  ReadInput >> rect;

  cout << "Salesman problem" << endl;
  cout << "Numero di città da visitare= " <<N_c << endl;
  cout << "Numero di step per temperatura = " << nstep << endl;
  cout << "Distribuzione città = "<<rect<<endl;

  ReadInput.close();

  Mappa.Set_par(N_c,rect,rnd);

  male=Chrom(N_c,rnd);

  return;
}

void Annealing(double beta){
	double p,alpha;
	accepted=0.;	
	attempted=0.;

	for(int i=0;i<nstep;i++){
		new_male=male;

		new_male.Pair_permutation();
		new_male.Shift();
		new_male.m_permutation();
		new_male.Inversion();

        p=exp(beta*(male.Cost_fn_L1(Mappa)-new_male.Cost_fn_L1(Mappa)));
		alpha=min(1.,p);
		if(alpha >= rnd->Rannyu()) {	
            male=new_male;
            accepted = accepted + 1.0;
        }
        attempted= attempted + 1.0;
    }
}






