/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"
using namespace std;

int main() {
  Input();            //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
    Reset(iblk);      //Reset block averages
    for(int istep=1; istep <= nstep; ++istep){
      Move(metro);
      Measure();
      Accumulate();   //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal();        //Write final configuration
  return 0;
}

// ****************DICHIARAZIONE DI TUTTE LE FUNZIONI******************

void Input(void) {
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

  ReadInput >> restart;
  ReadInput.close();
//initial or final configuration
  if(restart == 0){
      for (int i=0; i<nspin; ++i){
          if(rnd.Rannyu() >= 0.5) s[i] = 1;
          else s[i] = -1;
      }
  } else {
    ReadInput.open("config.final");
    for (int i=0; i<nspin; ++i){
      ReadInput >> s[i];
    }
    ReadInput.close();
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro){
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i){
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    attempted++;

    if(metro==1) //Metropolis
    {
        energy_old = Boltzmann( s[o], o );
        energy_new = Boltzmann( (-1) * s[o], o );
        p = min( 1. , exp(-beta * energy_new) / exp(-beta * energy_old)  );
        //cout << "prob Metropolis!   " << p <<endl;
        if(rnd.Rannyu() <= p){
            s[o] *= (-1.);
            accepted++;
        }
    }
    else {        //Gibbs sampling
        p = 1./( 1. + exp( - 2. * beta * (J * ( s[Pbc(o-1)] + s[Pbc(o+1)] )+h ) ));
        if(rnd.Rannyu() <= p){
            s[o] = 1.;
        } else {
            s[o] = -1.;
        }
    }
  }
}

double Boltzmann(int sm, int ip) {
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure() {
  int bin;
  double u = 0.0, u2 = 0.0, m = 0.0;
//cycle over spins
  for (int i=0; i<nspin; ++i) {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     u2 += pow(-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]) , 2 );
     m += s[i] ;
  }
  walker[iu] = u;
  walker[ic] = -( beta * beta ) * (u*u/(double)(nspin*nspin) - u2/(double)nspin );
  walker[im] = m;
  walker[ix] = beta * m * m;
}


void Reset(int iblk) { //Reset block averages
   if(iblk == 1) {
       for(int i=0; i<n_props; ++i){
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }
   for(int i=0; i<n_props; ++i) blk_av[i] = 0;

   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) { //Update block averages
   for(int i=0; i<n_props; ++i) blk_av[i] = blk_av[i] + walker[i];

   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
   ofstream Ene, Heat, Mag, Chi, final_ene, final_heat, final_chi, final_M, Gfinal_ene, Gfinal_heat, Gfinal_chi, Gfinal_M;
   const int wd=15;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    Heat.open("output.heat.0",ios::app);
    stima_c = blk_av[ic]/blk_norm; //Specific Heat
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();

    Chi.open("output.chi.0",ios::app);
    stima_x = blk_av[ix]/blk_norm/(double)nspin; //Magnetic susc.
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();
    
    Mag.open("output.mag.0",ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();

    cout << "----------------------------" << endl << endl;
    
    if (iblk == nblk && metro==1){
        final_ene.open("Mfinal.ene.0", ios::app);
        final_ene << setw(wd) << temp <<  setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
        final_heat.open("Mfinal.heat.0", ios::app);
        final_heat << setw(wd) << temp <<  setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
        final_chi.open("Mfinal.chi.0", ios::app);
        final_chi << setw(wd) << temp <<  setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
        final_M.open("Mfinal.mag.0", ios::app);
        final_M << setw(wd) << temp <<  setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    }else if(iblk == nblk && metro!=0){
        Gfinal_ene.open("Gfinal.ene.0", ios::app);
        Gfinal_ene << setw(wd) << temp <<  setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
        Gfinal_heat.open("Gfinal.heat.0", ios::app);
        Gfinal_heat << setw(wd) << temp <<  setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
        Gfinal_chi.open("Gfinal.chi.0", ios::app);
        Gfinal_chi << setw(wd) << temp <<  setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
        Gfinal_M.open("Gfinal.mag.0", ios::app);
        Gfinal_M << setw(wd) << temp <<  setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    }
    
        
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
