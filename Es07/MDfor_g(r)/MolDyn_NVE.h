/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_gdir;

// averages
double ave_pot, ave_kin, ave_temp, ave_etot;
double media_pot, media_pot2, media_kin, media_kin2;
double media_temp, media_temp2, media_etot, media_etot2;
double bin_size;
int blk_norm;
const int nbins=100;
double walker[nbins], blk_av[nbins];
double glob_av[nbins],glob_av2[nbins];


//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut,fs,eKb, sigma;

// simulation
int nstep, nblocks, iprint, seed,L,restart;
double delta, block;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void Accumulate(void);
void ConfXYZ(int);
void Measure(void);
void Blocco(double);
double Force(int, int);
double Pbc(double);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
