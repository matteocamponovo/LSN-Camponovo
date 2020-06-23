#include "es5.h"

int main (int argc, char *argv[]){
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

//----------------------------  STATI 1S e 2P IDROGENO  ------------------------------------------------------------

	Input();
    double* ave = new double[N];
    double* ave2 = new double[N];
    double* sum_prog = new double[N];
    double* su2_prog = new double[N];
    double* err_prog = new double[N];
    double* r = new double[M];
    
    for (int i=0; i<M; i++) r[i]=0;
    for (int i=0; i<N; i++) ave[i]=ave2[i]=sum_prog[i]=su2_prog[i]=err_prog[i]=0;
	conta=0;
	eq=0;
    delta=0;
	Equilibrazione();
	cout<<"ho modificato delta "<<eq<<" volte"<< " delta vale  " <<delta<<endl;
	L=M/N;
    ofstream out1;
    if (nome_stato=="1S"){
        out1.open("grafico3D1S.dat", ios::app);}
    if (nome_stato=="2P"){
        out1.open("grafico3D2P.dat", ios::app);}
    
	for (int i=0; i<M; i++) {
		Passo();
        out1<<x<<"    "<<y<<"    "<<z<<endl;
        r[i]=sqrt(x*x+y*y+z*z);
    }
    out1.close();
    
    for (int i=0; i<N; i++) {
        sum = 0;
        for (int j=0; j<L; j++) {
            int k = j+i*L;
            sum += r[k];
        }
        ave[i] = sum/L;
        ave2[i] = (ave[i])*(ave[i]);
    }

    for (int i=0; i<N; i++){
        for (int j=0; j<i+1; j++){
            sum_prog[i] += ave[j];
            su2_prog[i] += ave2[j];
        }
        sum_prog[i]/=(i+1);                              /* Cumulative average */
        su2_prog[i]/=(i+1);                              /* Cumulative square average */
        err_prog[i] = error(su2_prog,sum_prog,i);        /* Statistical uncertainty */
    }

    if (nome_stato=="1S"){
        ofstream fout1S("rmedio_1S.txt");
        for (int i=0; i<N; i++) fout1S<< i<< "   "<< sum_prog[i] << "    "<< err_prog[i]<< endl;
        fout1S.close();
    }
    if (nome_stato=="2P"){
        ofstream fout2P("rmedio_2P.txt");
        for (int i=0; i<N; i++) fout2P<< i << "   "<< sum_prog[i] << "    "<< err_prog[i]<< endl;
        fout2P.close();
    }
    /*cout<<"frazione di mosse accettate= "<<(double)conta/(double)M<<endl;*/

   rnd.SaveSeed();
   return 0;
}



double psi1S (double r) {
	return (1.0/M_PI)*exp(-2.0*r);
}
double psi2P (double x, double y, double z) {
	double radius=sqrt(x*x+y*y+z*z);
	return ((1.0/(32.0*M_PI))*radius*radius*exp(-radius)*(z*z/(radius*radius)));
}

void Input (){
	ifstream in("input.txt");
	in>>delta;
	in>>x;
	in>>y;
	in>>z;
	in>>M;
	in>>N;
	in>>nome_stato;
	in>>distribuzione;
}

void Passo () {
	if (distribuzione==1) {
		xtest=rnd.Rannyu(x-delta, x+delta);
		ytest=rnd.Rannyu(y-delta, y+delta);
		ztest=rnd.Rannyu(z-delta, z+delta);
	}
	if (distribuzione==2) {
		xtest=rnd.Gauss(x, delta);
		ytest=rnd.Gauss(y, delta);
		ztest=rnd.Gauss(z, delta);
	}
	if (nome_stato=="1S") {
		double q=sqrt(x*x+y*y+z*z);
		double s=sqrt(xtest*xtest+ytest*ytest+ztest*ztest);
		alpha=min(1.,psi1S(s)/psi1S(q));
	}
	if (nome_stato=="2P") {alpha=min(1.,psi2P(xtest, ytest, ztest)/psi2P(x, y, z));}	
	double p=rnd.Rannyu(0,1);
	if (p<=alpha) {
		x=xtest;
		y=ytest;
		z=ztest;
		conta++;
	}
}

void Equilibrazione(){
	double lanci=1000;
	for (int j=0; j<lanci; j++) {
			Passo();
		}
	double rapporto=(double)conta/(double)lanci;
	if (rapporto>=0.48 && rapporto<=0.52) {return;} 
	do{	
		conta=0;
		delta+=0.1;
		for (int j=0; j<lanci; j++) {
			Passo();
		}
		rapporto=(double)conta/(double)lanci;
        cout << rapporto << endl;
		eq++;
	}
	while (rapporto<=0.48 || rapporto>=0.52); 
}

double error(double * AV2, double * AV, int n){           /* Function for statistical uncertainty estimation */
    if (n==0)
        return 0;
    else
        return sqrt((AV2[n] - AV[n]*AV[n])/n);
}
