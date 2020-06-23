#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <math.h>
#include <algorithm>
#include "es9.h"

using namespace std;
    

int main() {
    
   Input();

   if (!Check()){
       cerr << "C'è un errore nella costruzione della popolazione!" << endl;
       return -1;
   }
    
    ofstream Best, Half;
    Best.open("Best.dat");
    Half.open("Half.dat");
    
    for(int i=0;i<N_gen;i++){

           Generation();

           if (!Check()){
                cerr << "C'è un errore nella costruzione della popolazione!" << endl;
                return -1;
           }
           cout<<"Generazione :"<<"\t"<<i<<endl;
           cout<<"Lunghezza percorso :"<<"\t"<<pop[0].Cost_fn_L1(Mappa)<<endl;
           cout<<"------------------------------------"<<endl;
           Best<<setprecision(9)<<i+1<<"\t"<<pop[0].Cost_fn_L1(Mappa)<<endl;
           Half<<setprecision(9)<<i+1<<"\t"<<Half_mean()<<endl;
    }
    
    pop[0].Save_conf(Mappa);
    Half.close();
    Best.close();
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
    
//Leggo informazioni in input
  ReadInput.open("input.dat");

  ReadInput >> N_p; //numero di cromosomi/individui/traiettorie

  ReadInput >> N_c; //numero di città

  ReadInput >> N_gen; //numero di volte che esegue il codice genetico

  ReadInput >> rect; //quadrato o circonferenza


  cout << "PROBLEMA DEL COMMESSO VIAGGIATORE" << endl;
  cout << "Numero di città da visitare = " <<N_c << endl;
  cout << "Numero di cromosomi nella popolazione = " << N_p << endl;
  cout << "Numero di generazioni = " << N_gen << endl;
  cout << "Distribuzione città = "<<rect<<endl;

  ReadInput.close();
  Mappa.Set_par(N_c,rect,rnd); //genero la prima mappa di città, ossia posiziono le
                               //coordinate sulla circonferenza o all'interno del quadrato.
  //pop.resize(N_p);
  for(int i=0; i<N_p; i++) pop.push_back(Chrom(N_c,rnd)); //carico tutti i cromosomi nella popolazione
  sort(pop.begin(), pop.end(), sort_fitness);            //li riordino in base al fitness
}

bool Check(){
    for (int i=0; i<pop.size(); i++){
        if(pop[i].Get_Element(0) != 1){
            
            cout <<pop[i].Get_Element(0) <<endl;
            return false;
        }
        for (int j=0; j<pop[i].Get_Size()-1; j++){
            for (int z=j+1; z<pop[i].Get_Size(); z++){
                if (pop[i].Get_Element(j)==pop[i].Get_Element(z)) return false;
            }
        }
    }
    return true;
}

void Generation(void){
  int m,f,cross;
  vector<Chrom> pop_new;
  for(int i=0; i<N_p/2; i++) {           //selection (for su N_p/2 perché ad ogi giro genero due inidividui nuovi)
        m= int(pow(rnd->Rannyu(0.,1.),r)*N_p);    //selezione madre del crossover
        f= int(pow(rnd->Rannyu(0.,1.),r)*N_p);    //padre madre del crossover
        while(m==f) f= int(pow(rnd->Rannyu(0.,1.),r)*N_p);
        Chrom parent_m=pop[m];
        Chrom parent_f=pop[f];
        
        if(P_c>rnd->Rannyu()){                  //crossover
            cross=(int)rnd->Rannyu(0.,N_c-1);   //punto di taglio dei cromosomi
            pop_new.push_back(Chrom(parent_m,parent_f,cross,rnd));  //generazione primo figlio
            pop_new.push_back(Chrom(parent_f,parent_m,cross,rnd));  //generazione secondo figlio
        }else{
            pop_new.push_back(parent_m);    //se non avviene crossover riaggiungo padre e madre
            pop_new.push_back(parent_f);
        }
    }

    for(int i=0;i<N_p;i++){
        pop_new[i].Pair_permutation();
        pop_new[i].Shift();
        pop_new[i].m_permutation();
        pop_new[i].Inversion();
    }
 
    sort(pop_new.begin(), pop_new.end(), sort_fitness);   //riordino la nuova pop in base al fitness
    
    if(P_e>rnd->Rannyu()){//elite
        pop_new.back()=pop.front();
        sort(pop_new.begin(), pop_new.end(), sort_fitness);
    }

    copy(pop_new.begin(),pop_new.end(),pop.begin());
}


double Half_mean(void){             //funzione per restituire il valore medio di L1 nella metà migliore di pop
                                 //cioè dove ci sono i percorsi minimizzano la funzione costo
    double mean_val=0;
    for(int i=0;i<N_p/2;i++) mean_val+=pop[i].Cost_fn_L1(Mappa);
    return mean_val/double(N_p/2);
}

bool sort_fitness( Chrom &former, Chrom &latter){              //mi permette di ordinare con questo criterio la pop
    return (former.Cost_fn_L1(Mappa) < latter.Cost_fn_L1(Mappa));
}  
