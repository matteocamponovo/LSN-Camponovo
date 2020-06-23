#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <ostream>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include "es10.2.h"

using namespace std;
    

int main(int argc, char *argv[]){
    
    MPI_Init(&argc,&argv);

    int my_rank,size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //MPI_Status stat;
    Init_par(my_rank);
    Input(my_rank);
    if (!Check()){
	for (int j=0; j<N_p; j++){
	  for (int i=0; i<N_c; i++) cout <<"pop"<<j<<"->"<<pop[j].Get_Element(i)<<"  ";}
       cerr << "C'è un errore nella costruzione della popolazione!" << endl;
       return -1;
    }
   
    ofstream Best;
    Best.open("Best_"+to_string(my_rank)+".dat");
    //for (int i = 0; i<N_c; i++) cout <<my_rank<<"->"<< pop[0].Get_Element(i)<<"   ";
   
    for(int i=0;i<N_gen;i++){

         Generation();

         if (!Check()){
		for (int i =0; i<N_c; i++) cout <<my_rank<<"->"<<pop[0].Get_Element(i)<<"  ";
                cerr << "C'è un errore nella popolazione!" << endl;
                return -1;
         }
        
         if ((i%N_migr == 0) && (size!=1)) migration(my_rank, size);
   
         Best<<setprecision(9)<<i+1<<"\t"<<pop[0].Cost_fn_L1(Mappa)<<endl;
    }
    
    pop[0].Save_conf(Mappa, my_rank, rect);
    Best.close();
    //for (int i =0; i<N_c; i++) cout<<my_rank<<"->"<< pop[0].Get_Element(i)<<"  ";

    double val = pop[0].Cost_fn_L1(Mappa);
    //local_loss.rank = my_rank;
    cout<<"Rank :"<<"\t"<<my_rank<<"  Lunghezza percorso :"<<"\t"<<val<<endl;
    cout<<"------------------------------------"<<endl;
    
    MPI_Reduce(&val,&Mintot,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
    if(0 == my_rank){
        cout<<"Ho trovato il percorso migliore in assoluto"<<endl;
        cout<<"La sua lunghezza è "<<Mintot<<endl;
       // pop[0].Save_conf_finale(Mappa);
    } 

    delete rnd;
    MPI_Finalize();
    return 0;
}

/*---------------------------------------------------------------------------------------------------------*/


void migration(int my_rank, int size){
    if (my_rank == 0){
        sort1 = int(rnd->Rannyu(0,size));
        sort2 = int(rnd->Rannyu(0,size));
        while (sort1==sort2) sort2 = int(rnd->Rannyu(0,size));
    }

    MPI_Bcast(&sort1,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&sort2,1,MPI_INT,0,MPI_COMM_WORLD);
    int *appo1 = new int [N_c];
    int *appo2 = new int [N_c];
   
    if (my_rank == sort1){
        for (int i=0; i<N_c; i++) appo1[i]=(pop[0].Get_Element(i));
        MPI_Send(appo1,32,MPI_INTEGER,sort2,1,MPI_COMM_WORLD);

    }if (my_rank == sort2){
	 MPI_Recv(appo1,32,MPI_INTEGER,sort1,1,MPI_COMM_WORLD, &stat);
    	for (int i=0; i<N_c; i++) appo2[i] = (pop[0].Get_Element(i)); 
	MPI_Send(appo2,32,MPI_INTEGER,sort1,2,MPI_COMM_WORLD);

    } if (my_rank == sort1){
	MPI_Recv(appo2,32,MPI_INTEGER,sort2,2,MPI_COMM_WORLD,&stat);
    }
   
    if(my_rank==sort1)	{
	 for (int i =0; i<N_c; i++){
            if (i==0) { pop[0].Load_Element0(appo2[i]);
            }else if (i!=0){ pop[0].Load_Element(appo2[i]);
            }
        }
    }

    if (my_rank == sort2){
        for (int i =0; i<N_c; i++){
            if (i==0){ pop[0].Load_Element0(appo1[i]);
            }else if (i!=0){ pop[0].Load_Element(appo1[i]);
            }
        }
    }
}

void Init_par(int my_rank){
    int p1,p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    //cout << p1 << "   " << p2 << endl;
    rnd =new Random;

    if(my_rank==0){
         ifstream input("seed.in");
         input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
         input.close();

         ifstream ReadInput;
         ReadInput.open("input.dat");
        
         ReadInput >> N_p; //numero di cromosomi/individui/traiettorie
         ReadInput >> N_c; //numero città
         ReadInput >> N_gen; //numero di generazioni
         ReadInput >> rect; //quadrato o circonferenza
        
         ReadInput.close();
    }
    MPI_Bcast(&seed,4,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&N_p,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&N_c,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&N_gen,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&rect,1,MPI_INT,0,MPI_COMM_WORLD);
     if (my_rank==0){     
    	cout << "PROBLEMA DEL COMMESSO VIAGGIATORE" << endl;
    	cout << "Numero di città da visitare = " <<N_c << endl;
    	cout << "Numero di cromosomi nella popolazione = " << N_p << endl;
    	cout << "Numero di generazioni = " << N_gen << endl;
    	cout << "Distribuzione città = "<<rect<<endl;
     }
    rnd->SetRandom(seed,p1,p2);
    Mappa.Set_par(N_c,rect,rnd);
}

void Input(int my_rank){
    ifstream Read;
    Read.open("Primes");
    
    if(!Read.is_open()){
        cerr <<"Error: file Primes does not exist" <<endl;
        exit(-1);
    }

    int p1,p2,count=0;
    while(count<=my_rank){
        Read >> p1 >> p2;
        count++;
    }
   if (my_rank == 0){ 
	p1=2892;
	p2=2587;
   }
   if (my_rank == 1){
	p1=2893;
	p2=195;
   } 
   if(my_rank == 2){
	p1=2894;
	p2=2965;
   }
   if(my_rank==3){
	p1=2896;
	p2=1015;
   }
    Read.close();

    rnd->SetRandom(seed,p1,p2);

    for(int i=0; i<N_p; i++) pop.push_back(Chrom(N_c,rnd)); //carico tutti i cromosomi nella popolazione
    sort(pop.begin(), pop.end(), sort_fitness);            //li riordino in base al fitness
}

bool Check(){
    
    for (unsigned int i=0; i<pop.size(); i++){
        if(pop[i].Get_Element(0) != 1){
	    cout <<"non c'è 1 all'inizio"<<endl;  
            cout <<pop[i].Get_Element(0) <<endl;
            return false;
        }
	//cout << i<< endl;
        for (int j=0; j<pop[i].Get_Size()-1; j++){
	  
            for (int z=j+1; z<pop[i].Get_Size(); z++){
                 
		  if (pop[i].Get_Element(j)==pop[i].Get_Element(z))
		//	cout <<"percorso "<<i<< "  posizione1 "<<j<<" posizione2 "<<z<<endl;
                	return false;
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
            cross=(int)rnd->Rannyu(1.,N_c-1);   //punto di taglio dei cromosomi
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

bool sort_fitness( Chrom former, Chrom latter){              //mi permette di ordinare con questo criterio la pop
    return (former.Cost_fn_L1(Mappa) < latter.Cost_fn_L1(Mappa));
} 
