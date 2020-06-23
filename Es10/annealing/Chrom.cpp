#include<iostream>
#include<algorithm>
#include"math.h"
#include"Chrom.h"
#include <fstream>
#include <ostream>
#include <iomanip>
using namespace std;

//costruttore semplice
Chrom::Chrom(){}

//costruttore inizializzazione
Chrom::Chrom(int N_all, Random *Srnd){
        m_N=N_all;
        Trnd=Srnd;
        for(int j=0;j<N_all;j++) vect.push_back(j+1);
        random_shuffle(vect.begin()+1,vect.end());
}
      
//distruttore
Chrom::~Chrom(){}

//metodo per calcolare la funzione costo
double Chrom::Cost_fn_L1(const Maps Cities){	
        double L1=Cities.Dist(vect.at(0),vect.at(m_N-1));
        for (unsigned int i=0; i<vect.size()-1; i++) L1+=Cities.Dist(vect.at(i),vect.at(i+1));
        return L1;
}

//************************************* METODI PER MUTAZIONI ****************************
void Chrom::Pair_permutation(void){
    if(P_m>Trnd->Rannyu()){
        int ind_1=int(Trnd->Rannyu(1.,m_N));
        int ind_2=int(Trnd->Rannyu(1.,m_N));
        while(ind_1==ind_2)	ind_2=int(Trnd->Rannyu(1.,m_N));
        swap(vect[ind_1],vect[ind_2]);
    }
}


void Chrom::Shift(void){
    if(P_m>Trnd->Rannyu()){
        int ind=int(Trnd->Rannyu(1.,m_N-1));
        rotate(vect.begin()+1,vect.end()-ind,vect.end());
	}
} //Rotates the order of the elements in the range [first,last), in such a way that
  //the element pointed by middle becomes the new first element.


void Chrom::m_permutation(void){
   
    if(P_m>Trnd->Rannyu()){                   //scelgo due intervalli di m contigui
		int m=int(Trnd->Rannyu(1.,m_N/2.));
		int n=int(Trnd->Rannyu(m_N/2,m_N-m));
        int ind=int(Trnd->Rannyu(1.,n-m));
        vector<int> vect_copy (m_N);
	
        copy(vect.begin()+ind,vect.begin()+ind+m,vect_copy.begin()+n);
        copy(vect.begin()+n,vect.begin()+n+m,vect_copy.begin()+ind);
        swap_ranges(vect_copy.begin()+ind,vect_copy.begin()+ind+m,vect.begin()+ind);
        swap_ranges(vect_copy.begin()+n,vect_copy.begin()+n+m,vect.begin()+n);
	}
}


void Chrom::Inversion(void){
    if(P_m>Trnd->Rannyu()){
        int m=int(Trnd->Rannyu(0.,m_N));
        int ind=int(Trnd->Rannyu(1.,m_N-m));
        reverse(vect.begin()+ind,vect.end()-m);
	}
}
/* ************************************************************************************************* */

int Chrom::Get_Size(void)const{
	return vect.size();
}

int Chrom::Get_Element(int i)const{
	return vect.at(i);
}

void Chrom::Print_conf(void){
  for(int i=0;i<m_N;i++) cout<<vect[i]<<" ";
  cout<<"\n";
}

void Chrom::Save_conf(const Maps Cities){
	ofstream Conf_cit("Maps");
    for(int i=0;i<m_N;i++)	Conf_cit<<setprecision(9)<<i<<"\t"<<Cities.GetPos(vect.at(i)).x<<"\t"<<Cities.GetPos(vect.at(i)).y<<endl;
        Conf_cit<<setprecision(9)<<m_N<<"\t"<<Cities.GetPos(vect.at(0)).x<<"\t"<<Cities.GetPos(vect.at(0)).y<<endl;
        Conf_cit.close();
}





