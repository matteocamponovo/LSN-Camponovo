#include"Maps.h"
#include<iostream>
#include<stdexcept>
#include<stdlib.h>
#include<sstream>
#include"math.h"

using namespace std;

//classe maps.
//il costruttore genera un vettore di strutture Pos
//i membri di pos sono la posizione (𝑥,𝑦) di ogni città
//la classe contiene anche un metodo per il calcolo delle
//distanze tra due città della mappa, che chiamo Dist

//Ncittà è il numero di città che andrò a posizionare
//disponi è per scegliere se posizionare su quadrato o circonferenza
//Trnd è un puntatore alla classe Random

Maps::Maps(){}

Maps::Maps(int Npts, int disp, Random *Trnd){
    m_N=Npts;
    Cities.resize(m_N);
 double theta;
    for(int j=0;j<m_N; j++){
    if(disp==0){//i punti sono generati su un quadrato lato 1
      Cities[j].x=Trnd->Rannyu(-1.,1.);
      Cities[j].y=Trnd->Rannyu(-1.,1.);
    }else{
      theta=Trnd->Rannyu(0.,2.*pi);
      Cities[j].x=cos(theta);
      Cities[j].y=sin(theta);
        }
    }
}

void Maps::Set_par(int Npts, int disp, Random *Trnd){
  m_N=Npts;
  
    Cities.resize(m_N);

 double theta;
    for(int j=0;j<m_N; j++){
    if(disp==0){//i punti sono generati su un quadrato lato 1
      Cities[j].x=Trnd->Rannyu(-1.,1.);
      Cities[j].y=Trnd->Rannyu(-1.,1.);
    }else{
      theta=Trnd->Rannyu(0.,2.*pi);
      Cities[j].x=cos(theta);
      Cities[j].y=sin(theta);
        }
    }
}
      

Pos Maps::GetPos(int i) const{
    if (i<=m_N) return Cities[i-1];
    else{
     cerr<<"l'intero assegnato e' maggiore della lunghezza della struct"<<endl;
    exit(-1);
    }
}


double Maps::Dist(int i, int l) const{
    double ds =pow(Cities[i-1].x-Cities[l-1].x,2.)+pow(Cities[i-1].y-Cities[l-1].y,2.);
  return sqrt(ds);
}



Maps::~Maps(){}








