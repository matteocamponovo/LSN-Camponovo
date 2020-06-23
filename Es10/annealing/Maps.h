#ifndef ___Maps_h__
#define ___Maps_h__
#include"random.h"
#include<vector>

struct Pos{
  double x;
  double y;
};

const double pi=3.1415927;
class Maps{
    public:
        Maps();
        Maps(int Npts,int dis, Random *Trnd);
        ~Maps();

    void Set_par(int Npts,int dis, Random *Trnd);
    Pos GetPos(int) const;
    double Dist(int, int) const;

  protected:
        int m_N;
        Random *Trnd;
        std::vector<Pos> Cities;
    
};

#endif
