#ifndef ___Chrom_h__
#define ___Chrom_h__

#include"random.h"
#include"Maps.h"
#include<vector>

class Chrom{
	public:
		Chrom(int, Random* );
		Chrom( Chrom&,  Chrom&,int, Random*);
	//	~Chrom();

    double Cost_fn_L1(const Maps);

    int Get_Size(void)const;
    int Get_Element(int)const;
    
    void Print_conf(void);
    void Pair_permutation(void);
    void Shift(void);
    void m_permutation(void);
    void Inversion(void);
    void Save_conf(const Maps, int, int);
    void Save_conf_finale(const Maps);
    void Load_Element0(int var);
    void Load_Element(int var);

    protected:
		int m_N;
        std::vector<int> vect;
        Random *Trnd;
        double P_m=0.1;
};

#endif
