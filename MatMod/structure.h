#ifndef __structure_h_
#define __structure_h_
#include "def.h"

class bsystem 
{
public:
	int N; 
	int n;
	int num_dof;
	std::vector<int> idxs;
	std::vector<double> i_dofs;
	std::vector<double> nodes;
	std::vector<beam> els;
	matrix K;
	matrix Q;
	matrix F;
	matrix p;
	matrix dsp;
	matrix s;
	matrix r;

	std::vector<double> Mgl, Ngl, DSgl, GRID;

	bsystem(init &F );
	void create_global_matrix( void );
	void create_right_vector( void );
	std::vector<double> error(void);
	void create_connection_matrix(void);
	void total_condence(void);
	void solve( void );
	void condense(int k);
	void save_result(void);
	void element_diagram_M_N( int idx_elem );
	void all_diagram_M_N(void);
	~bsystem( void );
};

#endif