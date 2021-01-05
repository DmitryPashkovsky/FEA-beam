
#ifndef __init_h_
#define __init_h_
#include "def.h"


class init
{
private:
	std::vector<double> n_dof;
	std::vector<int> x_dof;
	double x_start, x_end, q;
	int n_elems;
public:

	init( std::string fname );
	void info( void );
	std::vector<double> get_ndf( void );
	std::vector<int> get_xdf( void );
	double get_start( void );
	double get_end( void );
	double get_q(void);
	int get_n_el(void);
	~init( void );
}; 
#endif // !__init_h_
