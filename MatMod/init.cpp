#include "init.h"

init::init(std::string fname)
{
	std::string l;

	std::ifstream input(fname);
	if (input.is_open())
	{
		input >> x_start;
		input >> x_end;
		input >> q;
		input >> n_elems;

		while (!input.eof())
		{
			double x, f;
			input >> x;
			input >> f;
			x_dof.push_back(x);
			n_dof.push_back(f);
		}
	}
	input.close();
}

void init::info( void )
{
	std::cout << "INFORMATION ABOUT SYSTEM" << "\n";
	std::cout << "start: " << x_start << "\n";
	std::cout << "end: " << x_end << "\n";
	std::cout << "Fixed degree of freedom:\n";
	for (int i = 0; i < int(x_dof.size()); i++)
		std::cout << " " << x_dof[i] << " " << n_dof[i] << "\n";
}

std::vector<double> init::get_ndf( void )
{
	return n_dof;
}

std::vector<int> init::get_xdf( void )
{
	return x_dof;
}

double init::get_start( void )
{
	return x_start;
}

double init::get_end( void )
{
	return x_end;
}
double init::get_q(void)
{
	return q;
}

int init::get_n_el(void)
{
	return n_elems;
}

init::~init( void )
{

}
/* END OF FILE 'init.cpp' */