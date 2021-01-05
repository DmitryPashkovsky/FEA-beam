
#ifndef __beam_element_h_
#define __beam_element_h_
#include "def.h"



class beam
{
public:
	double l;
	double q;
	matrix K;
	matrix p;

	beam(double q_, double x1, double x2) : K(4, 4), p(4, 1)
	{
		q = q_;
		l = fabs(x2 - x1);
		K.M[0][0] = 12  / (l * l * l);
		K.M[0][1] = 6 / (l * l);
		K.M[0][2] = -12  / (l * l * l);
		K.M[0][3] = 6 / (l * l);

		K.M[1][0] = 6  / (l * l);
		K.M[1][1] = 4 / l;
		K.M[1][2] = -6  / (l * l);
		K.M[1][3] = 2  / l;

		K.M[2][0] = -12  / (l * l * l);
		K.M[2][1] = -6 / (l * l);
		K.M[2][2] = 12  / (l * l * l);
		K.M[2][3] = -6  / (l * l);

		K.M[3][0] = 6  / (l * l);
		K.M[3][1] = 2  / l;
		K.M[3][2] = -6  / (l * l);
		K.M[3][3] = 4 / l;
		///****
 		p.M[0][0] = q * l / 2;
		p.M[1][0] = q * l * l / 12;
		p.M[2][0] = q * l / 2;
		p.M[3][0] = -q * l * l / 12;
	}


	beam( const beam &B ) : K(4, 4), p(4, 1)
	{
		K = B.K;
		l = B.l;
		p = B.p;
		q = B.q;
	}

	~beam(void)
	{

	}
};

#endif