#include "structure.h"

bsystem::bsystem( init &F ) : N(F.get_n_el()), n(N + 1), num_dof(2 * (N + 1)), K(num_dof, num_dof), p(num_dof, 1), Q(num_dof, 4 * N), F(4 * N, 4 * N), dsp(num_dof,1), s(4 * N, 1), r(num_dof, 1)
{
	double q = F.get_q(), h = fabs(F.get_start() - F.get_end()) / n;
	std::cout << q << "\n";
	for (int i = 0; i < n; i++)
		nodes.push_back(i * h);
	
	idxs = F.get_xdf();
	i_dofs = F.get_ndf();


	for (int i = 0; i < N; i++)
		els.push_back(beam(q, nodes[i], nodes[i + 1]));

}

void bsystem::create_global_matrix( void )
{
	for (int i = 0; i < N; i++)
	{
		for (int k = 0; k < 4; k++)
		    for (int j = 0; j < 4; j++)
			     F.M[4 * i + k][4 * i + j] = els[i].K.M[k][j];
	}

	K = Q * F * Q.trans();

}
void  bsystem::condense(int k)
{
	for (int i = 0; i < K.m; i++)
	{
		K.M[k][i] = 0;
		K.M[i][k] = 0;
	}

	K.M[k][k] = 1;
	p.M[k][0] = 0;
}

void  bsystem::total_condence(void)
{
	for (int i = 0; i < int(idxs.size()); i++)
	{
		if (i_dofs[i] == 2)
		{
			this->condense(2 * idxs[i]);
			this->condense(2 * idxs[i] + 1);
		}
		else if (i_dofs[i] == 1.1)
			this->condense(2 * idxs[i]);
		else if (i_dofs[i] == 1.2)
			this->condense(2 * idxs[i] + 1);
		else
		{
		}

	}
}

void bsystem::create_right_vector( void )
{
	for (int i = 0; i < N; i++)
	{
		int bdof = 2 * i;
		
		p.M[bdof][0] -= els[i].p.M[0][0];
		p.M[bdof + 1][0] -= els[i].p.M[1][0];
		p.M[bdof + 2][0] -= els[i].p.M[2][0];
		p.M[bdof + 3][0] -= els[i].p.M[3][0];
	}
	//p.print();
}

void bsystem::create_connection_matrix( void )
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < 4; j++)
			Q.M[2 * i + j ][4 * i + j] = 1;
}

void bsystem::save_result(void)
{
	std::ofstream out1("displacement.txt");
	std::ofstream out2("grid.txt");

	std::ofstream out3("Moment.txt");
	std::ofstream out4("Forces.txt");
	std::ofstream out5("GRIDDIAGRAMM.txt");


	if (!out1.is_open())
	{
		std::cout << "displacement.txt is not opened !\n";
		return;
	}
	if (!out2.is_open())
	{
		std::cout << "grid.txt is not opened !\n";
		return;
	}
	if (!out3.is_open())
	{
		std::cout << "Moment.txt is not opened !\n";
		return;
	}
	if (!out4.is_open())
	{
		std::cout << "Forces.txt is not opened !\n";
		return;
	}
	if (!out5.is_open())
	{
		std::cout << "GRIDDIAGRAMM.txt is not opened !\n";
		return;
	}

	std::vector<double> dspp;
	for (int i = 0; i < n; i++)
		dspp.push_back(dsp.M[2 * i][0]);

	std::copy(dspp.begin(), dspp.end(), std::ostream_iterator<double>(out1, ","));
	std::copy(nodes.begin(), nodes.end(), std::ostream_iterator<double>(out2, ","));

	std::copy(Mgl.begin(), Mgl.end(), std::ostream_iterator<double>(out3, ","));
	std::copy(Ngl.begin(), Ngl.end(), std::ostream_iterator<double>(out4, ","));
	std::copy(GRID.begin(), GRID.end(), std::ostream_iterator<double>(out5, ","));
	
	out1.close();
	out2.close();
	out3.close();
	out4.close();
	out5.close();
	
}

void bsystem::solve(void)
{
	SLAE FEA(K, p);

	FEA.LDU_decomposition();
	FEA.reverse_substitution();
	dsp = FEA.x;
	s = F * Q.trans() * dsp;
	r = Q * s - p;

}


void bsystem::element_diagram_M_N( int idx_elem )
{
	int num_pts = 10;
	double l = els[idx_elem].l, q = els[idx_elem].q, h = l / num_pts, l2 = l * l, l3 = l2 * l;
	std::vector<double> X;
	std::vector<double> V, M, d;
	//double vt, vr;
	double d1 = dsp.M[2 * (idx_elem)][0],
		   d2 = dsp.M[2 * (idx_elem) + 1][0],
		   d3 = dsp.M[2 * (idx_elem) + 2][0],
		   d4 = dsp.M[2 * (idx_elem) + 3][0];

	std::cout << d1 << " " << d2 << " " << d3 << " " << d4 << "\n";
	for (int i = 0; i <= 10; i++)
		X.push_back(i * h);

	double Va = q * l / 2, Ra = -q * l * l / 24;

	// от нагрузки
	for (int i = 0; i <= 10; i++)
	{
		V.push_back(Va * (1 - 2 * X[i] / l));
		M.push_back(-Va * l / 6 * (6 * X[i] / l - 6 * X[i] * X[i] / (l * l) - 1));

	}

	//// 2
	double vt = 12 / (l * l * l) * d1;
	for (int i = 0; i <= 10; i++)
	{
		V[i] += vt;
		M[i] += (0.5 * l * vt) * (1 - 2 * X[i] / l);
	}
	 
	double vr = 6 / (l * l) * d2;
	for (int i = 0; i <= 10; i++)
	{
		V[i] += vr;
		M[i] += (vr * l / 3) * (2 - 3 * X[i] / l);
	}


	// 4
	vt = 12 / (l * l * l) * d3;
	for (int i = 0; i <= 10; i++)
    {
		V[i] -= vt;
		M[i] += (0.5 * l * vt) * (1 - 2 * (l - X[i]) / l);
    }

	// 5
	vr = 6 / (l * l) * d4;
	for (int i = 0; i <= 10; i++)
	{
		V[i] += vr;
		M[i] -= (vr * l / 3) * (2 - 3 * (l - X[i]) / l);
	}

	for (int i = 0; i <= 10; i++)
	{
		Mgl.push_back(M[i]);
		Ngl.push_back(V[i]);
	}
}


void bsystem::all_diagram_M_N(void)
{
	for (int i = 0; i < int(els.size()); i++)
	{
		this->element_diagram_M_N(i);
		
		if (i != int(els.size()) - 1)
		{
			Mgl.pop_back();
			Ngl.pop_back();
		}
	}
}


bsystem::~bsystem(void)
{
}