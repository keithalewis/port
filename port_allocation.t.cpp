// port_allocation.t.cpp - test allocation
#include <cassert>
#include "port_allocation.h"

using namespace port;

int allocation_test()
{
	{
		int n = 2;
		double ER[] = { 1.1,1.2 };
		double Cov[] = { 1, 2, 2, 13 };
		
		double L[3];
		port::packl(2, Cov, L);
		assert(L[0] == 1 && L[1] == 2 && L[2] == 13);
		cholesky(2, L, L);
		assert(L[0] == 1 && L[1] == 2 && L[2] == 3);

		allocation p(n, ER, Cov, false);
		//p.print(0);

		// L = [1 0; 2 3], L^-1 = [1 0; -2/3 1/3]
		// V = L L', V^-1 = L'^-1 L^-1
		// V^-1 1 = [1 -2/3; 0 1/3] [1 -1/3] = [1 + 2/9, -1/9]
		assert(fabs(p.V_1[0] - 11. / 9) < 1e-13);
		assert(fabs(p.V_1[1] + 1. / 9) < 1e-13);

		double x[4];
		double R, sigma;

		{
			x[0] = 0;
			x[1] = 0;
			assert(0 == quad(2, x, L));
			x[0] = 1;
			assert(1 == quad(2, x, L));
			x[1] = 1;
			assert(1 + 2 * 2 + 13 == quad(2, x, L));
			sigma = p.fmin(0, x, 0, 0);
			sigma = quad(2, x, L) / 2;
		}

		{
			R = ER[0];
			sigma = p.minimum(R, x);
			assert(fabs(1 - x[0]) < 1e-12);
			assert(fabs(0 - x[1]) < 1e-12);
			R = p.maximum(sigma, x);
			assert(fabs(ER[0] - R) < 1e-12);
		}
		{
			R = ER[1];
			sigma = p.minimum(R, x);
			assert(fabs(0 - x[0]) < 1e-12);
			assert(fabs(1 - x[1]) < 1e-12);
			R = p.maximize(sigma, x);
			assert(fabs(ER[1] - R) < 1e-12);
		}
		{
			double lu = 0;// DBL_MAX;
			double l[] = { -lu, -lu, -DBL_MAX, -DBL_MAX };
			double u[] = { DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX };

			R = ER[1];
			sigma = p.minimum(R, x);
			assert(fabs(0 - x[0]) < 1e-11);
			assert(fabs(1 - x[1]) < 1e-11);

			x[0] += .001;
			x[1] -= .001;
			//p.lower(l);
			//p.upper(u);
			sigma = p.minimize(R, x);
		}
		/*
		{
			double l[] = { .1, .1, -DBL_MAX, -DBL_MAX };
			R = ER[0];
			R = p.maximum(sigma, x);
			R = p.maximize(sigma, x, l);
			R = R;
		}
		*/
	}

	return 0;
}
//int allocation_test_ = allocation_test();