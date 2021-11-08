// fms_allocation.t.cpp - test allocation
#include <cassert>
#include "allocation.h"

using namespace port;


int allocation_test()
{
	{
		int n = 2;
		double ER[] = { 1.1,1.2 };
		double Cov[] = { 1, 2, 2, 13 };
		allocation p(n, ER, Cov, true);
		double x[4];
		double R, sigma;
		{
			R = ER[0];
			sigma = p.minimize(R, x);
			assert(fabs(1 - x[0]) < 1e-13);
			assert(fabs(0 - x[1]) < 1e-13);
			R = p.maximize(sigma, x);
			assert(fabs(ER[0] - R) < 1e-13);
		}
		{
			R = ER[1];
			sigma = p.minimize(R, x);
			assert(fabs(0 - x[0]) < 1e-13);
			assert(fabs(1 - x[1]) < 1e-13);
			R = p.maximize(sigma, x);
			assert(fabs(ER[1] - R) < 1e-13);
		}
		{
			R = ER[1];
			sigma = p.minimize(R, x);
			double l[] = { -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX };
			double u[] = { DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX };
			assert(fabs(0 - x[0]) < 1e-13);
			assert(fabs(1 - x[1]) < 1e-13);
			x[0] = 1;
			x[1] = 1;
			R = p.minimize(R, x, l, u);
			assert(fabs(0 - x[0]) < 1e-13);
			assert(fabs(1 - x[1]) < 1e-13);
		}
		{
			double l[] = { 1, 1, -DBL_MAX, -DBL_MAX };
			R = ER[0];
			R = p.maximize(sigma, x);
			sigma = p.minimize(R, x, l);
			//assert(fabs(1 - x[0]) < 1e-13);
			//assert(fabs(0 - x[1]) < 1e-13);
			//assert(fabs(ER[0] - R) < 1e-13);
		}
	}

	return 0;
}
int allocation_test_ = allocation_test();