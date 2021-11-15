// cport.cpp : Test PORT library
#include <cassert>
#include "port_dmn.h"

int d1mach_test()
{
	double d;
	int i = 1;
	d = D1MACH(&i);

	return 0;
}
int d1mach_test_ = d1mach_test();

int quad_test()
{
	{
		double L[] = { 1, 2, 3 };
		// LL' = [1 2; 2 13]
		double x[] = { 0, 0 };
		assert(0 == port::quad(2, x, L));
		x[0] = 1;
		assert(1 == port::quad(2, x, L));
		x[1] = 1;
		assert(1 + 4 + 13 == port::quad(2, x, L));

	}

	return 0;
}
int quad_test_ = quad_test();

int dmnf_test()
{
	port::RETURN_CODE ret;
	double x[] = { 1,2 };
	port::dmn p(2);

	auto f = [](int* N, double* X, int* /*NF*/, double* F, int* /*UI*/, double* /*UR*/, void* /*UF*/) {
		*F = 0;
		for (int i = 0; i < *N; ++i) {
			*F += X[i] * X[i];
		}
	};
	ret = p.solve(x, f);
	assert(port::RETURN_CODE::ABS_CONV == ret);

	return 0;
}
//int dmnf_test_ = dmnf_test();

int dmnfb_test()
{
	port::RETURN_CODE ret;
	double x[] = { 3, 4 };
	double l[] = { -.5, .5 };
	double u[] = { DBL_MAX, DBL_MAX };
	port::dmn p(2);

	auto f = [](int* N, double* X, int* /*NF*/, double* F, int* /*UI*/, double* /*UR*/, void* /*UF*/) {
		*F = 0;
		for (int i = 0; i < *N; ++i) {
			*F += (i + 1) * X[i] * (X[i]);
		}
	};
	p.lower(l);
	p.upper(u);
	ret = p.solve(x, f);
	assert(port::RETURN_CODE::REL_CONV == ret);

	return 0;
}
//int dmnfb_test_ = dmnfb_test();

int main()
{
	return 0;
}