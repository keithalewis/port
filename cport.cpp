// cport.cpp : Test PORT library
#include <cassert>
#include "dmn.h"


int d1mach_test()
{
	double d;
	int i = 1;
	d = D1MACH(&i);

	return 0;
}
int d1mach_test_ = d1mach_test();

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

	return 0;
}
//int dmnf_test_ = dmnf_test();

int dmnfb_test()
{
	port::RETURN_CODE ret;
	double x[] = { 1,2 };
	double l[] = { 0.5, 0.5 };
	double u[] = { DBL_MAX, DBL_MAX };
	port::dmn p(2);

	auto f = [](int* N, double* X, int* /*NF*/, double* F, int* /*UI*/, double* /*UR*/, void* /*UF*/) {
		*F = 0;
		for (int i = 0; i < *N; ++i) {
			*F += X[i] * X[i];
		}
	};
	p.lower(l);
	p.upper(u);
	ret = p.solve(x, f);

	return 0;
}
//int dmnfb_test_ = dmnfb_test();

int main()
{
	return 0;
}