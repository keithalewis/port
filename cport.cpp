// cport.cpp : Test PORT library
#include <cassert>
#include "dmn.h"
#include "dmnb.h"

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
	double x[] = { 1,2 };
	port::dmnf p(2, x);

	auto f = [](int* N, double* X, int* /*NF*/, double* F, int* /*UI*/, double* /*UR*/, void* /*UF*/) {
		*F = 0;
		for (int i = 0; i < *N; ++i) {
			*F += X[i] * X[i];
		}
	};
	p.solve(f);

	return 0;
}
int dmnf_test_ = dmnf_test();

int dmnfb_test()
{
	double x[] = { 1,2 };
	double l[] = { 0.5, 0.5 };
	double u[] = { DBL_MAX, DBL_MAX };
	port::dmnfb p(2, x, l, u);

	auto f = [](int* N, double* X, int* /*NF*/, double* F, int* /*UI*/, double* /*UR*/, void* /*UF*/) {
		*F = 0;
		for (int i = 0; i < *N; ++i) {
			*F += X[i] * X[i];
		}
	};
	p.solve(f);

	return 0;
}
int dmnfb_test_ = dmnfb_test();

int main()
{
	return 0;
}