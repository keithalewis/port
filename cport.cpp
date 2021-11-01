// cport.cpp : Test PORT library
#include <cassert>
#include "port.h"

extern "C" double D1MACH(int*);

int d1mach_test()
{
	double d;
	int i = 1;
	d = D1MACH(&i);

	return 0;
}
int d1mach_test_ = d1mach_test();

void f(int* N, double* X, int* /*NF*/, double* F, int* /*UI*/, double* /*UR*/, void* /*UF*/) {
	*F = 0;
	for (int i = 0; i < *N; ++i) {
		*F += X[i] * X[i];
	}
}

int drmn_test()
{
	double x[] = { 1,2 };
	port::drmn p(2, x);
#if 0
	auto f = [](int* N, double* X, int* /*NF*/, double* F, int* /*UI*/, double* /*UR*/, void* /*UF*/) {
		*F = 0;
		for (int i = 0; i < *N; ++i) {
			*F += X[i] * X[i];
		}
	};
#endif
	p.solve(f);

	return 0;
}
int drmn_test_ = drmn_test();

int main()
{
	return 0;
}