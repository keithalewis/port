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
		assert(1 + 2 * 2 + 13 == port::quad(2, x, L));
	}

	return 0;
}
int quad_test_ = quad_test();

int main()
{
	return 0;
}