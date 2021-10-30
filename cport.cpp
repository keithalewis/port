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

int drmn_test()
{
	double x[] = { 1,2 };
	port::drmn(2, x);

	return 0;
}
int drmn_test_ = drmn_test();

int main()
{
	return 0;
}