// port_dmn.t.cpp - test DMN* routines
#include <cassert>
#include "port_dmn.h"
#include <numeric>

using namespace port;

double f1(int n, const double* x, const double* a = nullptr)
{
	double y = 0;

	for (int i = 0; i < n; ++i) {
		y += (a ? a[i] : 1) * x[i] * x[i];
	}

	return y;
}

int dmnf_test()
{
	port::RETURN_CODE ret;
	double x[] = { 1,1,0 };
	port::dmn p(3);

	// min F(X) = x[0] + x[1] - lambda(x[0]^2 + x[1]^2 - R^2)
	// DF = (1 - 2lambda x[0], 1 - 2lambda x[1], -(x[0]^2 + x[1]^2 - R^2))
	auto f = [](int* N, double* X, int* NF, double* F, int* /*UI*/, double* UR, void* UF) {
		int n = *N;
		*NF = *NF;
		port::dmn* p = (port::dmn*)UF; p = p;
		double vf, vf0;
		vf = p->v[V::F];
		vf0 = p->v[V::F0];
		*F = std::accumulate(X, X + n - 1, 0.);
		double R = *UR;
		*F -= X[n-1] * (f1(n - 1, X) - R*R);
	};
	double r = 1;
	//p.step0(1e-4);
	ret = p.solve(x, f, 0, &r, &p);
	assert(port::RETURN_CODE::REL_CONV == ret);
	double err;
	err = x[0] + sqrt(r / 2);
	assert(fabs(err) < 1e-5);
	err = x[1] + sqrt(r / 2);
	assert(fabs(err) < 1e-5);

	return 0;
}
int dmnf_test_ = dmnf_test();

int dmnfb_test()
{
	port::RETURN_CODE ret;

	double x[] = { 1, 1 };
	double l[] = { -.5, .5 };
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
	//assert(port::RETURN_CODE::X_CONV == ret);

	return 0;
}
//int dmnfb_test_ = dmnfb_test();
