// port_dmn.t.cpp - test DMN* routines
#include <cassert>
#include "port_dmn.h"
#include <numeric>

using namespace port;

int dmnf_test()
{
	// x[0] + x[1] - x[2]*(x[0]^2 + x[1]^2 - R^2)/2
	auto f = [](int* N, double* X, int* NF, double* F, int* /*UI*/, double* UR, void* UF) {
		*N = *N;
		*NF = *NF;
		dmn* p = (dmn*)UF; p = p;
		double vf, vf0, vdstnrm, vlmax0;
		vf = p->v[V::F-1];
		vf0 = p->v[V::F0-1];
		vdstnrm = p->v[V::DSTNRM-1];
		vlmax0 = p->v[V::LMAX0-1];
		double vrfctol;
		vrfctol = p->v[RFCTOL-1];
		double vsctol;
		vsctol = p->v[V::SCTOL-1];
		double veta0;
		veta0 = p->v[V::ETA0-1];
		p->v[V::ETA0 - 1] /= 10;

		double R = *UR;
		*F = X[0] + X[1] - X[2] * (X[0] * X[0] + X[1] * X[1] - R * R)/2;
		//*F = (X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
		*F = *F;
	};
	auto g = [](int* N, double* X, int* NF, double* G, int* /*UI*/, double* UR, void* UF) {
		*N = *N;
		*NF = *NF;
		dmn* p = (dmn*)UF; p = p;
		double R = *UR;

		G[0] = 1 - X[2] * X[0];
		G[1] = 1 - X[2] * X[1];
		G[2] = (R * R - (X[0] * X[0] + X[1] * X[1]))/2;
		G[2] = G[2];
	};
	auto h = [](int* N, double* X, int* NF, double* H, int* /*UI*/, double* UR, void* UF) {
		*N = *N;
		*NF = *NF;
		dmn* p = (dmn*)UF; p = p;
		double R = *UR;

		H[0] = -X[2];
		H[1] = 0;
		H[2] = -X[2];
		H[3] = -X[0];
		H[4] = -X[1];
		H[5] = 0;
		H[0] = H[0];
	};


	port::dmn p(3);
	//p.v[V::ETA0-1] = 1e-9;
	double r = sqrt(1);
	double eps = 0.0001;
	double x[] = { -r/sqrt(2) + eps,-r/sqrt(2) + eps, -sqrt(2)/r + eps };
	//p.step0(1e-4);
	p.print(0);
	RETURN_CODE ret;
	ret = p.solve(x, f, 0, &r, &p);
	//ret = p.solve(x, f, g, 0, &r, &p);
	//assert(RETURN_CODE::REL_CONV == ret);
	double err;
	err = x[0] + r / sqrt(2);
	assert(fabs(err) < 1e-5);
	err = x[1] + r/sqrt(2);
	assert(fabs(err) < 1e-5);
	err = x[2] + sqrt(2) / r;
	assert(fabs(err) < 1e-5);

	return 0;
}
//int dmnf_test_ = dmnf_test();

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
		*F = *F;
	};
	p.lower(l);
	p.upper(u);
	ret = p.solve(x, f);
	//assert(port::RETURN_CODE::X_CONV == ret);

	return 0;
}
//int dmnfb_test_ = dmnfb_test();
