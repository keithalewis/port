// port_dmn.h - PORT optimization routines
#pragma once
#pragma warning(disable: 26451)
#include "port.h"
#include <vector>

namespace port {

	// minimize objective function
	struct dmn {
		int n, liv, lv;
		std::vector<int> iv;
		std::vector<double> v, d, b;
	
		// initialize upper and lower bounds
		void bounds()
		{
			b.resize(2 * n);
			for (int i = 0; i < n; ++i) {
				b[2*i] = -DBL_MAX;
				b[2*i + 1] = DBL_MAX;
			}
		}

		// reverse communication functions
		typedef void(*FGH)(int* N, double* X, int* NF, double* F, int* UI, double* UR, void* UF);

		dmn(int n)
			: n(n), liv(59 + 3*n), iv(liv), lv(78 + (n * (n + 27)) / 2), v(lv), d(n, 1)
		{
			//int alg = 2; // general optimization
			//DIVSET(&alg, &iv[0], &liv, &lv, &v[0]);
		}
		dmn(const dmn&) = delete;
		dmn& operator=(const dmn&) = delete;
		~dmn()
		{ }

		// control printing to console
		void print(int i = 2) // stdout
		{
			iv[20] = i ? I1MACH(&i) : 0;
		}

		// initial step size
		void step0(double LMAX0 = 1)
		{
			v[34] = LMAX0;
		}

		// set lower bounds
		void lower(double l)
		{
			if (b.size() == 0) {
				bounds();
			}
			for (int i = 0; i < n; ++i) {
				b[2 * i] = l;
			}
		}
		void lower(const double* l)
		{
			if (b.size() == 0) {
				bounds();
			}
			for (int i = 0; i < n; ++i) {
				b[2 * i] = l[i];
			}
		}
		// set upper bounds
		void upper(double u)
		{
			if (b.size() == 0) {
				bounds();
			}
			for (int i = 0; i < n; ++i) {
				b[2 * i + 1] = u;
			}
		}
		void upper(const double* u)
		{
			if (b.size() == 0) {
				bounds();
			}
			for (int i = 0; i < n; ++i) {
				b[2 * i + 1] = u[i];
			}
		}

		// function
		RETURN_CODE solve(double* x, FGH f, int* ui = 0, double* ur = 0, void* uf = 0)
		{
			if (b.size() == 0) {
				DMNF(&n, &d[0], &x[0], f, &iv[0], &liv, &lv, &v[0], ui, ur, uf);
			}
			else {
				DMNFB(&n, &d[0], x, &b[0], f, &iv[0], &liv, &lv, &v[0], ui, ur, uf);
			}

			return (RETURN_CODE)iv[0];
		}

		// function and gradient
		RETURN_CODE solve(double* x, FGH f, FGH g, int* ui = 0, double* ur = 0, void* uf = 0)
		{
			if (b.size() == 0) {
				DMNG(&n, d.data(), x, f, g, iv.data(), &liv, &lv, v.data(), ui, ur, uf);
			}
			else {
				DMNGB(&n, d.data(), x, b.data(), f, g, iv.data(), &liv, &lv, v.data(), ui, ur, uf);
			}

			return (RETURN_CODE)iv[0];
		}

		/*
		// function, gradient, and hessian
		RETURN_CODE solve(double* x, FGH f, FGH g, FGH h, int* ui = 0, double* ur = 0, void* uf = 0)
		{
			struct GH {
				FGH g, h
			} GH = { g, h };
			auto gh[] = [](...) {
				GH* pgh = *uf;
				...
			}
			if (b.size() == 0) {
				DMNH(&n, d.data(), x, f, gh, iv.data(), &liv, &lv, v.data(), ui, ur, uf);
			}
			else {
				DMNHB(&n, d.data(), x, b.data(), f, gh, iv.data(), &liv, &lv, v.data(), ui, ur, uf);
			}

			return (RETURN_CODE)iv[0];
		}
		*/
	};

} // namespace port
