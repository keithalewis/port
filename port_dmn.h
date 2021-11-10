// port_dmn.h - PORT optimization routines
#pragma once
#include "port.h"
#include <vector>

namespace port {

	// minimize unconstrained objective function
	struct dmn {
		int n, liv, lv;
		std::vector<int> iv;
		std::vector<double> v, d;
		std::vector<double> b; // optional lower and upper bounds
		
		// initialize upper and lower bounds
		void bounds()
		{
			b.resize(2 * n);
			double x = DBL_MAX;
			upper(&x, 0);
			x = -x;
			lower(&x, 0);
		}

		// reverse communication functions
		typedef void(*FGH)(int* N, double* X, int* NF, double* F, int* UI, double* UR, void* UF);

		dmn(int n)
			: n(n), liv(59 + 3*n), lv(78 + n * (n + 15)), iv(liv), v(lv), d(n, 1)
		{ }
		dmn(const dmn&) = delete;
		dmn& operator=(const dmn&) = delete;
		~dmn()
		{ }

		// set lower bounds
		void lower(const double* l, size_t stride = 1)
		{
			if (b.size() == 0) {
				bounds();
			}
			for (int i = 0; i < n; ++i) {
				b[2 * i] = l[i * stride];
			}
		}
		// set upper bounds
		void upper(const double* u, size_t stride = 1)
		{
			if (b.size() == 0) {
				bounds();
			}
			for (int i = 0; i < n; ++i) {
				b[2 * i + 1] = u[i * stride];
			}
		}

		// function
		RETURN_CODE solve(double* x, FGH f, int* ui = 0, double* ur = 0, void* dummy = 0)
		{
			if (b.size() == 0) {
				DMNF(&n, d.data(), x, f, iv.data(), &liv, &lv, v.data(), ui, ur, dummy);
			}
			else {
				DMNFB(&n, d.data(), x, b.data(), f, iv.data(), &liv, &lv, v.data(), ui, ur, dummy);
			}

			return (RETURN_CODE)iv[0];
		}

		// function and gradient
		RETURN_CODE solve(double* x, FGH f, FGH g, int* ui = 0, double* ur = 0, void* dummy = 0)
		{
			if (b.size() == 0) {
				DMNG(&n, d.data(), x, f, g, iv.data(), &liv, &lv, v.data(), ui, ur, dummy);
			}
			else {
				DMNGB(&n, d.data(), x, b.data(), f, g, iv.data(), &liv, &lv, v.data(), ui, ur, dummy);
			}

			return (RETURN_CODE)iv[0];
		}

		/*
		// function, gradient, and hessian
		RETURN_CODE solve(double* x, FGH f, FGH g, FGH h, int* ui = 0, double* ur = 0, void* dummy = 0)
		{
			if (b.size() == 0) {
				DMNH(&n, d.data(), x, f, g, h, iv.data(), &liv, &lv, v.data(), ui, ur, dummy);
			}
			else {
				DMNHB(&n, d.data(), x, b.data(), f, g, h, iv.data(), &liv, &lv, v.data(), ui, ur, dummy);
			}

			return (RETURN_CODE)iv[0];
		}
		*/
	};

} // namespace port
