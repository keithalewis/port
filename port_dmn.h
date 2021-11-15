// port_dmn.h - PORT optimization routines
#pragma once
#pragma warning(disable: 26451)
#include "port.h"
#include <vector>

namespace port {

	// minimize unconstrained objective function
	struct dmn {
		int n, liv, lv;
		std::vector<int> iv;
		std::vector<double> v, d, b;
		
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
			: n(n), liv(59), iv(liv), lv(71), v(lv), d(n, 1)
		{
			//int alg = 2; // general optimization
			//DIVSET(&alg, &iv[0], &liv, &lv, &v[0]);
		}
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
		RETURN_CODE solve(double* x, FGH f, int* ui = 0, double* ur = 0, void* uf = 0)
		{
			if (b.size() == 0) {
				double fx;
				DRMNF(&d[0], &fx, &iv[0], &liv, &lv, &n, &v[0], x);
				liv = iv[43]; // IVLAST
				iv.resize(liv);
				lv = iv[44]; // VLAST
				v.resize(lv);

				iv[0] = 0;
				DRMNF(&d[0], &fx, &iv[0], &liv, &lv, &n, &v[0], x);

				while (iv[0] <= 2) {
					int nf = iv[5]; // NFCALL
					f(&n, x, &nf, &fx, ui, ur, uf);
					DRMNF(&d[0], &fx, &iv[0], &liv, &lv, &n, &v[0], x);
					if (nf <= 0) {
						iv[1] = 1; // TOOBIG
					}
				}
			}
			else {
				double fx;
				DRMNFB(&b[0], &d[0], &fx, &iv[0], &liv, &lv, &n, &v[0], x);
				liv = iv[43]; // IVLAST
				iv.resize(liv);
				lv = iv[44]; // VLAST
				v.resize(lv);

				iv[0] = 0;
				DRMNFB(&b[0], &d[0], &fx, &iv[0], &liv, &lv, &n, &v[0], x);

				while (iv[0] <= 2) {
					int nf = iv[5]; // NFCALL
					f(&n, x, &nf, &fx, ui, ur, uf);
					DRMNFB(&b[0], &d[0], &fx, &iv[0], &liv, &lv, &n, &v[0], x);
					if (nf <= 0) {
						iv[1] = 1; // TOOBIG
					}
				}
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
