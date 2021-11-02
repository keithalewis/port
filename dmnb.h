// port.h - PORT library header
#pragma once
#include "cport.h"
#include <vector>

namespace port {

	// minimize constrained objective function
	class dmnfb {
		int n, liv, lv;
		std::vector<int> iv;
		std::vector<double> d, v, x, b;
	public:
		typedef void(*QF)(int* N, double* X, int* NF, double* F, int* UI, double* UR, void* UF);

		dmnfb(int n, const double* x, const double* l, const double* u)
			: n(n), liv(59 + n), lv(77 + n * (n + 23) / 2), iv(59 + n), d(n, 1), v(lv), x(x, x + n), b(2*n)
		{
			iv[0] = 0;
			// convert bounds to column major format
			for (int i = 0; i < n; ++i) {
				b[2 * i] = l[i];
				b[2 * i + 1] = u[i];
			}
		}
		dmnfb(const dmnfb&) = delete;
		dmnfb& operator=(const dmnfb&) = delete;
		~dmnfb()
		{ }
		void solve(QF qf, int* ui = 0, double* ur = 0, void* dummy = 0)
		{
			DMNFB(&n, d.data(), x.data(), b.data(), qf, iv.data(), &liv, &lv, v.data(), ui, ur, dummy);
		}
	};

	// class dmng // gradient
	// class dmnh // hession

} // namespace port
