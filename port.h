// port.h - PORT library header
#pragma once
#include "cport.h"
#include <vector>

namespace port {

	// minimize unconstrained objective function
	class dmnf {
		int liv, lv, n;
		int iv[59];
		std::vector<double> d, v, x;
	public:
		typedef void(*QF)(int* N, double* X, int* NF, double* F, int* UI, double* UR, void* UF);

		dmnf(int n, const double* x)
			: liv(59), lv(77 + n * (n + 17) / 2), n(n), d(n, 1), v(lv), x(x, x + n)
		{
			//int kind = 2;
			iv[0] = 0;
			//DIVSET(&kind, iv, &liv, &lv, v.data());
		}
		dmnf(const dmnf&) = delete;
		dmnf& operator=(const dmnf&) = delete;
		~dmnf()
		{ }
		void solve(QF qf)
		{
			DMNF(&n, d.data(), x.data(), qf, iv, &liv, &lv, v.data(), 0/*UI*/, 0/*UR*/, 0/*DUMMY*/);
		}
	};

} // namespace port
