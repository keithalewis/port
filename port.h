// port.h - PORT library header
#pragma once
#include <vector>

extern "C" {
	void DIVSET(int* ALG, int* IV, int* LIV, int* LV, double* V);
	void DRMNF(double* D, double* FX, int* IV, int* LIV, int* LV, int* N, double* V, double* X);
}

namespace port {

	// minimize unconstrained objective function
	class drmn {
		int liv, lv, n;
		int iv[45];
		std::vector<double> d, v, x;
	public:
		drmn(int n, const double* x)
			: liv(45), lv(77 + n * (n + 17) / 2), n(n), d(n, 1), v(lv), x(x, x + n)
		{
			int kind = 2;
			DIVSET(&kind, iv, &liv, &lv, v.data());
		}
	};

} // namespace port
