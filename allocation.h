// allocation.h - portfolio allocation
#pragma once
#include <vector>
#include "dmn.h"

namespace port {

	// x . y
	inline double dot(size_t n, const double* x, const double* y, size_t stride = 1)
	{
		double s = 0;

		for (size_t i = 0; i < n; ++i) {
			s += x[i] * y[i*stride];
		}

		return s;
	}

	// x' A x
	inline double quad(size_t n, const double* x, const double* A)
	{
		double q = 0;

		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j) {
				q += A[n * i + j] * x[i] * x[j];
			}
		}

		return q;
	}

	class allocation {
		int n;
		const double* ER;
		const double* Cov;
		std::vector<double> L; // LL' = Cov packed cholesky
		std::vector<double> V_x, V_ER; // V^-1 x, V^-1 ER
		double A, B, C, D;
	public:
		allocation(int n, const double* ER, const double* Cov, bool upper = false)
			: n(n), ER(ER), Cov(Cov), L((n*(n + 1))/2), V_x(n, 1), V_ER(ER, ER + n)
		{
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j <= i; ++j) {
					L[(i * (i + 1)) / 2 + j] = Cov[upper ? j + n*i : n*i + j];
				}
			}
			int n1 = 1;
			int irc;
			DL7SRT(&n1, &n, L.data(), L.data(), &irc);
			if (0 != irc) {
				throw irc;
			}

			// V^-1 x
			DL7ITV(&n, V_x.data(), L.data(), V_x.data());
			DL7IVM(&n, V_x.data(), L.data(), V_x.data());

			// V^-1 ER
			DL7ITV(&n, V_ER.data(), L.data(), V_ER.data());
			DL7IVM(&n, V_ER.data(), L.data(), V_ER.data());

			double one = 1;
			A = dot(n, V_x.data(), &one, 0);
			B = dot(n, V_x.data(), ER);
			C = dot(n, V_ER.data(), ER);
			D = A * C - B * B;
		}
		~allocation()
		{ }

		// minimize (1/2) xi' Cov xi given xi . 1 = 1 and xi . ER = R
		// D_xi F = Cov xi - lambda . 1 - mu ER
		// D^2_xi F = Cov;
		double minimize(double R, double* x, const double* l = 0, const double* u = 0) const
		{
			double lambda = (C - R * B) / D;
			double mu = (-B + R * A) / D;

			// xi = lambda V^-1 x + mu V^-1 EX
			for (int i = 0; i < n; ++i) {
				x[i] = lambda * V_x[i] + mu * V_ER[i];
			}
			x[n] = lambda;
			x[n + 1] = mu;

			if (!l and !u) { // closed form
				return sqrt((A * R * R - 2 * B * R + C) / D);
			}

			dmn p(n + 2);

			if (l) {
				p.lower(l);
				for (int i = 0; i < n; ++i) {
					if (x[i] < l[i]) {
						x[i] = l[i];
					}
				}
			}
			if (u) {
				p.upper(u);
				for (int i = 0; i < n; ++i) {
					if (x[i] > u[i]) {
						x[i] = u[i];
					}
				}
			}

			auto f = [](int* N, double* X, int* /*NF*/, double* F, int* /*UI*/, double* UR, void* UF) {
				int n = *N - 2;
				double lambda = X[n];
				double mu = X[n + 1];
				double R = *UR;
				allocation* p = (allocation*)UF;
				// F = (xi' Cov xi)/2 - lambda (xi'1 - 1) - mu(xi' ER - R) 
				*F = quad(n, X, p->Cov) / 2;
				double one = 1;
				*F -= lambda * (dot(n, X, &one, 0) - 1);
				*F -= mu * (dot(n, X, p->ER) - R);
			};

			auto g = [](int* N, double* X, int* /*NF*/, double* G, int* /*UI*/, double* UR, void* UF) {
				int n = *N - 2;
				double lambda = X[n];
				double mu = X[n + 1];
				double R = *UR;
				allocation* p = (allocation*)UF;

				// G = D_xi F = Cov xi - lambda . 1 - mu ER
				for (int i = 0; i < n; ++i) {
					G[i] = dot(n, X, p->Cov + i * n) - lambda - mu*p->ER[i];
				}
				double one = 1;
				G[n] = -(dot(n, X, &one, 0) - 1);
				G[n + 1] = -(dot(n, X, p->ER) - R);
			};

			RETURN_CODE ret;
			ret = p.solve(x, f, g, 0, &R, (void*)this);

			return 0; // p.v[?];
		}

		// maximize xi' ER given xi . 1 = 1 and xi' Cov xi = sigma^2
		// D_xi F = ER - lambda . 1 - mu Cov xi
		// D^2_xi F = - mu Cov
		// x = (xi, lambda, mu)
		double maximize(double sigma, double* x, const double* l = 0, const double* u = 0) const
		{
			double c = C - sigma * sigma * B * B;
			double b = B - sigma * sigma * A * B;
			double a = A - sigma * sigma * A * A;
			double d = sqrt(b * b - a * c);

			double lambda = (b + d) / a; // +- d ???
			double mu = B - lambda * A;

			// xi = (V_EX - lambda V_X)/mu
			for (int i = 0; i < n; ++i) {
				x[i] = (V_ER[i] - lambda * V_x[i]) / mu;
			}
			x[n] = lambda;
			x[n + 1] = mu;

			if (!l && !u) { // closed form
				return (C - lambda * B) / mu;
			}

			dmn p(n + 2);

			if (l) {
				p.lower(l);
				for (int i = 0; i < n; ++i) {
					if (x[i] < l[i]) {
						x[i] = l[i];
					}
				}
			}
			if (u) {
				p.upper(u);
				for (int i = 0; i < n; ++i) {
					if (x[i] > u[i]) {
						x[i] = u[i];
					}
				}
			}

			auto f = [](int* N, double* X, int* /*NF*/, double* F, int* /*UI*/, double* UR, void* UF) {
				int n = *N - 2;
				double lambda = X[n];
				double mu = X[n + 1];
				double sigma = *UR;
				allocation* p = (allocation*)UF;
				// F = xi' ER - lambda (xi'1 - 1) - mu (xi' Cov xi - sigma^2)/2
				*F = dot(n, X, p->ER);
				double one = 1;
				*F -= lambda * (dot(n, X, &one, 0) - 1);
				*F -= mu * (quad(n, X, p->Cov) - sigma * sigma) / 2;
			};

			auto g = [](int* N, double* X, int* /*NF*/, double* G, int* /*UI*/, double* UR, void* UF) {
				int n = *N - 2;
				double lambda = X[n];
				double mu = X[n + 1];
				double sigma = *UR;
				allocation* p = (allocation*)UF;

				// G = D_xi F = ER - lambda . 1 - mu Cov xi
				std::copy(p->ER, p->ER + n, G);
				for (int i = 0; i < n; ++i) {
					G[i] -= lambda;
					G[i] -= mu * dot(n, X, p->Cov + i * n);
				}
				double one = 1;
				G[n] = -(dot(n, X, &one, 0) - 1);
				G[n + 1] -= (quad(n, X, p->Cov) - sigma * sigma) / 2;
			};

			RETURN_CODE ret;
			ret = p.solve(x, f, g, 0, &sigma, (void*)this);

			return 0; // p.v[?];
		}
	};

} // namespace port

