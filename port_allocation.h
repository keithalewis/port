// port_allocation.h - portfolio allocation
#pragma once
#include <vector>
#include "port_dmn.h"

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

	// x' LL' x = ||L'x||^2 , L packed
	inline double quad(int n, double* x, double* L, double* work)
	{
		DL7TVM(&n, x, L, work);
		
		return dot(n, work, work);
	}

	// Cholesky decomposition of packed A into packed L. A = LL'
	inline int cholesky(int n, double* A, double* L)
	{
		int n1 = 1;
		int irc;
		// Cholesky decomposition
		DL7SRT(&n1, &n, L, A, &irc);

		return irc;
	}

	class allocation {
public:
		int n;
		const double* ER;
		const double* Cov; // !!!eliminate this
		std::vector<double> L; // packed cholesky LL' = Cov 
		std::vector<double> V_x, V_ER; // V^-1 x, V^-1 ER
		double A, B, C, D;
	public:
		allocation(int n, const double* ER, const double* Cov, bool upper = false)
			: n(n), ER(ER), Cov(Cov), L((n*(n + 1))/2), V_x(n, 1), V_ER(ER, ER + n)
		{
			// packed triangular symmetric
			if (upper) {
				packu(n, Cov, L.data());
			}
			else {
				packl(n, Cov, L.data());
			}
		
			int irc = cholesky(n, L.data(), L.data());
			// Cholesky decomposition L.data(), &irc);
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

		// return minimum volatiltiy given expected realized return
		double minimum(double R, double* x) const
		{
			double lambda = (C - R * B) / D;
			double mu = (-B + R * A) / D;

			// xi = lambda V^-1 x + mu V^-1 EX
			for (int i = 0; i < n; ++i) {
				x[i] = lambda * V_x[i] + mu * V_ER[i];
			}
			x[n] = lambda;
			x[n + 1] = mu;

			double sigma = sqrt((A * R * R - 2 * B * R + C) / D);

			return sigma;
		}

		// minimize (1/2) xi' Cov xi - lambda (xi . 1 - 1) - mu (xi . ER - R)
		// D_xi F = Cov xi - lambda . 1 - mu ER
		// D^2_xi F = Cov;
		double minimize(double R, double* x, const double* l = 0, const double* u = 0) const
		{
			if (!l and !u) { // closed form
				return minimum(R, x);
			}

			dmn p(n + 2);
			p.lower(l);
			p.upper(u);

			auto f = [](int* /*N*/, double* X, int* /*NF*/, double* F, int* /*UI*/, double* UR, void* UF) {
				allocation* p = (allocation*)UF;
				int n = p->n;
				double lambda = X[n];
				double mu = X[n + 1];
				double R = *UR;
				// F(xi, lambda, mu) = (xi' Cov xi)/2 - lambda (xi'1 - 1) - mu(xi' ER - R) 
				*F = quad(n, X, p->L.data(), X) / 2;
				double one = 1;
				*F -= lambda * (dot(n, X, &one, 0) - 1);
				*F -= mu * (dot(n, X, p->ER) - R);
			};

			auto g = [](int* /*N*/, double* X, int* /*NF*/, double* G, int* /*UI*/, double* UR, void* UF) {
				allocation* p = (allocation*)UF;
				int n = p->n;
				double lambda = X[n];
				double mu = X[n + 1];
				double R = *UR;

				// G = D_xi F = Cov xi - lambda . 1 - mu ER
				for (int i = 0; i < n; ++i) {
					G[i] = dot(n, X, p->Cov + i * n) - lambda - mu*p->ER[i];
				}
				double one = 1;
				G[n] = -(dot(n, X, &one, 0) - 1);
				G[n + 1] = -(dot(n, X, p->ER) - R);
			};

			//auto h = [](int* N, double* X, int* /*NF*/, double* H, int* /*UI*/, double* UR, void* UF) {
		
			RETURN_CODE ret;
			ret = p.solve(x, f, g, 0, &R, (void*)this);

			return 0; // p.v[?];
		}

		// x = (xi, lambda, mu)
		// max xi' EX - lambda (xi' x - 1) - mu/2 (xi' V xi - sigma^2)
		// 0 = EX - lambda x - mu V xi
		// xi = V^{-1}(EX - lambda x)/mu
		// 1 = xi' x = (B - lambda A)/mu so mu = B - lambda A
		// sigma^2 = xi' V xi 
		//         = (C - 2B^2 lambda + A lambda^2)/mu^2
		//         = (C - 2B^2 lambda + A lambda^2)/(B - lambda A)^2 
		// 0 = (C - 2B lambda + A lambda^2) - sigma^2(B^2 - 2 AB lambda + A^2 lambda^2)
		// 0 = (C - sigma^2 B^2) - 2(B - sigma^2 AB) lambda + (A - sigma^2 A^2) lambda^2
		// return xi' EX = (C - lambda B)/mu
		double maximum(double sigma, double* x) const
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

			return (C - lambda * B) / mu;
		}
		double maximize(double sigma, double* x, const double* l = 0, const double* u = 0) const
		{
			double R = maximum(sigma, x);

			if (!l && !u) { // closed form
				return R;
			}

			dmn p(n + 2);

			if (l) {
				p.lower(l);
			}
			if (u) {
				p.upper(u);
			}

			auto f = [](int* /*N*/, double* X, int* /*NF*/, double* F, int* /*UI*/, double* UR, void* UF) {
				allocation* p = (allocation*)UF;
				int n = p->n;
				double lambda = X[n];
				double mu = X[n + 1];
				double sigma = *UR;
				// F = xi' ER - lambda (xi'1 - 1) - mu (xi' Cov xi - sigma^2)/2
				*F = -dot(n, X, p->ER);
				double one = 1;
				*F += lambda * (dot(n, X, &one, 0) - 1);
				*F += mu * (quad(n, X, p->L.data(), X) - sigma * sigma) / 2;
			};

			auto g = [](int* /*N*/, double* X, int* /*NF*/, double* G, int* /*UI*/, double* UR, void* UF) {
				allocation* p = (allocation*)UF;
				int n = p->n;
				double lambda = X[n];
				double mu = X[n + 1];
				double sigma = *UR;

				// G = D_xi F = ER - lambda . 1 - mu Cov xi
				for (int i = 0; i < n; ++i) {
					G[i] = -p->ER[i];
					G[i] += lambda;
					G[i] += mu * dot(n, X, p->Cov + i * n);
				}
				double one = 1;
				G[n] = (dot(n, X, &one, 0) - 1);
				G[n + 1] = (quad(n, X, p->L.data(), X) - sigma * sigma) / 2;
			};

			RETURN_CODE ret;
			ret = p.solve(x, f, g, 0, &sigma, (void*)this);

			return 0; // p.v[?];
		}
	};

} // namespace port

