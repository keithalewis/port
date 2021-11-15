// port_allocation.h - portfolio allocation
#pragma once
#include <vector>
#include "port_dmn.h"

namespace port {


	class allocation {
	public:
		int n;
		const double* ER;
		const double* Cov; // !!!eliminate this
		std::vector<double> L; // packed cholesky LL' = Cov 
		std::vector<double> V_1, V_ER; // V^-1 1, V^-1 ER
		double A, B, C, D;
	public:
		allocation(int n, const double* ER, const double* Cov, bool upper = false)
			: n(n), ER(ER), Cov(Cov), L((n* (n + 1)) / 2), V_1(n, 1), V_ER(ER, ER + n)
		{
			// packed triangular symmetric
			if (upper) {
				packu(n, Cov, L.data());
			}
			else {
				packl(n, Cov, L.data());
			}

			// Cov = L L'
			int irc = cholesky(n, L.data(), L.data());
			// Cholesky decomposition L.data(), &irc);
			if (0 != irc) {
				throw irc;
			}

			// V^-1 x = L'^-1 L^-1 x
			DL7IVM(&n, V_1.data(), L.data(), V_1.data());
			DL7ITV(&n, V_1.data(), L.data(), V_1.data());

			// V^-1 ER = L'^-1 L^-1 ER
			DL7IVM(&n, V_ER.data(), L.data(), V_ER.data());
			DL7ITV(&n, V_ER.data(), L.data(), V_ER.data());

			double one = 1;
			A = dot(n, V_1.data(), &one, 0);
			B = dot(n, V_1.data(), ER);
			C = dot(n, V_ER.data(), ER);
			D = A * C - B * B;
		}
		~allocation()
		{ }

		// Return minimum volatiltiy given expected realized return
		// minimize (1/2) xi' Cov xi - lambda (xi . 1 - 1) - mu (xi . ER - R)
		// D_xi F = Cov xi - lambda . 1 - mu ER, D^2_xi F = Cov;
		// 
		// 0 = V xi - lambda 1 - mu ER
		// xi = V^-1(lambda 1 + mu ER) = lambda V^-1 1 + mu V^-1 ER
		// 
		// 1 = 1' xi  = lambda 1' V^-1 1 + mu 1' V^-1 ER   = lambda A + mu B
		// R = ER' xi = lambda ER' V^-1 1 + mu ER' V^-1 ER = lambda B + mu C
		// 
		// [ 1 ] = [ A B ] [ lambda ]
		// [ R ]   [ B C ] [ mu ]
		// 
		// [ lambda ] = 1/(AC - B^2) [ C -B ] [ 1 ] = [ (C - RB)/D  ]
		// [ mu     ]                [ -B A ] [ R ]   [ (-B + RA)/D ]
		//
		// (xi' V) xi = (lambda 1' + mu ER') (lambda V^-1 1 + mu V^-1 ER)
		//            = lambda^2 A + 2 lambda mu B + mu^2 C

		double minimum(double R, double* x) const
		{
			double lambda = (C - R * B) / D;
			double mu = (-B + R * A) / D;

			// xi = lambda V^-1 x + mu V^-1 EX
			for (int i = 0; i < n; ++i) {
				x[i] = lambda * V_1[i] + mu * V_ER[i];
			}
			x[n] = lambda;
			x[n + 1] = mu;

			double sigma = (A * R * R - 2 * B * R + C) / D;
			sigma = sqrt(sigma);

			return sigma;
		}

		double minimize(double R, double* x, const double* l = 0, const double* u = 0) const
		{
			if (!l and !u) { // closed form
				return minimum(R, x);
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
				double R = *UR;
				// F(xi, lambda, mu) = (xi' Cov xi)/2 - lambda (xi'1 - 1) - mu(xi' ER - R) 
				*F = quad(n, X, p->L.data()) / 2;
				double one = 1;
				*F -= lambda * (dot(n, X, &one, 0) - 1);
				*F -= mu * (dot(n, X, p->ER) - R);
				*F = *F;
			};

			auto g = [](int* /*N*/, double* X, int* /*NF*/, double* G, int* /*UI*/, double* UR, void* UF) {
				allocation* p = (allocation*)UF;
				int n = p->n;
				double lambda = X[n];
				double mu = X[n + 1];
				double R = *UR;

				// G = D_xi F = Cov xi - lambda . 1 - mu ER
				DL7TVM(&n, G, p->L.data(), X); // G = L' xi
				DL7VML(&n, G, p->L.data(), G); // G = L G = L L' xi
				for (int i = 0; i < n; ++i) {
					G[i] += -lambda - mu * p->ER[i];
				}

				double one = 1;
				G[n] = -(dot(n, X, &one, 0) - 1);
				G[n + 1] = -(dot(n, X, p->ER) - R);
			};

			auto h = [](int* /*N*/, double* /*X*/, int* /*NF*/, double* H, int* /*UI*/, double* /*UR*/, void* UF) {
				allocation* p = (allocation*)UF;
				int n = p->n;
				DL7SQR(&n, H, p->L.data());
			};

			RETURN_CODE ret;
			ret = p.solve(x, f, 0, &R, (void*)this);

			return 0; // p.v[?];
		}

		// Maximize expected return given volatility
		// maximiz xi' EX - lambda (xi' x - 1) - mu/2 (xi' V xi - sigma^2)
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
			double a = A - sigma * sigma * A * A;
			double b = B - sigma * sigma * A * B;
			double c = C - sigma * sigma * B * B;
			double d = sqrt(b * b - a * c);

			double lambda = (b + d) / a; // +- d ???
			double mu = B - lambda * A;

			// xi = (V_EX - lambda V_X)/mu
			for (int i = 0; i < n; ++i) {
				x[i] = (V_ER[i] - lambda * V_1[i]) / mu;
			}
			x[n] = lambda;
			x[n + 1] = mu;

			return (C - lambda * B) / mu;
		}
		double maximize(double sigma, double* x, const double* l = 0, const double* u = 0) const
		{
			if (!l && !u) { // closed form
				return 	maximum(sigma, x);
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
				*F += mu * (quad(n, X, p->L.data()) - sigma * sigma) / 2;
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
				G[n + 1] = (quad(n, X, p->L.data()) - sigma * sigma) / 2;
			};

			RETURN_CODE ret;
			ret = p.solve(x, f, g, 0, &sigma, (void*)this);

			return 0; // p.v[?];
		}
	};

} // namespace port
