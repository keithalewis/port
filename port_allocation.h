// port_allocation.h - portfolio allocation
#pragma once
#include <valarray>
#include <vector>
#include "port.h"

namespace port {

	class allocation {
	public:
		int N;
		std::vector<double> ER; // expected realized returns
		std::vector<double> L; // packed cholesky LL' = Cov 
		std::vector<double> V_1, V_ER; // V^-1 1, V^-1 ER
		double A, B, C, D;

		// Assumes V_1 = 1 and V_ER = ER on entry
		void init()
		{
			// V^-1 1 = L'^-1 L^-1 1
			DL7IVM(&N, &V_1[0], &L[0], &V_1[0]);
			DL7ITV(&N, &V_1[0], &L[0], &V_1[0]);

			// V^-1 ER = L'^-1 L^-1 ER
			DL7IVM(&N, &V_ER[0], &L[0], &V_ER[0]);
			DL7ITV(&N, &V_ER[0], &L[0], &V_ER[0]);

			double one = 1;
			A = dot(N, &V_1[0], &one, 0);
			B = dot(N, &V_1[0], &ER[0]);
			C = dot(N, &V_ER[0], &ER[0]);
			D = A * C - B * B;
		}

		// xi' 1
		double xi_1(const double* xi)
		{
			double one = 1;
			return dot(N, xi, &one, 0);
		}
		// xi' ER
		double xi_ER(const double* xi)
		{
			return dot(N, xi, &ER[0]);
		}
		// xi' Cov xi
		double Var(const double* xi)
		{
			return quad(N, xi, &L[0]);
		}
		// F(xi, lambda, mu) = (xi' Cov xi)/2 - lambda (xi'1 - 1) - mu(xi' ER - R) 
		double fmin(double R, const double* xi, double lambda, double mu) 
		{
			return Var(xi) / 2 - lambda * (xi_1(xi) - 1) - mu * (xi_ER(xi) - R);
		}
		// xi_ = a xi + b 1, 1 = xi_' 1, R = xi_' ER
		static void normalize(std::valarray<double>& xi, double R, const double* ER)
		{
			auto N = xi.size();
			double _1xi = xi.sum(); // 1' xi
			double ER_xi = dot(N, ER, &xi[0]); // ER' xi
			double one = 1;
			double ER_1 = dot(N, ER, &one, 0); // ER' 1
			double d_ = _1xi * ER_1 - ER_xi * N;
			double a_ = (ER_1 - N * R) / d_;
			double b_ = (-ER_xi + _1xi * R) / d_;
			
			xi *= a_;
			xi += b_;
		}
		// G = DF(xi, lambda, mu) = (Cov xi - lambda 1 - mu ER, -xi' 1  + 1, -xi' ER + R)
		void gmin(double R, const double* xi, double lambda, double mu, double* g)
		{
			// g = L' xi
			DL7TVM(&N, g, &L[0], xi);
			// g = L g
			DL7VML(&N, g, &L[0], g);
			
			for (int i = 0; i < N; ++i) {
				g[i] -= lambda;
				g[i] -= mu * ER[i];
			}

			double one = 1;
			g[N] = -dot(N, xi, &one, 0) + 1;
			g[N + 1] = -dot(N, xi, &ER[0]) + R;
		}
		// F(xi, lambda, mu) = xi' ER - lambda (xi'1 - 1) - mu/2 (xi' Cov xi - sigma^2)
		double fmax(double sigma, double* xi, double lambda, double mu) 
		{
			return xi_ER(xi) - lambda * (xi_1(xi) - 1) - mu * (Var(xi) - sigma*sigma)/2;
		}
	public:
		/*
		// L is lower cholesky factor
		allocation(int N, const double* ER_, const double* L_)
			: dmn(N + 2), N(N), ER(ER_, ER_ + N), L(L_, L_ + (N* (N + 1)) / 2), V_1(N, 1), V_ER(ER_, ER_ + N)
		{
			init();
		}
		*/

		allocation(int N, const double* ER_, const double* Cov, bool upper)
			: N(N), ER(ER_, ER_ + N), L((N* (N + 1)) / 2), V_1(N, 1), V_ER(ER_, ER_ + N)
		{
			// packed triangular symmetric
			if (upper) {
				packu(N, Cov, &L[0]);
			}
			else {
				packl(N, Cov, &L[0]);
			}

			// Cov = L L'
			int irc = cholesky(N, &L[0], &L[0]);
			if (0 != irc) {
				throw irc;
			}

			init();
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
			if (x) {
				for (int i = 0; i < N; ++i) {
					x[i] = lambda * V_1[i] + mu * V_ER[i];
				}
				x[N] = lambda;
				x[N + 1] = mu;
			}

			double sigma = (A * R * R - 2 * B * R + C) / D;
			sigma = sqrt(sigma);

			return sigma;
		}
		// F(x + h) = F(x) + DF(x)h + o(||h||)
		// x' = x - gamma DF(x) >= 0
		double stepmin(double R, double* x)
		{
			double lambda = x[N];
			double mu = x[N + 1];

			std::valarray<double> g(N + 2);
			gmin(R, x, lambda, mu, &g[0]);

			// largest gamma for which xi >= 0
			double gamma = DBL_MAX;
			for (int i = 0; i < N; ++i) {
				gamma = std::min(gamma, x[i] / g[i]);
			}
			for (int i = 0; i < N + 2; ++i) {
				x[i] -= gamma * g[i];
			}
			lambda = x[N];
			mu = x[N + 1];

			return gamma;
		}
		double minimize(double R, double* x)
		{
			double lambda = x[N];
			double mu = x[N + 1];
			double f0 = fmin(R, x, lambda, mu);

			std::valarray<double> g(N + 2);
			gmin(R, x, lambda, mu, &g[0]);

			// largest gamma for which xi >= 0
			double gamma = DBL_MAX;
			for (int i = 0; i < N; ++i) {
				gamma = std::min(gamma, x[i]/g[i]);
			}
			for (int i = 0; i < N + 2; ++i) {
				x[i] -= gamma * g[i];
			}
			lambda = x[N];
			mu = x[N + 1];
			double f = fmin(R, x, lambda, mu);

			return f - f0;
		}
#if 0
		double minimize(double R, double* x)
		{
			/*
			if (b.size() == 0) { // closed form
				return minimum(R, x);
			}
			*/
			auto f = [](int* /*N*/, double* X, int* /*NF*/, double* F, int* /*UI*/, double* UR, void* UF) {
				allocation* p = (allocation*)UF;

				int N = p->N;
				double lambda = X[N];
				double mu = X[N + 1];
				double R = *UR;
				
				*F = p->fmin(R, X, lambda, mu);
			};

			auto g = [](int* /*N*/, double* X, int* /*NF*/, double* G, int* /*UI*/, double* UR, void* UF) {
				allocation* p = (allocation*)UF;

				int N = p->N;
				double lambda = X[N];
				double mu = X[N + 1];
				double R = *UR;

				p->gmin(R, X, lambda, mu, G);
			};

			auto h = [](int* /*N*/, double* /*X*/, int* /*NF*/, double* H, int* /*UI*/, double* /*UR*/, void* UF) {
				allocation* p = (allocation*)UF;
				int N = p->N;
				DL7SQR(&N, H, &p->L[0]); // !!! not correct, need lambda, mu hessian
			};

			RETURN_CODE ret;
			ret = dmn::solve(x, f, 0, &R, (void*)this);
			//ret = dmn::solve(x, f, g, 0, &R, (void*)this);

			return v[10];
		}
		// Maximize expected return given volatility
		// maximize xi' EX - lambda (xi' x - 1) - mu/2 (xi' V xi - sigma^2)
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
			double a_ = A - sigma * sigma * A * A;
			double b_ = B - sigma * sigma * A * B;
			double c_ = C - sigma * sigma * B * B;
			double d_ = sqrt(b_ * b_ - a_ * c_);

			double lambda = (b_ + d_) / a_; // +- d ???
			double mu = B - lambda * A;

			// xi = (V_EX - lambda V_X)/mu
			if (x) {
				for (int i = 0; i < N; ++i) {
					x[i] = (V_ER[i] - lambda * V_1[i]) / mu;
				}
				x[N] = lambda;
				x[N + 1] = mu;
			}

			return (C - lambda * B) / mu;
		}
		double maximize(double sigma, double* x, const double* l = 0, const double* u = 0) const
		{
			if (!l && !u) { // closed form
				return 	maximum(sigma, x);
			}

			dmn p(N + 2);

			if (l) {
				p.lower(l);
			}
			if (u) {
				p.upper(u);
			}

			auto f = [](int* /*N*/, double* X, int* /*NF*/, double* F, int* /*UI*/, double* UR, void* UF) {
				allocation* p = (allocation*)UF;
				int N = p->N;
				double lambda = X[N];
				double mu = X[N + 1];
				double sigma = *UR;
				
				*F = p->fmax(sigma, X, lambda, mu);
			};

			auto g = [](int* /*N*/, double* X, int* /*NF*/, double* G, int* /*UI*/, double* UR, void* UF) {
				allocation* p = (allocation*)UF;
				int N = p->N;
				double lambda = X[N];
				double mu = X[N + 1]; mu = mu;
				double sigma = *UR;

				// G = D_xi F = ER - lambda . 1 - mu Cov xi
				for (int i = 0; i < N; ++i) {
					G[i] = -p->ER[i];
					G[i] += lambda;
					//!!!G[i] += mu * dot(N, X, p->Cov + i * N);
				}
				double one = 1;
				G[N] = (dot(N, X, &one, 0) - 1);
				G[N + 1] = (quad(N, X, &p->L[0]) - sigma * sigma) / 2;
			};

			RETURN_CODE ret;
			ret = p.solve(x, f, 0, &sigma, (void*)this);

			return -p.v[10]; // optimal value
		}
#endif 0
	};

} // namespace port
