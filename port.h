// cport.h - C declarations of PORT functions
#pragma once
#include <vector>

extern "C" {
	double D1MACH(int*);
	int I1MACH(int*);
	void DIVSET(int* ALG, int* IV, int* LIV, int* LV, double* V);

	// W = A X + Y
	void DV2AXY(int* N, double* W, const double* A, const double* X, const double* Y);
	
	// compute X = LY
	void DL7VML(const int* N, double* X, const double* L, const double* Y);
	// compute X = L'Y
	void DL7TVM(const int* N, double* X, const double* L, const double* Y);
	// solve LX = Y
	void DL7IVM(int* N, double* X, double* L, double* Y);
	// solve L'X = Y
	void DL7ITV(int* N, double* X, double* L, double* Y);
	
	// Cholesky factor A = LL'
	void DL7SRT(const int* N1, const int* N, double* L, const double* A, int* IRC);
	// Set A to lower triangle of LL'
	void DL7SQR(int* N, double* A, const double* L);
	// Set A to lower triangle of L'L
	void DL7TSQ(int* N, double* A, const double* L);

	// Compute finite difference gradient by Stewart*s scheme
	void DS7GRD(double* ALPHA, double* D, double* ETA0, double* FX, double* G, int* IRC, int* N, double* W, double* X);
	
	// minimization routines
	void DMNF(int* N, double* D, double* X, void* F, int* IV, int* LIV, int* LV, double* V, int* UI, double* UR, void* DUMMY);
	void DMNG(int* N, double* D, double* X, void* F, void* G, int* IV, int* LIV, int* LV, double* V, int* UI, double* UR, void* DUMMY);
	void DMNH(int* N, double* D, double* X, void* F, void* GH, int* IV, int* LIV, int* LV, double* V, int* UI, double* UR, void* DUMMY);
	void DMNFB(int* N, double* D, double* X, double* B, void* F, int* IV, int* LIV, int* LV, double* V, int* UIPARM, double* URPARM, void* UFPARM);
	void DMNGB(int* N, double* D, double* X, double* B, void* F, void* G, int* IV, int* LIV, int* LV, double* V, int* UIPARM, double* URPARM, void* UFPARM);
	void DMNHB(int* N, double* D, double* X, double* B, void* F, void* GH, int* IV, int* LIV, int* LV, double* V, int* UIPARM, double* URPARM, void* UFPARM);
	void DRMNF(double* D, double* FX, int* IV, int* LIV, int* LV, int* N, double* V, double* X);
	void DRMNFB(double* B, double* D, double* FX, int* IV, int* LIV, int* LV, int* N, double* V, double* X);
}

namespace port {

	enum IV {
		COVMAT = 26,
		COVPRT = 14,
		COVREQ = 15,
		D = 27,
		DRADPR = 101,
		DTOL = 59,
		DTYPE = 16,
		INITS = 25,
		LASTIV = 44,
		LASTV = 45,
		MC = 83,
		ME = 86,
		ME1 = 87,
		MODE = 35,
		MODEL = 5,
		MXFCAL = 17,
		NEXTIV = 46,
		NEXTV = 47,
		NFCALL = 6,
		NFCOV = 52,
		NFGCAL = 7,
		NGCALL = 30,
		NGCOV = 53,
		NITER = 31,
		OUTLEV = 19,
		PARPRT = 20,
		PC = 90,
		PRUNIT = 21,
	};

	enum V {
		AFCTOL = 31,
		COSMIN = 47,
		D0INIT = 40,
		DELTA0 = 44,
		DFAC = 41,
		DGNORM = 1,
		DINIT = 38,
		DLTFDC = 42,
		DLTFDJ = 43,
		DSTNRM = 2,
		DTINIT = 39,
		ETA0 = 42,
		F = 10,
		F0 = 13,
		FUZZ = 45,
		LMAX0 = 35,
		LMAXS = 36,
		NREDUC = 6,
		PREDUC = 7,
		RADIUS = 8,
		RCOND = 53,
		RELDX = 17,
		RFCTOL = 32,
		RLIMIT = 46,
		SCTOL = 37,
		STPPAR = 5,
		XCTOL = 33,
		XFTOL = 34,
	};

	// First value of IV after regression or optimization is run
	enum class RETURN_CODE {
		X_CONV = 3,    // — X - convergence: the current iterate appears to be a scaled distance(see §5) of at most
		               // V(XCTOL) = V(33) from a locally optimal point.
		REL_CONV = 4,  // — relative function convergence : the current objective function value f(x) appears to differ from
		               // a locally optimal value by at most |f(x)| V(RFCTOL) 
		BOTH_CONV = 5, // — both X - and relative function convergence(3 and 4 combined).
		ABS_CONV = 6,  // — absolute function convergence : |f(x)| < V(AFCTOL) = V(31).
		               // This test is only of interest in problems where f(x) = 0 means 
					   // a "perfect fit", such as nonlinear least - squares problems.
		//
		// Error returns from which restarts(§13) are possible :
		//
		SING_CONV = 7,       // — singular convergence : x may have too many free components.See §5.
		FALSE_CONV = 8,      // — false convergence : the gradient ∇ f(x) may be computed incorrectly, the other stopping tolerances may be too tight, or either f or ∇ f may be discontinuous near the current iterate x.
		FUNC_EVAL_LIMIT = 9, // — function evaluation limit : no convergence after IV(MXFCAL) = IV(17) evaluations of f(x).
		ITER_LIM = 10,       // — iteration limit : no convergence after IV(MXITER) = IV(18) iterations.
		STOPX = 11,          // — STOPX returned.TRUE. : you supplied a system - dependent STOPX(see §12) routine and hit the BREAK key.
		// 12–13 — impossible : these are input IV(1) values only. (12 means allocate storage within IV and V and
		// start the algorithm; this is the default IV(1) value supplied by[D]IVSET. 13 means just allocate storage and return.See §4a for an example.)
		STORAGE = 14,        // — storage has been allocated(after a call with IV(1) = 13 — see, for example, §4a below).
		//
		// Error returns that preclude restarts :
		//
		LIV_SMALL = 15, // — LIV too small.
		LV_SMALL = 16,  // — LV too small.
		RESTART = 17,   // — restart attempted(§13) with problem dimensions changed.
		NEG_SCALE = 18, // — d has a negative component and IV(DTYPE) ≤ 0 : see §4.
		//19–43 — V(IV(1)) is out of range.
		//44–62 — reserved.
		FUNC_COMPUTE = 63, // — f(x) cannot be computed at the initial x.
		BAD_PARAM = 64,    // — bad parameters on an internal call(should not occur).
		BAD_GRAD = 65,     // — the gradient could not be computed at x.
		BAD_INPUT = 66,    // — bad input array — if this return is relevant, see the associated PORT reference sheet
		BAD_FIRST = 67,    // — bad first parameter(KIND in §2) to IVSET.
		//68–69 — bugs encountered(should not occur).
		INIT_S = 70,       // — couldn’t get initial S matrix by finite differences(regression routines only, and only when
			               // IV(INITS) is at least 3).[The S matrix is an approximation to part of the Hessian matrix, ∇^2 f
		// 71–79 — reserved.
		BAD_IV1 = 80,         // — IV(1) was out of range(e.g.exceeded 14).
		BAD_DIM = 81,         // — bad problem dimensions(e.g.a nonpositive number of variables or , for regression routines, number of observations).
		BAD_BOUNDS = 82,      // — inconsistent bounds
		BAD_BOUNDS_GEN = 83,  // — inconsistent general bounds
		BAD_CONSTRAINT = 84,  // — some row of the constraint matrix, C, is all zeros.
		INC_CONSTRAINTS = 85, // — inconsistent constraints 
	};

	// pack lower triangle of a into l
	inline void packl(int n, const double* a, double* l)
	{
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j <= i; ++j) {
				l[(i * (i + 1)) / 2 + j] = a[n * i + j];
			}
		}
	}
	// unpack l into lower triangle of a
	inline void unpackl(int n, const double* l, double* a)
	{
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j <= i; ++j) {
				a[n * i + j] = l[(i * (i + 1)) / 2 + j];
			}
		}
	}
	// pack upper triangle of a into l
	inline void packu(int n, const double* a, double* l)
	{
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j <= i; ++j) {
				l[(i * (i + 1)) / 2 + j] = a[i + n * j];
			}
		}
	}
	// unpack l into upper triangle of a
	inline void unpacku(int n, const double* l, double* a)
	{
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j <= i; ++j) {
				a[i + n * j] = l[(i * (i + 1)) / 2 + j];
			}
		}
	}
	// unpack l into a
	inline void unpack(int n, const double* l, double* a)
	{
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j <= i; ++j) {
				a[j + n * i] = a[i + n * j] = l[(i * (i + 1)) / 2 + j];
			}
		}
	}

	// x . y
	inline double dot(size_t n, const double* x, const double* y, size_t stride = 1)
	{
		double s = 0;

		for (size_t i = 0; i < n; ++i) {
			s += x[i] * y[i * stride];
		}

		return s;
	}

	// x' LL' x = ||L'x||^2 , L packed
	inline double quad(int n, const double* x, const double* L)
	{
		std::vector<double> x_(x, x + n);

		DL7TVM(&n, x_.data(), L, x_.data());

		return dot(n, x_.data(), x_.data());
	}

	// Cholesky decomposition of packed A into packed L. A = LL'
	inline int cholesky(int n, const double* A, double* L)
	{
		int n1 = 1;
		int irc;
		// Cholesky decomposition
		DL7SRT(&n1, &n, L, A, &irc);

		return irc;
	}

}