// xll_port.cpp - Call PORT library
#include "port.h"
#include "xll/xll/xll.h"

#ifdef CATEGORY
#undef CATEGORY
#endif
#define CATEGORY "PORT"

using namespace xll;

AddIn xai_pack(
	Function(XLL_FPX, "xll_pack", "PACK")
	.Arguments({
		Arg(XLL_FPX, "A", "is a square matrix."),
		Arg(XLL_BOOL, "_upper", "is an optional argument indicating upper trangle of A is used. Default is lower.")
		})
	.FunctionHelp("Pack lower or upper triangle of A.")
	.Category(CATEGORY)
);
_FPX* WINAPI xll_pack(_FPX* pa, BOOL upper)
{
#pragma XLLEXPORT
	static FPX l;

	int n = pa->rows;
	if (n != pa->columns) {
		XLL_ERROR(__FUNCTION__ ": matrix must be square");

		return nullptr;
	}

	l.resize(1, (n * (n + 1)) / 2);

	if (upper) {
		port::packu(n, pa->array, l.array());
	}
	else {
		port::packl(n, pa->array, l.array());
	}

	return l.get();
}

AddIn xai_unpack(
	Function(XLL_FPX, "xll_unpack", "UNPACK")
	.Arguments({
		Arg(XLL_FPX, "L", "is a packed matrix."),
		})
	.FunctionHelp("Unpack L into lower or upper triangle of A.")
	.Category(CATEGORY)
);
_FPX* WINAPI xll_unpack(_FPX* pl)
{
#pragma XLLEXPORT
	static FPX a;
	int m = size(*pl);

	// m = n(n+1)/2
	// n^2 + n - 2m = 0
	// b^2 - 4ac = 1 + 8m
	auto d = sqrt(1 + 8 * m);
	int n = static_cast<int>((-1 + d) / 2);
	a.resize(n, n);

	port::unpack(n, pl->array, a.array());

	return a.get();
}

AddIn xai_d1mach(
	Function(XLL_DOUBLE, "xll_d1mach", "D1MACH")
	.Arguments({
		Arg(XLL_LONG, "index", "is the index."),
		})
	.FunctionHelp("Call PORT D1MACH.")
	.Category(CATEGORY)
	.Documentation(R"(
<pre>
C  DOUBLE-PRECISION MACHINE CONSTANTS
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C  D1MACH( 5) = LOG10(B)
</pre>
)")
);
double WINAPI xll_d1mach(int i)
{
#pragma XLLEXPORT
	return 1 <= i and i <= 5 ? D1MACH(&i) : XLL_NAN;
}

AddIn xai_dl7srt(
	Function(XLL_FPX, "xll_dl7srt", "DL7SRT")
	.Arguments({
		Arg(XLL_FPX, "A", "is a matrix."),
		Arg(XLL_LONG, "_n1", "is an optional row to start from. Default is 1.")
		})
	.FunctionHelp("Compute the Cholesky decomposition A = LL'.")
	.Category(CATEGORY)
	.Documentation(R"(
<pre>
C  ***  COMPUTE ROWS N1 THROUGH N OF THE CHOLESKY FACTOR  L  OF
C  ***  A = L*(L**T),  WHERE  L  AND THE LOWER TRIANGLE OF  A  ARE BOTH
C  ***  STORED COMPACTLY BY ROWS (AND MAY OCCUPY THE SAME STORAGE).
C  ***  IRC = 0 MEANS ALL WENT WELL.  IRC = J MEANS THE LEADING
C  ***  PRINCIPAL  J X J  SUBMATRIX OF  A  IS NOT POSITIVE DEFINITE --
C  ***  AND  L(J*(J+1)/2)  CONTAINS THE (NONPOS.) REDUCED J-TH DIAGONAL.
</pre>
)")
);
_FPX* WINAPI xll_dl7srt(_FPX* pa, int n1)
{
#pragma XLLEXPORT
	if (n1 == 0) {
		n1 = 1;
	}
	int n = size(*pa);
	int irc;

	DL7SRT(&n1, &n, pa->array, pa->array, &irc);

	return irc == 0 ? pa : nullptr;
}

