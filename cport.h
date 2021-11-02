// cport.h - C declarations of PORT functions
#pragma once

extern "C" {
	double D1MACH(int*);
	void DIVSET(int* ALG, int* IV, int* LIV, int* LV, double* V);
	//void DRMNF(double* D, double* FX, int* IV, int* LIV, int* LV, int* N, double* V, double* X);
	void DMNF(int* N, double* D, double* X, void* QF, int* IV, int* LIV, int* LV, double* V, int* UI, double* UR, void* DUMMY);
	void DMNFB(int* N, double* D, double* X, double* B, void* CALCF, int* IV, int* LIV, int* LV, double* V, int* UIPARM, double* URPARM, void* UFPARM);
}
