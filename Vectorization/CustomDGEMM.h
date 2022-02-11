#pragma once
#include <iostream>
#include <immintrin.h>

void setMtxSize(int& m, int& n, int& k, int& p, int& x, int& y, int& z, int& t);
double* mtxGenerate(const int& m, const int& n, const int& x, const int& y);
double* naiveCustomVectorizedDGEMM(double* mtx_a, double* mtx_b, int& m, int& n, int& k, int& p, int& x, int& y, int& z, int& t);
double* naiveCustomUnvectorizedDGEMM(double* mtx_a, double* mtx_b, int& m, int& n, int& k, int& p, int& x, int& y, int& z, int& t);
double* naiveHandmadeVectorizedDGEMM(double* mtx_a, double* mtx_b, const int& m, const int& n, const int& k, const int& p, const int& x, const int& y, const int& z, const int& t);