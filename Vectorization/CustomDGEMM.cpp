#include "CustomDGEMM.h"

void setMtxSize(int& m, int& n, int& k, int& p, int& x, int& y, int& z, int& t)
{
	std::cout << "Enter amount of columns of external matrix A:\n";
	std::cin >> m;
	std::cout << "Enter amount of rows of external matrix A:\n";
	std::cin >> n;
	std::cout << "Enter amount of columns of external matrix B:\n";
	std::cin >> k;
	std::cout << "Enter amount of rows of external matrix B:\n";
	std::cin >> p;
	std::cout << "Enter amount of columns of internal matrix A:\n";
	std::cin >> x;
	std::cout << "Enter amount of rows of internal matrix A:\n";
	std::cin >> y;
	std::cout << "Enter amount of columns of internal matrix B:\n";
	std::cin >> z;
	std::cout << "Enter amount of rows of internal matrix B:\n";
	std::cin >> t;
}

double* mtxGenerate(const int& m, const int& n, const int& x, const int& y) {
	double* mtx = new double[m * n * x * y];
	for (int i = 0; i < m * n * x * y; i++) {
		mtx[i] = (double)rand() / RAND_MAX;
	}
	return mtx;
}

//my own implementation of general matrix multiplication with real numbers
double* naiveCustomVectorizedDGEMM(double* mtx_a, double* mtx_b, int& m, int& n, int& k, int& p, int& x, int& y, int& z, int& t) {
	double* mtx_c = new double[m * p * y * z];
	for (int h1 = 0; h1 < m; h1++) {
		for (int l2 = 0; l2 < p; l2++) {
			//create zero-matix for sum
			double* mtx_sum = new double[x * t];
			for (int o = 0; o < x * t; o++) {
				mtx_sum[o] = 0;
			}
			for (int kn = 0; kn < k; kn++) {
				//multiplexation of internal matrices
				for (int i1 = x; i1 >= 0; i1--) {
					for (int j2 = 0; j2 < t; j2++) {
						double el_sum = 0.0;
						for (int yz = 0; yz < y; yz++) {
							el_sum += mtx_a[h1 * x * y * n + x * y * kn + i1 * y + yz] * mtx_b[kn * p * z * t + z * t * l2 + t * yz + j2];
						}
						mtx_sum[t * i1 + j2] += el_sum;
					}
				}

			}
			for (int i = 0; i < x; i++) {
				for (int j = 0; j < t; j++) {
					mtx_c[h1 * x * t * p + x * t * l2 + t * i + j] = mtx_sum[i * t + j];
				}
			}
		}
	}
	return mtx_c;
}


double* naiveCustomUnvectorizedDGEMM(double* mtx_a, double* mtx_b, int& m, int& n, int& k, int& p, int& x, int& y, int& z, int& t) {
	double* mtx_c = new double[m * p * y * z];
	for (int h1 = 0; h1 < m; h1++) {
		for (int l2 = 0; l2 < p; l2++) {
			//create zero-matix for sum
			double* mtx_sum = new double[x * t];
#pragma loop(no_vector)
			for (int o = 0; o < x * t; o++) {
				mtx_sum[o] = 0.0;
			}
			for (int kn = 0; kn < k; kn++) {
				//multiplexation of internal matrices
				for (int i1 = 0; i1 < x; i1++) {
					for (int j2 = 0; j2 < t; j2++) {
						double el_sum = 0.0;
#pragma loop(no_vector)
						for (int yz = 0; yz < y; yz++) {
							el_sum += mtx_a[h1 * x * y * n + x * y * kn + i1 * y + yz] * mtx_b[kn * p * z * t + z * t * l2 + t * yz + j2];
						}
						mtx_sum[t * i1 + j2] += el_sum;
					}
				}
			}
			for (int i = 0; i < x; i++) {
				for (int j = 0; j < t; j++) {
					mtx_c[h1 * x * t * p + x * t * l2 + t * i + j] = mtx_sum[i * t + j];
				}
			}
		}
	}
	return mtx_c;
}

double* naiveHandmadeVectorizedDGEMM(double* mtx_a, double* mtx_b, const int& m, const int& n, const int& k, const int& p, const int& x, const int& y, const int& z, const int& t) {
	double* mtx_c = new double[m * p * y * z];
	double tmp[4];
	memset(tmp, 0.0, 4);
	memset(mtx_c, 0.0, m * p * y * z);
	for (int h1 = 0; h1 < m; h1++) {
		for (int l2 = 0; l2 < p; l2++) {
			for (int kn = 0; kn < k; kn++) {
				//multiplexation of internal matrices
				for (int i1 = 0; i1 < x; i1++) {
					for (int j2 = 0; j2 < t; j2++) {
						__m256d row_a = _mm256_load_pd(mtx_a + h1 * x * y * n + x * y * kn + i1 * y);
						tmp[0] = mtx_b[kn * p * z * t + z * t * l2 + t * 0 + j2];
						tmp[1] = mtx_b[kn * p * z * t + z * t * l2 + t * 1 + j2];
						tmp[2] = mtx_b[kn * p * z * t + z * t * l2 + t * 2 + j2];
						tmp[3] = mtx_b[kn * p * z * t + z * t * l2 + t * 3 + j2];
						__m256d row_b = _mm256_load_pd(tmp);
						__m256d mul = _mm256_mul_pd(row_a, row_b);
						_mm256_storeu_pd(tmp, mul);
						mtx_c[h1 * x * t * p + x * t * l2 + t * i1 + j2] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
					}
				}

			}
		}
	}
	return mtx_c;
}
