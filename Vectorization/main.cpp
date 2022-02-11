#include <iostream>
#include <chrono>
#include "Strassens.h"
#include "CustomDGEMM.h"

#define N 512
#define X 4

using namespace std;


int main() {
	char c;
	char s;
	int m, n, k, p; // size of external matrices
	int x, y, z, t; // size of internal matrices

	cout << "Do you want to set matrix size by yourself? (y/n)\n";
	cin >> c;
	if (c == 'y') {
		setMtxSize(m, n, k, p, x, y, z, t);
		cout << m << " " << n << " " << k << " " << p << endl;
		cout << x << " " << y << " " << z << " " << t << endl;
		if (n != k || y != z) {
			cout << "Mismatched matrices!";
			return 0;
		}
	}
	else {
		m = n = k = p = N;
		x = y = z = t = X;
	}
	cout << "Do you want test Strassen algorithm? (y/n)\n";
	cin >> s;

	// mtx[l][h][i][j] = arr[h*n*x*y + x*y*l + y*i + j]

	double* mtx_a = mtxGenerate(m, n, x, y);
	double* mtx_b = mtxGenerate(k, p, z, t);

	//1
	chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

	double* unvec_mtx_c = naiveCustomUnvectorizedDGEMM(mtx_a, mtx_b, m, n, k, p, x, y, z, t);

	chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
	chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	cout << "Processing using naive alogrithm with auto vectorization took: " << time_span.count() << " seconds.\n";
	delete[] unvec_mtx_c;
	//2
	t1 = chrono::high_resolution_clock::now();
	double* vec_mtx_c = naiveCustomVectorizedDGEMM(mtx_a, mtx_b, m, n, k, p, x, y, z, t);

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	cout << "Processing using naive alogrithm without auto vectorization took: " << time_span.count() << " seconds.\n";


	//3
	t1 = chrono::high_resolution_clock::now();
	double* custom_vec_mtx_c = naiveHandmadeVectorizedDGEMM(mtx_a, mtx_b, m, n, k, p, x, y, z, t);

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	cout << "Processing using naive alogrithm with custom vectorization took: " << time_span.count() << " seconds.\n";
	for (int i = 0; i < n * m * x * y; i++) {
		if (custom_vec_mtx_c[i] - vec_mtx_c[i] > 0 && custom_vec_mtx_c[i] - vec_mtx_c[i] > 0.001) {
			cout << "Wrong answer!";
		}
	}

	delete[] mtx_a;
	delete[] mtx_b;
	delete[] vec_mtx_c;
	delete[] custom_vec_mtx_c;

	if (s == 'y') {
		cout << "Now let's use Strassens' alogrithm..." << endl;

		double** A = initializeMatrix(N * X);
		double** B = initializeMatrix(N * X);
		mtxGenerate2D(A, N * X);
		mtxGenerate2D(B, N * X);


		double** C = initializeMatrix(N * X);

		t1 = chrono::high_resolution_clock::now();
		C = strassenMultiply(A, B, N * X);

		t2 = chrono::high_resolution_clock::now();
		time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
		cout << "Processing using Strassen's alogrithm with auto vectorization took: " << time_span.count() << " seconds.\n";

		for (int i = 0; i < N * X; i++)
			delete[] A[i];
		delete[] A;

		for (int i = 0; i < N * X; i++)
			delete[] B[i];
		delete[] B;

		for (int i = 0; i < N * X; i++)
			delete[] C[i];
		delete[] C;
	}	

	
	return 0;
}