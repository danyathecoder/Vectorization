#pragma once
#include <iostream>

double** initializeMatrix(int n);
void mtxGenerate2D(double** M, int n);
void printMatrix(double** M, int n);
double** add(double** M1, double** M2, int n);
double** subtract(double** M1, double** M2, int n);
double** strassenMultiply(double** A, double** B, int n);