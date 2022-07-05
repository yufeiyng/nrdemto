#include "functions.h"
#include "constant.h"

#include <math.h>
#include <corecrt_math_defines.h>

double ackley(double* x) {

	double sigma1 = 0.0, sigma2 = 0.0;
	double a = 20.0, b = 0.2, c = 2 * M_PI;

	for (int j = 0; j < DIMENSION; j++) {
		sigma1 += x[j] * x[j];
		sigma2 += cos(c * x[j]);
	}

	double result = -1.0 * a * exp(-1.0 * b * sqrt(sigma1 / (double)DIMENSION)) - exp(sigma2 / (double)DIMENSION) + a + exp(1.0);
	return result;
}

double griewank(double* x) {

	double sigma = 0.0, pi = 1.0;

	for (int j = 0; j < DIMENSION; j++) {
		sigma += (x[j] * x[j]) / 4000.0;
		pi *= cos(x[j] / sqrt(j + 1.0));
	}

	double result = sigma - pi + 1;
	return result;
}

double rastrigin(double* x) {

	double sigma = 0.0;
	for (int j = 0; j < DIMENSION; j++) {
		sigma += x[j] * x[j] - 10.0 * cos(2.0 * M_PI * x[j]);
	}
	double result = 10.0 * DIMENSION + sigma;
	return result;
}

double schwefel(double* x) {

	double sigma = 0.0;
	for (int j = 0; j < DIMENSION; j++) {
		sigma += x[j] * sin(sqrt(fabs(x[j])));
	}
	double result = 418.9829 * (double)DIMENSION - sigma;
	return result;
}

double sphere(double* x) {

	double sigma = 0.0;
	for (int j = 0; j < DIMENSION; j++) {
		sigma += x[j] * x[j];
	}
	return sigma;
}

double rosenbrock(double* x) {

	double sigma = 0.0;
	for (int j = 0; j < DIMENSION - 1; j++) {
		sigma += 100 * (x[j + 1] - x[j] * x[j]) * (x[j + 1] - x[j] * x[j]) + (x[j] - 1.0) * (x[j] - 1.0);
	}
	return sigma;
}

double weierstrass(double* x) {

	double sigma1 = 0.0;
	double sigma2 = 0.0;
	for (int j = 0; j < DIMENSION; j++) {
		for (int k = 0; k < 20; k++) {
			sigma1 += pow(0.5, k) * cos(2 * M_PI * pow(3, k) * (x[j] + 0.5));
		}
	}
	for (int k = 0; k < 20; k++) {
		sigma2 += pow(0.5, k) * cos(2 * M_PI * pow(3, k) * 0.5);
	}

	double result = sigma1 - DIMENSION * sigma2;
	return result;
}