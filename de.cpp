#include "de.h"

double evaluate(double* solution, int dimension, int* evaluations) {

	(*evaluations)++;
	return generailizedpenalized2(solution, dimension);
}

void setSolution(double* solution1, double* solution2, int dimension) {

	for (int j = 0; j < dimension; j++) {
		solution1[j] = solution2[j];
	}
}

void initialPop(double** pop, double upper_bound, double lower_bound, int poplation_size, int dimension) {

	for (int i = 0; i < poplation_size; i++) {
		for (int j = 0; j < dimension; j++) {
			pop[i][j] = lower_bound + (double)rand() / RAND_MAX * (upper_bound - lower_bound);
		}
	}
}

void mutation(double f, double** pop, double** mut, double upper_bound, double lower_bound, int population_size, int dimension) {
	
	int r1, r2, r3;
	for (int i = 0; i < population_size; i++) {

		do {
			r1 = rand() % population_size;
			r2 = rand() % population_size;
			r3 = rand() % population_size;
		} while (r1 == i || r2 == i || r3 == i || r1 == r2 || r2 == r3);

		for (int j = 0; j < dimension; j++) {
			mut[i][j] = pop[r1][j] + f * (pop[r2][j] - pop[r3][j]);
			
			//	边界处理方案来自于 Differential Evolution: A Practical Approach to Global Optimization
			if (mut[i][j] < lower_bound)
				mut[i][j] = (lower_bound + pop[i][j]) / 2.0;
			if (mut[i][j] > upper_bound)
				mut[i][j] = (upper_bound + pop[i][j]) / 2.0;
		}
	}
}

void crossover(double cr, double** pop, double** mut, double** trial, int population_size, int dimension) {
	
	for (int i = 0; i < population_size; i++) {

		int r = rand() % dimension;
		for (int j = 0; j < dimension; j++) {

			if (j == r || (double)rand() / RAND_MAX <= cr)
				trial[i][j] = mut[i][j];
			else
				trial[i][j] = pop[i][j];
		}
	}
}

void evaluatePop(BinarySpacePartition* bsp, double** trial, double* fitness_trial, int* evaluations, int np, int dimension) {
	
	for (int i = 0; i < np; i++) {
		BinarySpacePartition* node = searchNode(bsp, trial[i]);
		if (node->axis == -2) {
			fitness_trial[i] = evaluate(trial[i], dimension, evaluations);
			setSolution(node->solution, trial[i], dimension);
			node->fitness = fitness_trial[i];
			continue;
		}
		if (differSolution(node->solution, trial[i], dimension) > 100) {
			fitness_trial[i] = evaluate(trial[i], dimension, evaluations);
			insertNode(node, trial[i], fitness_trial[i], dimension);
		}
		else {
			fitness_trial[i] = node->fitness;
		}
	}
}

void selection(double** pop, double* fitness_pop, double** trial, double* fitness_trial, 
	           double* global_fitness, double* global_solution, int population_size, int dimension) {
	
	for (int i = 0; i < population_size; i++) {

		if (fitness_trial[i] < fitness_pop[i]) {
			setSolution(pop[i], trial[i], dimension);
			fitness_pop[i] = fitness_trial[i];
		}

		if (fitness_pop[i] <= *global_fitness) {
			setSolution(global_solution, pop[i], dimension);
			*global_fitness = fitness_pop[i];
		}
	}
}

void printPop(double** pop, int population_size, int dimension) {
	for (int i = 0; i < population_size; i++) {
		for (int j = 0; j < dimension; j++) {
			printf("%.4f ", pop[i][j]);
		}
		printf("\n");
	}
}

void randomReplace(double** mutant, double* data, int population_size, int dimension) {

	int rnd = rand() % population_size;
	setSolution(mutant[rnd], data, dimension);
}
