#include <time.h>
#include <chrono>

#include "de.h"

int main() {

	srand((unsigned)time(NULL));

	FILE* fp;
	fopen_s(&fp, "nrde.txt", "w");

	const double f = 0.5;
	const double cr = 0.9;

	const int np = 50;
	const int d = 30;

	const double upper_bound = 50;
	const double lower_bound = -50;

	double** pop = (double**)malloc(np * sizeof(double*));
	double** mut = (double**)malloc(np * sizeof(double*));
	double** trial = (double**)malloc(np * sizeof(double*));
	double* fitness_pop = (double*)malloc(np * sizeof(double));
	double* fitness_trial = (double*)malloc(np * sizeof(double));

	for (int i = 0; i < np; i++) {
		pop[i] = (double*)malloc(d * sizeof(double));
		mut[i] = (double*)malloc(d * sizeof(double));
		trial[i] = (double*)malloc(d * sizeof(double));
	}

	double* critical = (double*)malloc(sizeof(double));	//THE CRITICAL OF REVISITING
	*critical = 100;

	const int max_run_times = 20;
	const int max_evaluations = 30000;

	for (size_t r = 0; r < max_run_times; r++) {

		double* best_solution = (double*)malloc(d * sizeof(double));
		double* best_fitness = (double*)malloc(sizeof(double));
		*best_fitness = DBL_MAX;

		int solution_repeat_times = 0;
		int* evaluations = (int*)malloc(sizeof(int));
		*evaluations = 0;

		BinarySpacePartition* bsp = (BinarySpacePartition*)malloc(sizeof(BinarySpacePartition));
		initialNode(bsp, d);
		bsp->axis = -2; // THE ROOT HAS NOT SOTRE ANY SOLUTION

		int* cur_tree_size = (int*)malloc(sizeof(int));
		*cur_tree_size = 0;

		initialPop(pop, upper_bound, lower_bound, np, d);
		evaluatePop(bsp, pop, fitness_pop, evaluations, np, d);

		while (*evaluations < max_evaluations) {

			bool flag = true;
			double* last_best = (double*)malloc(d * sizeof(double));
			setSolution(last_best, best_solution, d);

			mutation(f, pop, mut, upper_bound, lower_bound, np, d);
			crossover(cr, pop, mut, trial, np, d);
			evaluatePop(bsp, trial, fitness_trial, evaluations, np, d);
			selection(pop, fitness_pop, trial, fitness_trial, best_fitness, best_solution, np, d);
		}

		printf("%.4le\n", *best_fitness);
		fprintf(fp, "%.4le\n", *best_fitness);

		freeBSP(bsp);
		free(best_solution);


	}

	free(pop);

	free(mut);
	free(trial);

	free(fitness_pop);
	free(fitness_trial);

	return 0;
}