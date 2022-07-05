#include "demto.h"
#include "nr_scheme.h"

#include <stdlib.h>
#include <time.h>
#include <random>
#include <fstream>
#include <iostream>

int main() {

	srand((unsigned int)time(NULL));

	std::default_random_engine generator;
	std::cauchy_distribution<double> f(0.3, 0.1);
	std::normal_distribution<double> cr(0.4, 0.1);

	clock_t start, stop;

	int no_problem;

	std::cout << "Please choose the problem number (1~9):\n";
	std::cin >> no_problem;
	double lower[T];
	double upper[T];
	setBound(lower, upper, no_problem);

	start = clock();

	std::ofstream ofs_t1;
	std::ofstream ofs_t2;
	ofs_t1.open("nrmtde_t1.csv", std::ios::out);
	ofs_t2.open("nrmtde_t2.csv", std::ios::out);

	for (size_t r = 0; r < MAX_RUN_TIMES; r++) {

		int repeat_times[T] = { 0,0 };
		double last_indv[T][DIMENSION];
		int* evaluations = (int*)malloc(T * sizeof(int));
		Population* pop = (Population*)malloc(T * sizeof(Population));
		Population* trial = (Population*)malloc(T * sizeof(Population));
		BinarySpacePartition** bsp = (BinarySpacePartition**)malloc(T * sizeof(BinarySpacePartition*));

		for (size_t t = 0; t < T; t++) {
			evaluations[t] = 0;
			bsp[t] = (BinarySpacePartition*)malloc(sizeof(BinarySpacePartition));
			initialNode(bsp[t]);
			bsp[t]->axis = -2;
			initializePop(pop, lower[t], upper[t], t);
			evaluatePop(pop, bsp[t], evaluations, t, no_problem);
			setSolution(last_indv[t], pop[t].ind[pop[t].best_index]);
		}

		while (evaluations[0] < MAX_EVALUATIONS || evaluations[1] < MAX_EVALUATIONS) {
		//for (int g = 0; g < MAX_GENERATIONS; g++) {

			for (int t = 0; t < T; t++) {
				if (evaluations[t] > MAX_EVALUATIONS) continue;
				if ((double)rand() / (RAND_MAX + 1.0) < RMP) {
					transfer(pop, t);
				}
				evolution(pop, trial, f(generator), cr(generator), t);
				evaluatePop(trial, bsp[t], evaluations, t, no_problem);
				selection(pop, trial, t);

				int best_index = pop[t].best_index;

				if (euclidean(last_indv[t], pop[t].ind[best_index]) < 1.0 / evaluations[t]) {
					repeat_times[t] += 1;
				}
				else {
					repeat_times[t] = 0;
					setSolution(last_indv[t], pop[t].ind[best_index]);
				}

				int r = rand() % POP_SIZE;
				
				BinarySpacePartition* best_node = searchNode(bsp[t], pop[t].ind[pop[t].best_index]);
				updateIndividual(bsp[t], best_node, pop, evaluations, t, r, no_problem);
				if (repeat_times[t] > 100) {
					
					BinarySpacePartition* shallow_node = levelOrder(bsp[t]);
					updateIndividual(bsp[t], shallow_node, pop, evaluations, t, r, no_problem); 
					
				}

				// Convergence
				/*if (t == 0) {
					ofs_t1 << pop[t].best_fitness << std::endl;
				}
				else {
					ofs_t2 << pop[t].best_fitness << std::endl;
				}	*/			
			}
		}

		for (int t = 0; t < T; t++) {
			printf("%.4lE\n", pop[t].best_fitness);

			// Boxplot
			if (t == 0) {
				ofs_t1 << pop[t].best_fitness << std::endl;
			}
			else {
				ofs_t2 << pop[t].best_fitness << std::endl;
			}

			freeBSPTree(bsp[t]);
		}
		std::cout << std::endl;

		free(bsp);
		free(pop);
		free(trial);
		free(evaluations);
	}

	stop = clock();
	float cpu_time = (float)(stop - start);
	printf("nrdemto time (CPU) : %f ms \n", cpu_time);

	ofs_t1.close();
	ofs_t2.close();

	return 0;

}