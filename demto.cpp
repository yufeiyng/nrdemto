#include "demto.h"
#include "functions.h"
#include "nr_scheme.h"

#include <stdlib.h>
#include <float.h>
#include <stdio.h>

void setBound(double* lower, double* upper, int no_problem) {
	switch (no_problem)
	{

	case 1:	//complete intersection with high similarity, Griewankand Rastrigin
		lower[0] = -100;
		upper[0] = 100;

		lower[1] = -50;
		upper[1] = 50;
		break;

	case 2:	//complete intersection with medium similarity, Ackleyand Rastrigin
		lower[0] = -50;
		upper[0] = 50;

		lower[1] = -50;
		upper[1] = 50;
		break;

	case 3:	//complete intersection with low similarity, Ackleyand Schwefel
		lower[0] = -50;
		upper[0] = 50;

		lower[1] = -50;
		upper[1] = 50;
		break;

	case 4:	//partially intersection with high similarity, Rastriginand Sphere
		lower[0] = -50;
		upper[0] = 50;

		lower[1] = -100;
		upper[1] = 100;
		break;

	case 5:	//partially intersection with medium similarity, Ackleyand Rosenbrock
		lower[0] = -50;
		upper[0] = 50;

		lower[1] = -50;
		upper[1] = 50;
		break;

	case 6:	//partially intersection with low similarity, Ackley and Weierstrass
		lower[0] = -50;
		upper[0] = 50;

		lower[1] = -0.5;
		upper[1] = 0.5;
		break;

	case 7:	//no intersection with high similarity, Rosenbrockand Rastrigin
		lower[0] = -50;
		upper[0] = 50;

		lower[1] = -50;
		upper[1] = 50;
		break;

	case 8: //% no intersection with medium similarity, Griewankand Weierstrass
		lower[0] = -100;
		upper[0] = 100;

		lower[1] = -0.5;
		upper[1] = 0.5;
		break;

	case 9:// % no overlap with low similarity, Rastriginand Schwefel
		lower[0] = -50;
		upper[0] = 50;

		lower[1] = -500;
		upper[1] = 500;
		break;

	default:
		break;
	}
}

double evaluate(double* data, int* evaluations, int t, int no_problem) {
	evaluations[t]++;
	double result = 0.0;
	switch (no_problem)
	{
	case 1:
		if (t == 0)
			result = griewank(data);
		else
			result = rastrigin(data);
		break;
	case 2:
		if (t == 0)
			result = ackley(data);
		else
			result = rastrigin(data);
		break;
	case 3:
		if (t == 0)
			result = ackley(data);
		else
			result = schwefel(data);
		break;
	case 4:
		if (t == 0)
			result = rastrigin(data);
		else
			result = sphere(data);
		break;
	case 5:
		if (t == 0)
			result = ackley(data);
		else
			result = rosenbrock(data);
		break;
	case 6:
		if (t == 0)
			result = ackley(data);
		else
			result = weierstrass(data);
		break;
	case 7:
		if (t == 0)
			result = rosenbrock(data);
		else
			result = rastrigin(data);
	case 8:
		if (t == 0)
			result = griewank(data);
		else
			result = weierstrass(data);
		break;
	case 9:
		if (t == 0)
			result = rastrigin(data);
		else
			result = schwefel(data);
		break;
	default:
		printf("error problem number!");
		break;
	}

	return result;
}

void setSolution(double* solution1, double* solution2) {

	for (int j = 0; j < DIMENSION; j++) {
		solution1[j] = solution2[j];
	}
}

void initializePop(Population* pop, double cur_lower, double cur_upper,  int t) {
	
	pop[t].lower = cur_lower;
	pop[t].upper = cur_upper;
	pop[t].best_index = 0;
	pop[t].best_fitness = DBL_MAX;

	for (int i = 0; i < POP_SIZE; i++) {
		for (int j = 0; j < DIMENSION; j++) {
			pop[t].ind[i][j] = pop[t].lower + (double)rand() / RAND_MAX * (pop[t].upper - pop[t].lower);
		}
	}
}

void evaluatePop(Population* pop, BinarySpacePartition* cur_tree, int* evaluations, int t, int no_problem) {

	for (int i = 0; i < POP_SIZE; i++) {

		BinarySpacePartition* curr_node = searchNode(cur_tree, pop[t].ind[i]);

		if (curr_node->axis == -2) {
			pop[t].fitness[i] = evaluate(pop[t].ind[i], evaluations, t, no_problem);
			insertNode(curr_node, pop[t].ind[i], pop[t].fitness[i]);
			continue;
		}

		BinarySpacePartition* real_node = curr_node;
		while (real_node->flag_vir) {
			real_node = real_node->parent;
		}

		if (euclidean(real_node->data, pop[t].ind[i]) < 1.0 / evaluations[t]) {
			pop[t].fitness[i] = curr_node->fitness;
			updateIndividual(cur_tree, curr_node, pop, evaluations, t, i, no_problem);
		}
		else {
			pop[t].fitness[i] = evaluate(pop[t].ind[i], evaluations, t, no_problem);
			insertNode(curr_node, pop[t].ind[i], pop[t].fitness[i]);
		}
	}
}

void evolution(Population* pop, Population* trial, double f, double cr, int t) {

	for (int i = 0; i < POP_SIZE; i++) {

		int r1, r2, r3;
		do {
			r1 = rand() % POP_SIZE;
			r2 = rand() % POP_SIZE;
			r3 = rand() % POP_SIZE;
		} while (r1 == i || r2 == i || r3 == i || r1 == r2 || r2 == r3);

		int r = rand() % DIMENSION;

		for (int j = 0; j < DIMENSION; j++) {
			if (j == r || ((double)rand() / (RAND_MAX + 1.0)) < cr) {
				trial[t].ind[i][j] = pop[t].ind[r1][j] + f * (pop[t].ind[r2][j] - pop[t].ind[r3][j]);
				if (trial[t].ind[i][j] > pop[t].upper || trial[t].ind[i][j] < pop[t].lower) {
					trial[t].ind[i][j] = pop[t].lower + (double)rand() / RAND_MAX * (pop[t].upper - pop[t].lower);
				}
			}
			else {
				trial[t].ind[i][j] = pop[t].ind[i][j];
			}
		}
	}
}

void selection(Population* pop, Population* trial, int t) {

	for (int i = 0; i < POP_SIZE; i++) {

		if (trial[t].fitness[i] < pop[t].fitness[i]) {
			setSolution(pop[t].ind[i], trial[t].ind[i]);
			pop[t].fitness[i] = trial[t].fitness[i];
		}

		if (pop[t].fitness[i] < pop[t].best_fitness) {
			pop[t].best_index = i;
			pop[t].best_fitness = pop[t].fitness[i];
		}

	}
}

void transfer(Population* pop, int t) {

	int aT = (t + 1) % T;
	int r = rand() % POP_SIZE;
	int best_index = pop[t].best_index;
	setSolution(pop[t].ind[r], pop[aT].ind[best_index]);
}