#pragma once
#include <stdio.h>

#include "benchmark.h"
#include "non_revisit.h"
#include "bsp_tree.h"

double evaluate(double* solution, int dimension, int* evaluations);

void setSolution(double* solution1, double* solution2, int dimension);

void initialPop(double** pop, double upperBound, double lowerBound, int poplulation_size, int dimension);

void mutation(double f, double** pop, double** mut, double upper_bound, double lower_bound, int population_size, int dimension);

void crossover(double cr, double** pop, double** mut, double** trial, int population_size, int dimension);

void evaluatePop(BinarySpacePartition* bsp, double** pop, double* fitness_trial, int* evaluations, int np, int dimension);

void selection(double** pop, double* fitness_pop, double** trial, double* fitness_trial, double* global_fitness, double* global_solution, int population_size, int dimension);

void printPop(double** pop, int population_size, int dimension);

void randomReplace(double** mutant, double* data, int population_size, int dimension);