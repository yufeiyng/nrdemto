#pragma once

#include "bsp_tree.h"
#include "population.h"

void setBound(double* lower, double* upper, int no_problem);

double evaluate(double* data, int* evaluations, int t, int no_problem);

void setSolution(double* solution1, double* solution2);

void initializePop(Population* pop, double cur_lower, double cur_upper, int t);

void evaluatePop(Population* pop, BinarySpacePartition* cur_tree, int* evaluations, int t, int no_problem);

void evolution(Population* pop, Population* trial, double f, double cr, int t);

void selection(Population* pop, Population* trial, int t);

void transfer(Population* pop, int t);