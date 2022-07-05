#pragma once

#include "constant.h"

struct Population
{
	double ind[POP_SIZE][DIMENSION];
	double fitness[POP_SIZE];

	int best_index;
	double best_fitness;

	double lower;
	double upper;
};