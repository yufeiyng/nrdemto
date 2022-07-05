#include "nr_scheme.h"
#include "demto.h"

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <queue>

void initialNode(BinarySpacePartition* cur_node) {

	cur_node->axis = -1;
	cur_node->fitness = DBL_MAX;
	cur_node->flag_vir = false;

	cur_node->parent = nullptr;
	cur_node->child_left = nullptr;
	cur_node->child_right = nullptr;
}

BinarySpacePartition* searchNode(BinarySpacePartition* cur_tree, double* cur_data) {

	BinarySpacePartition* cur_node = cur_tree;

	while (cur_node->axis > -1) {
		if (cur_data[cur_node->axis] <= (cur_node->lower + cur_node->upper) / 2.0)
			cur_node = cur_node->child_left;
		else
			cur_node = cur_node->child_right;
	}
	return cur_node;
}

BinarySpacePartition* levelOrder(BinarySpacePartition* cur_tree) {
	std::queue<BinarySpacePartition*> q;
	if (cur_tree != nullptr) {
		q.push(cur_tree);
	}
	while (!q.empty()) {
		if (q.front()->child_left != nullptr) {
			q.push(q.front()->child_left);
		}
		if (q.front()->child_right != nullptr) {
			q.push(q.front()->child_right);
		}
		if (q.front()->axis == -1) return q.front();
		q.pop();
	}
}

void insertNode(BinarySpacePartition* cur_node, double* cur_data, double cur_fitness) {

	if (cur_node->axis == -2) {	
		cur_node->data = (double*)malloc(DIMENSION * sizeof(double));
		setSolution(cur_node->data, cur_data);
		cur_node->fitness = cur_fitness;
		cur_node->axis = -1;
		return;
	}

	double distance[DIMENSION];
	cur_node->axis = 0;
	for (int j = 0; j < DIMENSION; j++) {
		distance[j] = fabs(cur_node->data[j] - cur_data[j]);
		if (distance[j] > distance[cur_node->axis]) {
			cur_node->axis = j;
		}
	}
	cur_node->child_left = (BinarySpacePartition*)malloc(sizeof(BinarySpacePartition));
	cur_node->child_right = (BinarySpacePartition*)malloc(sizeof(BinarySpacePartition));

	initialNode(cur_node->child_left);
	initialNode(cur_node->child_right);

	cur_node->child_left->data = (double*)malloc(DIMENSION * sizeof(double));
	cur_node->child_right->data = (double*)malloc(DIMENSION * sizeof(double));
	
	if (cur_node->data[cur_node->axis] < cur_data[cur_node->axis]) {
		setSolution(cur_node->child_left->data, cur_node->data);
		setSolution(cur_node->child_right->data, cur_data);

		cur_node->child_left->fitness = cur_node->fitness;
		cur_node->child_right->fitness = cur_fitness;
	}
	else {
		setSolution(cur_node->child_left->data, cur_data);
		setSolution(cur_node->child_right->data, cur_node->data);

		cur_node->child_left->fitness = cur_fitness;
		cur_node->child_right->fitness = cur_node->fitness;
	}
	cur_node->child_left->parent = cur_node;
	cur_node->child_right->parent = cur_node;

	cur_node->lower = cur_node->child_left->data[cur_node->axis];
	cur_node->upper = cur_node->child_right->data[cur_node->axis];

	cur_node->flag_vir = true;
	free(cur_node->data);
	cur_node->fitness = DBL_MAX;
}

double euclidean(double* s1, double* s2) {

	double sigma = 0.0;
	for (int j = 0; j < DIMENSION; j++)
		sigma += pow((s1[j] - s2[j]), 2);
	return sqrt(sigma);
}

void updateIndividual(BinarySpacePartition* cur_tree, BinarySpacePartition* cur_node, Population* pop, int* evaluations, int t, int i, int no_problem) {

	BinarySpacePartition* parent = cur_node->parent;

	double new_indiv[DIMENSION];
	for (int j = 0; j < DIMENSION; j++) {
		new_indiv[j] = parent->child_left->data[j] + (double)rand() / RAND_MAX * (parent->child_right->data[j] - parent->child_left->data[j]);
	}

	BinarySpacePartition* new_node = searchNode(cur_tree, new_indiv);
	if (euclidean(new_node->data, new_indiv) < 1.0 / evaluations[t]) {
		return;
	}
	else {
		double new_fitness = evaluate(new_indiv, evaluations, t, no_problem);
		insertNode(new_node, new_indiv, new_fitness);

		setSolution(pop[t].ind[i], new_indiv);
		pop[t].fitness[i] = new_fitness;
	}
}

void freeBSPTree(BinarySpacePartition* cur_node) {

	if (cur_node == NULL) {
		return;
	}
	freeBSPTree(cur_node->child_left);
	freeBSPTree(cur_node->child_right);
	free(cur_node);
}