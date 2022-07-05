#pragma once

#include "bsp_tree.h"
#include "population.h"

void initialNode(BinarySpacePartition* cur_node);

BinarySpacePartition* searchNode(BinarySpacePartition* cur_tree, double* cur_data);

BinarySpacePartition* levelOrder(BinarySpacePartition* cur_tree);

void insertNode(BinarySpacePartition* cur_node, double* cur_data, double cur_fitness);

double euclidean(double* s1, double* s2);

void updateIndividual(BinarySpacePartition* cur_tree, BinarySpacePartition* cur_node, Population* pop, int* evaluations, int t, int i, int no_problem);

void freeBSPTree(BinarySpacePartition* cur_tree);