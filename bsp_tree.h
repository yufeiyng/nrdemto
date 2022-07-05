#pragma once

struct BinarySpacePartition
{
	double* data;
	double fitness;

	int axis;
	double lower;
	double upper;

	bool flag_vir;

	BinarySpacePartition* parent, * child_left, * child_right;
};