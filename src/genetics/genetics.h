/*
 * genetics.h
 *
 *  Created on: 18/mag/2015
 *      Author: nicola
 */

#ifndef GENETICS_H_
#define GENETICS_H_

#include "ufl.h"

using namespace std;

Solution generate_ramdom_solution(Instance& inst, int parent_nzcnt ,double sparseness, bool integer);
void selection(vector<Solution>& population, double harshness);
void selection(vector<Solution>& population, int size);
void selection(vector<Solution>& population, int size, int elite);


int mutate(int index, int gene);
int mutate(int index, int gene, Instance& inst, Solution& sol);

Solution reverse(Instance& inst, Solution& solution, bool compute_fitness);

Solution* crossover(Instance& inst, Solution a, Solution b, double mutation, Solution& gene_pool, int n_cut_points);
Solution* crossover(Instance& inst, Solution a, Solution b, double mutation);
Solution* crossover(Instance& inst, Solution a, Solution b, double mutation, Solution& gene_pool);

Solution* crossover(Instance& inst, Solution a, Solution b, double mutation, Solution& gene_pool, int n_cut_points, bool compute_fitness);

Solution* kCutsCrossover(Instance& inst, 	Solution a, 			Solution b, 		double mutation,
											Solution& gene_pool, 	int n_cut_points, 	bool compute_fitness, double& duration);

Solution* kCutsCrossover(Instance& inst, 	Solution a, 			Solution b, 		double mutation,
											Solution& gene_pool, 	int n_cut_points, 	bool compute_fitness);

Solution* simpleCrossover(Instance& inst, Solution a, Solution b, double mutation, Solution& gene_pool);

Solution* slowCrossover(Instance& inst, Solution a, Solution b, double mutation);

void genetic_algorithm(Instance& instance,  Solution& best_so_far);


#endif /* GENETICS_H_ */
