/*
 * ufl.h
 *
 *  Created on: 23/mar/2015
 *      Author: nicola
 */


//#include <ilcplex/cplex.h>
//#include <ilcplex/ilocplex.h>
//#include <ilconcert/ilomodel.h>

#include <string.h>

#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>

#include <map>
#include <utility>
#include <algorithm>

#define USERCALLBACK_MAX_ITER   40000 // or INTMAX
#define VERBOSE					0

#ifndef UFL_H_
#define UFL_H_

using namespace std;

class Key
{
  public:
    Key(double xstar, int i, int n_fac)
    {
      this->xstar = xstar;
      this->i = i;
      this->n_fac = n_fac;
    }
    double xstar;
    int i;
    int n_fac;
	bool operator< (const Key& a) const
	{
		int dbl_cmp = memcmp((void*) (&this->xstar), (void*)(&a.xstar), n_fac*sizeof(double) );
		if(dbl_cmp == 0){
			return this->i < a.i;
		}
		return dbl_cmp < 0;
	}
};

struct Instance {
    long n_facilities;
    long n_clients;
    double *fixed_costs;
    double *costs; // c[i,j] = inst->cost[i*inst->n_facilities+j]
    double *t_costs;
    double zBest;
    double* xBest;

    double giveUp;

    bool hammingAvailable;
    bool isBestAvailable;
    double rootLowerBound;

    double bestHamming;
};

struct Solution {
    double z_opt;
    double *xStar;
    int non_zero_count;
    int age;
    char type;
};

struct YStar {
	int index;
	double lpValue;
};

struct Map {
	int index;
	double value;
};

inline int ypos(int j, Instance *inst) {
    return j;
}

inline int wpos(int i, Instance *inst) {
	return inst->n_facilities + i;
}

inline int xpos(int i, int j, Instance *inst) {
    return inst->n_facilities + i*inst->n_facilities + j;
}

int count_non_zero(Instance inst, Solution sol);

void prettyPrintXStar(double * xStar, int n_fac);
void prettyPrint(Solution& sol, Instance& instance);

double calcZ(double* xstar, Instance& inst, int n_fac, int n_cli, double* fixed_costs, double* costs);

double intFunctionI(Instance& instance, double *xstar, int i);

void solve_problem (Instance& instance, Solution& sol);

#endif /* UFL_H_ */
