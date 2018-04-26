/*
 * ufl.cpp
 *
 *  Created on: 23/mar/2015
 *      Author: nicola
 */

#include "ufl.h"
#include "genetics.h"
#include "limits"
#include "float.h"
#include <sys/time.h>
#include <math.h>

//#include <ilcplex/cplex.h>
//#include <ilcplex/ilocplex.h>
//#include <ilconcert/ilomodel.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <exception>
#include <stdexcept>

#include <algorithm>

#include "cuda_functions.h"

using namespace std;

#define VERBOSE 					0
#define VERBOSE_CALLBACK			0	//Very verbose (loop couts)
#define VERBOSE_LESS_CALLBACK		0

/*
#define BENDERS_ON					1
#define USERCUT_CALLBACK_ON			1
#define LAZY_CALLBACK_ON 			1
#define HEURISTIC_CALLBACK_ON		1
#define FORCE_LAZY_CONSTRAINTS 		0
*/

bool BENDERS_ON			=			1;
bool HARD_FIXING		=			0;
bool PROXIMITY_SEARCH	=			0;
bool LOCAL_BRANCHING	=			0;

bool GENETICS		=				0;

bool USERCUT_CALLBACK_ON	=		1;

//CAN'T TOUCH THIS
bool LAZY_CALLBACK_ON 	=			1;

bool HEURISTIC_CALLBACK_ON	=		1;
bool FORCE_LAZY_CONSTRAINTS 	=	0;

bool PLAIN_SOLVE	=				0;
#define THRUST_ON					0
bool CUDA_FITNESS_ON			=	1;

//#define OMP_THREADS					8
int N_THREADS = 					8;

/*	Matheuristics	*/
//Function pointer
//int (*cb)(CALLBACK_HEURISTIC_ARGS) = NULL;

/* Matheuristics */

#define USE_BETTER_SOL          1       // Use the better sol. instead of the root node

//Stuff
double TOLERANCE = pow(10, -5);

///////////////////Data structures for DP///////////////////////////

map< Key, double> funcImap;

std::vector<cudaStream_t> streams;

//Cuda Device pointers
double* d_fixc;
double* d_costs;
double* d_t_costs;



/*
 * Unused
 */
double* reductionFunction2(Instance& instance, int i){

	const clock_t begin_time = clock();

	double* minLinkingCosts = new double[instance.n_facilities];

	//Initialize all to -1
	for(int i = 0; i < instance.n_facilities; i++){
		minLinkingCosts[i] = DBL_MAX;
	}

	for(int i = 0;i < instance.n_facilities;i++){
		for(int j = 0 ; j < instance.n_clients 	; j++){

			if(xpos(i, j, &instance) < minLinkingCosts[i]){
				minLinkingCosts[i] = xpos(i, j, &instance);
			}

		}
	}

	cout << "reduction done in: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC *1000 << " ms\n";

	return minLinkingCosts;
}

/*
 * Unused
 */
double reductionFunction (Instance& instance, double *xstar){

	const clock_t begin_time = clock();

	double wtrue = 0.0;
	double minLinkingCosts = DBL_MAX;

	for(int j = 0 ; j < instance.n_facilities ; j++){

		if( fabs(1.0 - xstar[ypos(j, &instance)]) < TOLERANCE ){

			for(int i = 0 ; i < instance.n_clients - 1 ; i++ ){
				if(instance.costs[xpos(i, j, &instance)] < minLinkingCosts){

					minLinkingCosts = instance.costs[xpos(i, j, &instance)];
					cout << "(" << i+1 << "," << j+1 << "): " << minLinkingCosts << endl;

				}
			}

			wtrue = minLinkingCosts;
			minLinkingCosts = DBL_MAX;
		}

	}

	std::cout << "reduction done in: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC *1000 << " ms\n";
	std::cout << "Wtrue = " << wtrue << "\n";

	return wtrue;
}

double calcZ(double* xstar, Instance& inst, int n_fac, int n_cli, double* fixed_costs, double* costs){

	double z = 0;
	double Wi = 0;

	//Fitness computing
	for(int j = 0 ; j < n_fac ; j++){
		if(xstar[j] == 1.0){
			z += inst.fixed_costs[j];
		}

	}

	for(int i = 0 ; i < n_cli ; i++){

		Wi = intFunctionI(inst, xstar, i);
		//cout << "W(" << i << "):" <<  Wi_0;
		//cout << "W(" << i << "):" <<  Wi_1;
		xstar[wpos(i, &inst)] = Wi;
		z += Wi;

	}

	return z;
}

void CPUfitness(vector<Solution>& sol, Instance& inst, double* fixed_costs, double* costs, double& duration){


	struct timeval start, end;

	gettimeofday(&start, NULL);

	for(int i = 0 ; i < sol.size() ; i++ ){
		//sol[i].z_opt = calcZ(sol[i].xStar, inst, fixed_costs, costs);
	}

	gettimeofday(&end, NULL);
	duration += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

}

double DP_IntfunctionI(Instance& instance, double *xstar, int i){

	clock_t begin_time;

	double* ystar = (double*)malloc(sizeof(double)*instance.n_facilities);
	memcpy(ystar, xstar, sizeof(double)*instance.n_facilities);

	if(VERBOSE_CALLBACK)
		begin_time = clock();

	Key k (*ystar, i, instance.n_facilities);

	double minLinkingCosts = DBL_MAX;

/*
	////////////////////////////////////////////////////////////////////////////////////
	minLinkingCosts = DBL_MAX;
	//for every facility
	for(int j = 0 ; j < instance.n_facilities; j++ ){

		//cout << "xstar[ypos(" << j <<  ", &instance)]: " <<  xstar[ypos(j, &instance)] << endl;

		//(Opened facility)
		if( ( ystar[ypos(j, &instance)] ) > TOLERANCE ){

			//cout << "M:(" << i+1 << "," << j+1 << "): " << instance.costs[i*instance.n_facilities + j] << endl;
			if(instance.costs[i*instance.n_facilities + j] < minLinkingCosts){

				minLinkingCosts = instance.costs[i*instance.n_facilities + j];
				//cout << "found minimum: " << minLinkingCosts << endl;
				//cout << "M:(" << i+1 << "," << j+1 << "): " << instance.costs[i*instance.n_facilities + j] << endl;
			}

		}
	}
*/
	//////////////////////////////////////////////////////////////////////////////////


	if( funcImap.count(k) == 0  ){	//IF NOT FOUND

		cout << "cache miss, writing";

		minLinkingCosts = DBL_MAX;
		//for every facility
		for(int j = 0 ; j < instance.n_facilities; j++ ){

			//cout << "xstar[ypos(" << j <<  ", &instance)]: " <<  xstar[ypos(j, &instance)] << endl;

			//(Opened facility)
			if( ( ystar[ypos(j, &instance)] ) > TOLERANCE ){

				//cout << "M:(" << i+1 << "," << j+1 << "): " << instance.costs[i*instance.n_facilities + j] << endl;
				if(instance.costs[i*instance.n_facilities + j] < minLinkingCosts){

					minLinkingCosts = instance.costs[i*instance.n_facilities + j];
					//cout << "found minimum: " << minLinkingCosts << endl;
					//cout << "M:(" << i+1 << "," << j+1 << "): " << instance.costs[i*instance.n_facilities + j] << endl;
				}

			}
		}

//		cout << " funcImap:" << funcImap[k] << " and minconst:  " <<  minLinkingCosts <<endl;
		funcImap.insert(pair<Key, double>(k, minLinkingCosts));
//		cout << "After, funcImap: " << funcImap[k] << " but, minconst:  " <<  minLinkingCosts << endl;
	} else
		cout << "cache hit ";
		//cout << "by the map: " << "xstar: " << xstar << " i: " << i <<funcImap[make_pair(xstar, i)] << endl;
//	} else {
//		if()

//		cout << minLinkingCosts << " funcImap[make_pair(x*:" ;
//		for (long j = 0; j<instance.n_facilities; j++) {
//
//			if ( ystar[ypos(j, &instance)] == 0.0)
//				continue; // skip if zero
//
//			cout << j << " ";
//			//cout << "Y(" << j+1 << ") = " << pop[k].xStar[ypos(j, &instance)] << " ";
//		}
//
//		cout << "," << i << ")]: " << funcImap[make_pair(*ystar, i)] << endl;
//	}

/*
			cout << "minLinkingCosts: " << minLinkingCosts << " funcImap[make_pair(x*:" ;

			for (long j = 0; j<instance.n_facilities; j++) {

				if ( ystar[ypos(j, &instance)] == 0.0)
					continue; // skip if zero

				cout << j << " ";
				//cout << "Y(" << j+1 << ") = " << pop[k].xStar[ypos(j, &instance)] << " ";
			}

			cout << "," << i << ")]: " << funcImap[make_pair(*ystar, i)] << endl;
*/
//			funcImap[make_pair(*ystar, i)] = minLinkingCosts;
//			exit(0);
//	}
//	cout << "********************************************************" << endl ;
//	funcImap[k] = minLinkingCosts;

	if(VERBOSE_CALLBACK){

		std::cout << "f(i) done in: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC * 1000 << " ms\n";
		std::cout << "minLinkCost = " << minLinkingCosts << "\n";
	}

	free(ystar);
	cout << " funcImap:" << funcImap[k] << " and minconst:  " <<  minLinkingCosts <<endl;

	return funcImap[k];
}

double intFunctionI(Instance& instance, double *xstar, int i){

	clock_t begin_time;

	if(VERBOSE_CALLBACK)
		clock_t begin_time;
		begin_time = clock();

	double minLinkingCosts = DBL_MAX;

	////////////////////////WORKING//////////////////////////////
	//for every facility
	for(int j = 0 ; j < instance.n_facilities; j++ ){

		if( ( xstar[ypos(j, &instance)] ) > TOLERANCE ){
			if(instance.costs[i*instance.n_facilities + j] < minLinkingCosts){
				{

					//(Opened facility)
					if(instance.costs[i*instance.n_facilities + j] < minLinkingCosts){

						minLinkingCosts = instance.costs[i*instance.n_facilities + j];
					}
				}
			}
		}
	}

	if(VERBOSE_CALLBACK){
		std::cout << "f(i) done in: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC * 1000 << " ms\n";
		std::cout << "minLinkCost = " << minLinkingCosts << "\n";
	}

	return minLinkingCosts;
}

double IntfunctionI_GOOD(Instance& instance, double *xstar, int i){

	clock_t begin_time;

	if(VERBOSE_CALLBACK)
		begin_time = clock();

	double minLinkingCosts = DBL_MAX;

	//for every facility
	for(int j = 0 ; j < instance.n_facilities; j++ ){

		//cout << "xstar[ypos(" << j <<  ", &instance)]: " <<  xstar[ypos(j, &instance)] << endl;

		//(Opened facility)
		if( ( xstar[ypos(j, &instance)] ) > TOLERANCE ){

			//cout << "M:(" << i+1 << "," << j+1 << "): " << instance.costs[i*instance.n_facilities + j] << endl;
			if(instance.costs[i*instance.n_facilities + j] < minLinkingCosts){

				minLinkingCosts = instance.costs[i*instance.n_facilities + j];
				//cout << "found minimum: " << minLinkingCosts << endl;
				//cout << "M:(" << i+1 << "," << j+1 << "): " << instance.costs[i*instance.n_facilities + j] << endl;
			}

		}
	}

	if(VERBOSE_CALLBACK){

		std::cout << "f(i) done in: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC * 1000 << " ms\n";
		std::cout << "minLinkCost = " << minLinkingCosts << "\n";
	}


	return minLinkingCosts;
}

/*
 * This works ALWAYS for EVERY case
 */
double functionI_GOOD(Instance& instance, double *xstar, int i){

	clock_t begin_time;

	if(VERBOSE_CALLBACK)
		begin_time = clock();

	double minLinkingCosts = DBL_MAX;

	//for every facility
	for(int j = 0 ; j < instance.n_facilities; j++ ){

		//cout << "xstar[ypos(" << j <<  ", &instance)]: " <<  xstar[ypos(j, &instance)] << endl;

		//(Opened facility)
		if( ( xstar[ypos(j, &instance)] ) > TOLERANCE ){

			//cout << "M:(" << i+1 << "," << j+1 << "): " << instance.costs[i*instance.n_facilities + j] << endl;

			if(instance.costs[i*instance.n_facilities + j] < minLinkingCosts)
				minLinkingCosts = instance.costs[i*instance.n_facilities + j];
				//cout << "found minimum: " << minLinkingCosts << endl;
				//cout << "M:(" << i+1 << "," << j+1 << "): " << instance.costs[i*instance.n_facilities + j] << endl;


		}
	}


	if(VERBOSE_CALLBACK){

		std::cout << "f(i) done in: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC * 1000 << " ms\n";
		std::cout << "minLinkCost = " << minLinkingCosts << "\n";
	}


	return minLinkingCosts;
}

double functionI(Instance& instance, double *xstar, int i){

	clock_t begin_time;

	if(VERBOSE_CALLBACK)
		begin_time = clock();

	double minLinkingCosts = DBL_MAX;
	double lp;

	//for every facility
	#pragma omp parallel private(lp)
	{
		lp = DBL_MAX;
		#pragma omp for
		for(int j = 0 ; j < instance.n_facilities; j++ ){

			//cout << "xstar[ypos(" << j <<  ", &instance)]: " <<  xstar[ypos(j, &instance)] << endl;

			//(Opened facility)
			if( ( xstar[ypos(j, &instance)] ) > TOLERANCE ){

				if(instance.costs[i*instance.n_facilities + j] < lp)
					lp = instance.costs[i*instance.n_facilities + j];

			}
		}
		if(lp < minLinkingCosts ){
			#pragma omp critical											//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			{
				if(lp < minLinkingCosts){
					minLinkingCosts = lp;
					//cout << "found minimum: " << minLinkingCosts << endl;
					//cout << "M:(" << i+1 << "," << j+1 << "): " << instance.costs[i*instance.n_facilities + j] << endl;
				}
			}
		}

	}

	if(VERBOSE_CALLBACK){

		std::cout << "f(i) done in: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC * 1000 << " ms\n";
		std::cout << "minLinkCost = " << minLinkingCosts << "\n";
	}



	return minLinkingCosts;		//DP_IntfunctionI(instance, xstar, i);
}

/*
 * Simple func_obj calculator
 */
double calc_obj_val(double* xstar, Instance inst){

	double objVal;

	for(int j = 0; j < inst.n_facilities ; j++){

		if(xstar[ypos(j, &inst)] > TOLERANCE){
			objVal += inst.fixed_costs[j];
		}
	}

	for(int i = 0 ; i < inst.n_clients ; i++){
		objVal += intFunctionI(inst, xstar, i);
	}

	return objVal;
}

int greater_value(const void* a, const void* b){   // comparison function

    const YStar *arg1 = *(YStar**)(a);
    const YStar *arg2 = *(YStar**)(b);

    //std::cout << "Comparing " << " idx " << arg1->index << " with idx  " << arg2->index << std::endl;

    if(fabs((double)(arg1->lpValue - arg2->lpValue)) < TOLERANCE) return 0; // EQUAL
    if(arg1->lpValue < arg2->lpValue) return 1; // for DESC ordering
    return -1; // (arg1.lpValue > arg2.lpValue)
}

int lesser_value(const void* a, const void* b){   // comparison function

    const Map *arg1 = *(Map**)(a);
    const Map *arg2 = *(Map**)(b);

    //std::cout << "Comparing " << " idx " << arg1->index << " with idx  " << arg2->index << std::endl;

    if(fabs((double)(arg1->value - arg2->value)) < TOLERANCE) return 0; // EQUAL
    if(arg1->value < arg2->value) return -1; // for ASC ordering
    return 1; // (arg1.lpValue < arg2.lpValue)
}

bool compare_sol(Solution a, Solution b){
	if(a.z_opt < b.z_opt){
		return true;
	}else {
		return false;
	}
}


int lesser_sol_value(const void* a, const void* b){   // comparison function

    const Solution *arg1 = *(Solution**)(a);
    const Solution *arg2 = *(Solution**)(b);

    //std::cout << "Comparing " << " idx " << arg1->index << " with idx  " << arg2->index << std::endl;

    if(fabs((double)(arg1->z_opt - arg2->z_opt)) < TOLERANCE) return 0; // EQUAL
    if(arg1->z_opt < arg2->z_opt) return -1; // for ASC ordering
    return 1; // (arg1.lpValue < arg2.lpValue)
}


void init(Instance& inst) {

	if(CUDA_FITNESS_ON){

		//Loading common data structures on GPU
		cout << "CudaFitness: Loading common data structures on GPU" << endl;

		cudaMalloc((void**) &d_fixc	, sizeof(double)*inst.n_facilities);
		CUDAErrorCheck();
		cudaMalloc((void**)	&d_costs	, sizeof(double)*inst.n_facilities*inst.n_clients);
		CUDAErrorCheck();
//		cudaMalloc((void**) &d_t_costs	, sizeof(double)*instance.n_facilities*instance.n_clients);
//		CUDAErrorCheck();

		cudaMemcpy(	d_fixc, inst.fixed_costs, sizeof(double)*inst.n_facilities, cudaMemcpyHostToDevice);
		CUDAErrorCheck();
		cudaMemcpy( d_costs, inst.costs, inst.n_facilities*inst.n_clients*sizeof(double), cudaMemcpyHostToDevice );
		CUDAErrorCheck();

		cudaThreadSynchronize();
	}


    if(BENDERS_ON){
    	printf("Using Benders' decomposition\n");
    } else {
    	printf("Using Standard model\n");
    }

}

int count_non_zero(Instance& inst, Solution& sol){

	int nz_cnt = 0;

	for (long j = 0; j < inst.n_facilities; j++){
		if ( sol.xStar[ypos(j, &inst)] == 0.0)
			continue; // skip if zero
		nz_cnt++;
	}

	return nz_cnt;
}

void prettyPrintXStar(double * xStar, int n_fac){

	for (long j = 0; j<n_fac; j++) {

		if ( xStar[j] == 0.0)
			continue; // skip if zero

		cout << j << " ";
	}
	cout << endl;
	//cout << " : " << sol.z_opt << endl;
}

void prettyPrintXStarWithValue(double * xStar, int n_fac){

	for (long j = 0; j<n_fac; j++) {

		if ( xStar[j] == 0.0)
			continue; // skip if zero

		cout << j << ":" << xStar[j] << " ";
	}
	//cout << " : " << sol.z_opt << endl;
}

void prettyPrint(Solution& sol, Instance& instance){

	for (long j = 0; j<instance.n_facilities; j++) {

		if ( sol.xStar[ypos(j, &instance)] == 0.0)
			continue; // skip if zero

		cout << j << " ";
	}
	cout << " : " << sol.z_opt << endl;
}

bool sameAs(Instance inst, Solution a, Solution b){

	if(a.non_zero_count != b.non_zero_count)
		return false;
	else{
		for(int j = 0 ; j < inst.n_facilities ; j++){
			if(a.xStar[j] != b.xStar[j])
				return false;
		}
		return true;
	}

}

void genetic_algorithm(Instance& instance,  Solution& best_so_far){

	//double non_zero_ratio = 0.5;

/*
	best_so_far.xStar = malloc(sizeof(double)*ncols);
	best_so_far.z_opt = (double)malloc(sizeof(double));
	int ncols = instance.n_facilities + instance.n_clients;
*/

	best_so_far = generate_ramdom_solution(instance, 1, 1, 1);

	cout << "zBest: " << best_so_far.z_opt << endl;

	for (long j = 0; j<instance.n_facilities; j++) {

		if ( best_so_far.xStar[ypos(j, &instance)] == 0.0)
			continue; // skip if zero

		cout << "Y(" << j+1 << ") = " << best_so_far.xStar[ypos(j, &instance)] << endl;
	}

	srand(54321);
	//srand(time(NULL));

	///////////////////////////////////////////	////////////////////////////////////////////////////////////////////////////


	printf("starting genetic algorithm . . .\n");

	//Genetic Algorithm params
	int n_ind = 2500;							//Max number of individuals
	int elite = n_ind*0.20;
	double mutation = 2;						//Average number of facilities that flip bit on mutation
	int DOOMSDAY_TIMER = 300;//50;//100;		//Number of unsuccesful iterations before mass extinction
	int farrows = 210*5;						//Farrow = "figliata", 2 new sons for every farrow

	//Data structures & internal vars
	Solution gene_pool;
	int n_iter = 1;
	gene_pool.xStar = (double*) calloc(instance.n_facilities, sizeof(double));
	vector<Solution> pop;
	double* fitnesses;
	bool new_sol = false;
	int nz_cnt = 0;
	int MAX_ITERATIONS = 50;
	int doomsday_clock = DOOMSDAY_TIMER;

	//Time measurement
	double total_time = 0.0;
	double cuda_total = 0.0f;
	double avg_time = 0;
	double max_time = 0;
	double min_time = FLT_MAX;
	struct timeval start, end;
	double cpu_duration = 0;
	double cuda_avg = 0.0f;

	//INIT
	int iteration = 1;
	pop.push_back(best_so_far);
	cout << "root node zBest: " << best_so_far.z_opt << endl;
	for (long j = 0; j<instance.n_facilities; j++){
		if ( best_so_far.xStar[ypos(j, &instance)] == 0.0)
			continue; // skip if zero
		nz_cnt++;
		cout << j+1 << " ";
	}
	cout << endl;
	pop[0].non_zero_count = nz_cnt;
	pop[0].type = 'O';


	while( iteration < MAX_ITERATIONS ){//doomsday_clock > 100 ){//&& i < 5 ){//i < 5) {//

		gettimeofday(&start, NULL);

		cpu_duration = 0;

		//Mass extinction
		if(doomsday_clock < 1){
			//selection(pop, 10, 5);	//10,3
			doomsday_clock = DOOMSDAY_TIMER;
			srand(12345);	//
			//srand(time(NULL));
		}

		int n_rand_ind = 0;

		//Fill with ramdom individuals
		if( iteration == 1)
			printf("Populating with random individuals\n");

		while(n_ind > pop.size()){

			Solution rand_individual = generate_ramdom_solution(instance, pop[0].non_zero_count, 1, false);//0, 1, true); //pop[0].non_zero_count, 1, true);

			bool already_present = false;

			for(int k = 0 ; k < pop.size() ; k++){
				if(sameAs(instance, rand_individual, pop[k])){
					already_present = true;
					//cout << "%";
				}
			}

			if(!already_present){
				rand_individual.type = 'R';
				pop.push_back(rand_individual);
				n_rand_ind++;
				//cout << "newcomer cameth\n" << endl;
			}
			//cout << pop.size() << ",";

			if( iteration < 2 ){
				if(pop.size()%25 == 0 ){
					clog << "*";
				}
			}
		}
		if(iteration == 1)
			clog << "DONE!" << endl;

		if(pop[0].z_opt + TOLERANCE < best_so_far.z_opt){

			printf("trovata soluzione migliore\n");
			cout << pop[0].z_opt << endl;

			for (long j = 0; j<instance.n_facilities; j++) {

				if ( pop[0].xStar[ypos(j, &instance)] == 0.0)
					continue; // skip if zero

				cout << "Y(" << j << ") = " << pop[0].xStar[ypos(j, &instance)] << endl;
			}

			best_so_far = pop[0];
			doomsday_clock = DOOMSDAY_TIMER;

		} else
			doomsday_clock--;

//          for(int k = 0 ; k < pop.size() ; k++ ){
//				cout << k << "-esimo nodo zBest: " << pop[k].z_opt << endl;
//              for (long j = 0; j<instance.n_facilities; j++)
//                  cout << pop[k].xStar[ypos(j, &instance)] << ", ";
//              cout << endl;
//          }

		////////////////////SORTING//////////////////////////////

		Solution* children = new Solution[2];
		vector<Solution> sons;

		std::sort(pop.begin(), pop.end(), compare_sol);

		///////////////////SELECTION////////////////////////////
		int deaths = 0;

		selection(pop, n_ind/2, elite);

//		pop.resize(n_ind/2);
//
//         while(pop.size()>floor(n_ind/2)){
//
//              int idx = (rand()%(pop.size()));
//              if( idx < n_ind/2  ) //&& rand() < RAND_MAX*0.2)
//                  continue;
//
//              //if(pop[idx].age > 5 || idx > elite ){
//                  pop.erase(pop.begin() + idx);
//                  deaths++;
//              //}
////                idx = (rand()%(pop.size()-elite))+elite;
////                if(pop[idx].non_zero_count > pop[1].non_zero_count*4 ){
////                    pop.erase(pop.begin() + idx);
////                    deaths++;
////                    }
//          }
//          random_shuffle(pop.begin(), pop.end());

		///////////////////CROSSOVER//////////////////



		struct timeval cpu_start, cpu_end;
		gettimeofday(&cpu_start, NULL);

		#pragma omp parallel for
		for(int k = 1 ; k < farrows ; k++ ){

			int idx1 = rand()%(pop.size());		//2*(rand()%(pop.size()/2));
			int idx2 = rand()%(pop.size()); 	//2*(rand()%(pop.size()/2)) + 1;

//					//Alternative breeding strategy
//					int idx1 = rand()%(pop.size();
//					int idx2 = rand()%(pop.size()/2) ;
//
//
//					int idx1 = 2*k;
//					int idx2 = 2*k + 1;
																					////////////////!!!!!!!!!!!!/////////////////
			children = kCutsCrossover(instance, pop[ idx1 ], pop[ idx2 ], mutation, gene_pool, pop[0].non_zero_count, !CUDA_FITNESS_ON);

			#pragma omp critical
			{
				children[0].type = 'A';
				children[1].type = 'A';
				if(children[0].non_zero_count < 128)
					sons.push_back(children[0]);
				else{
					cout << "Errore del trick" << endl;
					exit(0);
				}
				if(children[1].non_zero_count < 128)
					sons.push_back(children[1]);
				else{
					cout << "Errore del trick" << endl;
					exit(0);
				}
			}

		}


		#pragma omp parallel for
		for(int k = 1 ; k < int((farrows/5)) ; k++ ){

			int idx1 = 2*(rand()%(elite/2));
			int idx2 = 2*(rand()%(elite/2)) + 1;

//					//Alternative breeding strategy
//					int idx1 = rand()%(pop.size();
//					int idx2 = rand()%(pop.size()/2) ;
//
//					int idx1 = 2*k;
//					int idx2 = 2*k + 1;
																				//////////////////!!!!!!!!!!!!!!!///////////////////
			children = kCutsCrossover(instance, pop[ idx1 ], pop[ idx2 ], mutation, gene_pool, pop[0].non_zero_count, !CUDA_FITNESS_ON);

			#pragma omp critical
			{
				children[0].type = 'E';
				children[1].type = 'E';
				if(children[0].non_zero_count < 128)
					sons.push_back(children[0]);
				else{
					cout << "Errore del trick" << endl;
					exit(0);
				}
				if(children[1].non_zero_count < 128)
					sons.push_back(children[1]);
				else{
					cout << "Errore del trick" << endl;
					exit(0);
				}
			}

		}

		gettimeofday(&cpu_end, NULL);
		cpu_duration += ((cpu_end.tv_sec  - cpu_start.tv_sec) * 1000000u + cpu_end.tv_usec - cpu_start.tv_usec) / 1.e6;


//		for(int k = 1 ; k < 30 ; k++ ){
//
//			int idx1 = 0;
//			int idx2 = 0;
//
////					//Alternative breeding strategy
////					int idx1 = rand()%(pop.size();
////					int idx2 = rand()%(pop.size()/2) ;
////
////					int idx1 = 2*k;
////					int idx2 = 2*k + 1;
//
//			children = kCutsCrossover(instance, pop[ idx1 ], pop[ idx2 ], 6, gene_pool, pop[idx1].non_zero_count, !CUDA_FITNESS_ON);
//
//			//#pragma omp critical
//			{
//				children[0].type = 'o';
//				children[1].type = 'o';
//				sons.push_back(children[0]);
//				sons.push_back(children[1]);
//			}
//
//		}

/*
		//All with all ( alternatively all elite with all elite )
		//#pragma omp parallel for num_threads(8)
		for(int p = 0 ; p < pop.size() ; p++ ){
			for(int m = 0 ; m < pop.size() ; m++){

				children = crossover(instance, pop[ m ], pop[ p ], mutation, gene_pool);

				//#pragma omp critical
				{
					sons.push_back(children[0]);
					sons.push_back(children[1]);
				}

			}
		}
*/

		double cuda_duration;

		if(CUDA_FITNESS_ON){

			struct timeval cuda_start, cuda_end;
			gettimeofday(&cuda_start, NULL);

			//cout << "sons.size()" << sons.size() << endl;
			fitnesses = (double*)calloc(sons.size(), sizeof(double)); //(double*)malloc(sizeof(double)*n_ind);
			cudaFitness2(sons, fitnesses, d_fixc, d_costs, instance.n_facilities, instance.n_clients);

			gettimeofday(&cuda_end, NULL);
			cuda_duration = ((cuda_end.tv_sec  - cuda_start.tv_sec) * 1000000u + cuda_end.tv_usec - cuda_start.tv_usec) / 1.e6;

/*			for(int k = 0 ; k < sons.size() ; k++ ){

				cout << k << ")" << "GPU fitness: " << fitnesses[k] << " CPU fitness: " << sons[k].z_opt << endl;

			}
*/

			free(fitnesses);

		} /*else {

			//fitnesses = (double*)calloc(sons.size(), sizeof(double)); //(double*)malloc(sizeof(double)*n_ind);
			CPUfitness(sons, instance, cpu_duration);

		}*/


		//Add sons to population only if there isn't another identical
		int sons_survived = 0;
		int aliases = 0;

		for(int s = 0 ; s < sons.size() ; s++){

			if(sons[s].non_zero_count > 128){
				cout << "Cuda: Errore del kernel";
				exit(0);
			}

			bool already_present = false;
			//cout << endl;
			for(int k = 0 ; k < pop.size() ; k++){
				if(sameAs(instance, sons[s], pop[k])){
					already_present = true;
					aliases++;
					break;
				}
			}

			if(!already_present ){//&& sons[s].non_zero_count < pop[1].non_zero_count*3 ){

				pop.push_back(sons[s]);
				sons_survived++;
			} else {
				//If already present add a random individual
				Solution random = generate_ramdom_solution(instance, pop[0].non_zero_count, 1, true);//reverse(instance, pop[0], true);//
				random.type = 'R';
				pop.push_back( random );
				n_rand_ind++;
			}
		}

		gettimeofday(&end, NULL);
		double duration = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

		if(duration > max_time)
			max_time = duration;
		if(duration < min_time)
			min_time = duration;

		total_time = total_time + duration;
		cuda_total = cuda_total + cuda_duration;
		avg_time = (double) total_time/iteration;
		cuda_avg = (double) cuda_total/iteration;

		/////////////////////////////////DEBUG//BLOCK///////////////////////////////////////////
		cout << endl << "Iteration n. "<< iteration << " champion " << pop[0].z_opt << " " << pop[0].type <<  endl;
		cout << "Aliases: " << aliases << " Number of Sons: " << sons.size() << " Newcomers: " << n_rand_ind << " Extinction countdown: " << doomsday_clock <<endl;
		printf("All done in: %3.2f, ", duration);
		printf("Avg. time per iteration: %3.2f, Min: %3.2f, Max: %3.2f\n", avg_time, min_time, max_time);

		if(!CUDA_FITNESS_ON){
			printf("CPU took: %3.3f\n", cpu_duration);
		} else {
			printf("GPU took: %3.3f, ", cuda_duration);
			printf("GPU average: %3.3f\n", cuda_avg);
		}

		for( int k = 0 ; k < 1; k++){//pop.size(); k++){//1 ; k++ ){			//pop.size(); k++){

			cout << "Ind:" << k << " >> ";

			for (long j = 0; j<instance.n_facilities; j++) {

				if(pop[k].non_zero_count > 128){
					cout << "Cuda: Errore del kernel";
					exit(0);
				}

				if ( pop[k].xStar[ypos(j, &instance)] == 0.0)
					continue; // skip if zero

				cout << j << " ";
				//cout << "Y(" << j+1 << ") = " << pop[k].xStar[ypos(j, &instance)] << " ";
			}
			cout << " : " << pop[k].z_opt << " " << pop[k].type << endl;
		}
		///////////////////////////////////////////////////////////////////////////////////////////////

		//Gene pool check
//				for(int j = 0 ; j < instance.n_facilities; j++){
//					cout << j << ": " << gene_pool.xStar[j] << endl;
//				}

		iteration++;

		//END ITERATION
	}

}

void solve_problem (Instance& instance, Solution& sol){

	int error = 0;

	BENDERS_ON 				=	true;
	GENETICS 				= 	true;

	//Load instance in GPU memory
	init(instance);

	int ncols;

	if(BENDERS_ON)
		ncols = instance.n_facilities + instance.n_clients;
	else
		ncols = instance.n_facilities + instance.n_facilities*instance.n_clients;

	Solution best_so_far;
	best_so_far.xStar = (double*)malloc(sizeof(double)*ncols);
	best_so_far.z_opt = DBL_MAX;


	if(GENETICS){	//USE_GENETICS

		double non_zero_ratio = 0.5;

		genetic_algorithm(instance, best_so_far);

		cout << "Final: " << best_so_far.z_opt << endl;

		exit(0);
/*

		////////////////////////////////////////

		error = CPXmipopt(env, lp);

		///////////////////////////////////////////

		//Generate random solution
		//Third recipe

		int n_tries = 1000;

		for(int k = 0 ; k < n_tries ; k++){

			//Generate random sol
			Solution rand_sol;
			rand_sol.xStar = (double*)malloc(sizeof(double)*ncols);
			rand_sol.z_opt = DBL_MAX;

			for( int j = 0 ; j < instance.n_facilities ; j++ ){

				//Coin toss
				if( rand() < RAND_MAX*non_zero_ratio ){

					rand_sol.xStar[ypos(j, &instance)] = 1;
					rand_sol.z_opt += instance.fixed_costs[j];



				}

			}

			for (int i = 0 ; i < instance.n_clients ; i++){

				rand_sol.z_opt += functionI(instance, rand_sol.xStar, i);

			}

			cout << "random zBest: " << rand_sol.z_opt << endl;

			for (long j = 0; j<instance.n_facilities; j++) {

				if ( rand_sol.xStar[ypos(j, &instance)] == 0.0)
					continue; // skip if zero

				cout << "Yrandom(" << j+1 << ") = " << rand_sol.xStar[ypos(j, &instance)] << endl;
			}

			//////////////////////


			if(CPXaddmipstarts(env, lp, 1, nzcnt, beg, varindices, values, effortlevel , NULL))
				throw std::runtime_error("Cannot get add start.");


			CPXsetintparam(env, CPX_PARAM_NODELIM, 10);
			CPXsetdblparam(env, CPX_PARAM_TILIM, 10.0);

			error = CPXmipopt(env, lp);         // RE-SOLVE the model

			//////////////////////////////////

			if (CPXgetx(env, lp, rand_sol.xStar, 0, CPXgetnumcols(env, lp)-1))
				throw std::runtime_error("Cannot get solution variables.");

			if (CPXgetobjval(env, lp, &rand_sol.z_opt))
					throw std::runtime_error("Cannot get the solution.");



			cout << "actual zBest: " << rand_sol.z_opt << endl;

			for (long j = 0; j<instance.n_facilities; j++) {

				if ( rand_sol.xStar[ypos(j, &instance)] == 0.0)
					continue; // skip if zero

				cout << "Yrand(" << j+1 << ") = " << rand_sol.xStar[ypos(j, &instance)] << endl;
			}

			if(rand_sol.z_opt < best_so_far.z_opt){

				printf("New heuristic found better solution!\n");

				best_so_far = rand_sol;

				int nzcnt = 0;
				vector<int> var_indices;
				vector<double> _values;
				int* varindices;
				double* values;

				for( int j = 0 ; j < instance.n_facilities ; j++ ){

					//Coin toss
					if( rand() < RAND_MAX*non_zero_ratio ){//rand()%2 == 1 ){

						//cout << "Lancio " << j << ": Testa" << endl;
						//rand_sol[ypos(j, &instance)] = 1;
						nzcnt++;
						var_indices.push_back(j);
						_values.push_back(1.0);

					} else {

						//cout << "Lancio " << j << ": Croce" << endl;
						//nzcnt++;
						//_values.push_back(0.0);
					}


					} else {

						//cout << "Facility " << j <<" gia' attiva, non lancio niente" << endl;
						nzcnt++;
						var_indices.push_back(j);
						_values.push_back(1.0);
					}

				}

				varindices = &var_indices[0];
				values = &_values[0];
				int *beg = (int*)malloc(sizeof(int)*1);//new int[1];
				beg[0] = 0;

				cout << "Sending solution" << endl;

				int* effortlevel = (int*)malloc(sizeof(int)*1);//new int[1];
				effortlevel[0] = CPX_MIPSTART_AUTO;



				cout << "Solution sent" << endl;


			}

			if (CPXgetx(env, lp, best_so_far.xStar, 0, CPXgetnumcols(env, lp)-1))
				throw std::runtime_error("Cannot get solution variables.");

			cout << "zBest: " << best_so_far.z_opt << endl;

			for (long j = 0; j<instance.n_facilities; j++) {

				if ( best_so_far.xStar[ypos(j, &instance)] == 0.0)
					continue; // skip if zero

				cout << "Y(" << j+1 << ") = " << best_so_far.xStar[ypos(j, &instance)] << endl;
			}

			//cout << "exiting" << endl;
			//exit(0);

		}
*/

		//GENETICS END

	}

	/* Add a callback for Matheuristic, wrap the existing one */
}

