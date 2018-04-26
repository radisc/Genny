#include <thrust/sort.h>
#include <thrust/transform.h>
#include <thrust/host_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/device_allocator.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <algorithm>
#include <cstdlib>
#include <float.h>
#include <ufl.h>

#include <cuda.h>
#include <cuda_device_runtime_api.h>
#include <cuda_runtime_api.h>


//struct YStar {
//	int index;
//	double lpValue;
//};
//
//struct Map {
//	int index;
//	double value;
//};
//
//struct Instance {
//    long n_facilities;
//    long n_clients;
//    double *fixed_costs;
//    double *costs; // c[i,j] = inst->cost[i*inst->n_facilities+j]
//};


//Used to check if there are any errors launching the kernel
void CUDAErrorCheck(){
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess){
		printf("CUDA error : %s (%d)\n", cudaGetErrorString(error), error);
		exit(0);
	}
}

//inline int ypos(int j, const Instance *inst) {
//    return j;
//}
//

//__host__ __device__
//double cudaFunctionI(const double *xstar, int i, 	const thrust::device_vector< double >&	d_fixcs	,
//													const thrust::device_vector< double >& 	d_costs	,
//													const int	d_nfac	,
//													const int	d_ncli	){
//
//	double minLinkingCosts = DBL_MAX;
//
//	//for every facility
//	for(int j = 0 ; j < d_nfac; j++ ){
//
//		//cout << "xstar[ypos(" << j <<  ", &instance)]: " <<  xstar[ypos(j, &instance)] << endl;
//
//		//(Opened facility)
//		if( ( xstar[j] ) > 1e-5 ){
//
//			//cout << "M:(" << i+1 << "," << j+1 << "): " << instance.costs[i*instance.n_facilities + j] << endl;
//
//			if(d_costs[i*d_nfac + j] < minLinkingCosts)
//				minLinkingCosts = d_costs[i*d_nfac + j];
//				//cout << "found minimum: " << minLinkingCosts << endl;
//				//cout << "M:(" << i+1 << "," << j+1 << "): " << instance.costs[i*instance.n_facilities + j] << endl;
//
//
//		}
//	}
//
//	return minLinkingCosts;
//}

//double cudaFunctionI(const double *xstar, int i,	thrust::device_vector< double >&	d_fixcs	,
//													thrust::device_vector< double >& 	d_costs	,
//													int			d_nfac	,
//													int			d_ncli	){
//
//	double minLinkingCosts = DBL_MAX;
//
//	//for every facility
//	for(int j = 0 ; j < d_nfac; j++ ){
//
//		//cout << "xstar[ypos(" << j <<  ", &instance)]: " <<  xstar[ypos(j, &instance)] << endl;
//
//		//(Opened facility)
//		if( ( xstar[j] ) > 1e-5 ){
//
//			//cout << "M:(" << i+1 << "," << j+1 << "): " << instance.costs[i*instance.n_facilities + j] << endl;
//
//			if(d_costs[i*d_nfac + j] < minLinkingCosts)
//				minLinkingCosts = d_costs[i*d_nfac + j];
//				//cout << "found minimum: " << minLinkingCosts << endl;
//				//cout << "M:(" << i+1 << "," << j+1 << "): " << instance.costs[i*instance.n_facilities + j] << endl;
//
//
//		}
//	}
//
//	return minLinkingCosts;
//}

//// square<T> computes the square of a number f(x) -> x*x
//struct calc_fitness
//{
//    const Instance inst;
//
//    calc_fitness(Instance _instance) : inst(_instance) {}
//
//	__host__ __device__
//        float operator()(const double xstar) const {
//
//    		double z_value = 0.0;
//
//    		double objVal;
//
//			for(int j = 0; j < inst.n_facilities ; j++){
//
//				if(xstar[ypos(j, &inst)] > 0.0 ){
//					objVal += inst.fixed_costs[j];
//				}
//			}
//
//			for(int i = 0 ; i < inst.n_clients ; i++){
//				objVal += cudaFunctionI(inst, &xstar, i);
//			}
//
//			return objVal;
//    	}
//};

/*
void loadInst(Instance* inst){

	//Init and fill with inst.costs
	d_costs = inst->costs;//(inst->costs, inst->costs + (inst->n_clients * inst->n_facilities));

	//Init with costs cardinality
	d_costs_trspnd = inst->costs;

	//Init and fill with fixed_costs
	d_fixed_costs = inst->fixed_costs;

	//Manually "traspose"
	for( int i = 0 ; i < inst->n_clients ; i++){
		for( int j = 0 ; j < inst->n_facilities ; j++ ){
			d_costs_trspnd[j*inst->n_clients + i] = inst->costs[i*inst->n_facilities + j];
		}
	}

	//d_fixed_costs = inst->fixed_costs;

}
*/

//void transpose(double* lin_matrix, double* lin_t_matrix , int n_cols, int n_rows ){
//
//	if(lin_t_matrix == NULL){
//		lin_t_matrix = (double*)malloc(sizeof(double)*n_cols*n_rows );
//	}
//
//	for( int i = 0 ; i < n_rows ; i++){
//		for( int j = 0 ; j < n_cols ; j++){
//			lin_t_matrix[j*n_rows + i] = lin_matrix[i*n_cols + j];
//		}
//	}
//
//}

void loadInstance(Instance& instance,	double& d_fixcs	,
										double& d_costs,
										double& d_t_costs,
										int		d_nfac	,
										int 	d_ncli	){


	cudaMalloc((void**) &d_fixcs	, sizeof(double)*instance.n_facilities);
	CUDAErrorCheck();
	cudaMalloc((void**)	&d_costs	, sizeof(double)*instance.n_facilities*instance.n_clients);
	CUDAErrorCheck();
	cudaMalloc((void**) &d_t_costs	, sizeof(double)*instance.n_facilities*instance.n_clients);
	CUDAErrorCheck();


	cudaMemcpy(	&d_fixcs, instance.fixed_costs, sizeof(double)*d_nfac, cudaMemcpyHostToDevice);
	CUDAErrorCheck();
	cudaMemcpy( &d_costs, instance.costs, d_nfac*d_ncli*sizeof(double), cudaMemcpyHostToDevice );
	CUDAErrorCheck();

	cudaThreadSynchronize();
}

//class calc_fitness{
//
//	const thrust::device_vector< double >& 	d_fixcs;
//	const thrust::device_vector< double >& 	d_t_costs;
//	const int			d_nfac;
//	const int			d_ncli;
//
//public:
//	calc_fitness( 	thrust::device_vector< double >& fixcs_d, thrust::device_vector< double >& costs_t_d,
//					int		nfac_d		, int ncli_d ) :
//		d_fixcs (fixcs_d),
//		d_t_costs (costs_t_d),
//		d_nfac (nfac_d),
//		d_ncli (ncli_d) {}
//
//	__host__ __device__
//	double operator()(const double xstar) const {
//
//		double objVal = 0.0;
//
//		for(int j = 0; j < d_nfac ; j++){
//
//			if(xstar[j] > 0.0 ){
//				objVal += d_fixcs[j];
//			}
//		}
//
//		for(int i = 0 ; i < d_ncli ; i++){
//
//
//			//objVal += cudaFunctionI(xstar, i, d_fixcs, d_costs, d_nfac, d_ncli);
//		}
//
//		return objVal;
//
//	}
//};

//__global__
//void compute_fitness(double* fitness, double* xStar, int n_ind, double* d_fixc	,
//																double* d_costs,
//																int n_fac		,
//																int n_cli ){
//
//	__shared__ double xstar[3000];
//
//	int k = blockDim.x*blockIdx.x + threadIdx.x;
//	fitness[k] = 0.0;
//
//	if(k < n_ind){
//
//		double Wi;
//		double Cij;
//		for(int j = 0 ; j < n_fac ; j++){
//			//fac = xStar[k*n_fac + j];
//			if(xStar[k*n_fac + j] > 1e-5){
//				fitness[k] += d_fixc[j];//*xStar[k*n_fac + j];	//fitness[k] += d_fixc[j];
//			}
//			xstar[j] = xStar[k*n_fac + j];
//		}
//
//		for(int i = 0 ; i < n_cli ;  i++){
//
//			//Mini-fI
//			Wi = DBL_MAX;
//			for( int j = 0 ; j < n_fac ; j++){
//				if(xstar[j] > 1e-5){
//					Cij = d_costs[i*n_fac + j];
//					if(Cij  < Wi ){
//						Wi = Cij;
//					}
//				}
//			}
//
//			fitness[k] += Wi;
//		}
//	}
//}

//__global__
//void compute_fitness(double* fitness, double* xStar, int n_ind, double* d_fixc	,
//																double* d_costs,
//																int n_fac		,
//																int n_cli ){
//
//
//	int k = blockDim.x*blockIdx.x + threadIdx.x;
//	fitness[k] = 0.0;
//
//	if(k < n_ind){
//
//		double Wi;
//		double Cij;
//		for(int j = 0 ; j < n_fac ; j++){
//			//fac = xStar[k*n_fac + j];
//			if(xStar[k*n_fac + j] > 1e-5){
//				fitness[k] += d_fixc[j];//*xStar[k*n_fac + j];	//fitness[k] += d_fixc[j];
//			}
//		}
//
//		for(int i = 0 ; i < n_cli ;  i++){
//
//			//Mini-fI
//			Wi = DBL_MAX;
//			for( int j = 0 ; j < n_fac ; j++){
//				if(xStar[k*n_fac + j] > 1e-5){
//					Cij = d_costs[i*n_fac + j];
//					if(Cij  < Wi ){
//						Wi = Cij;
//					}
//				}
//			}
//
//			fitness[k] += Wi;
//		}
//	}
//}


__global__
void compute_fitness(double* fitness, double* xStar, int n_ind, double* d_fixc, double* d_costs, int n_fac , int n_cli ){

	double Wi;
	double Cij;
	int c = 0;
	int opened[128] = {0};

	int k = blockDim.x*blockIdx.x + threadIdx.x;


	if(k < n_ind){

		fitness[k] = 0.0;

		for(int j = 0 ; j < n_fac ; j++){
			if(xStar[k*n_fac + j] > 1e-5){
				opened[c] = j;
				c++;
				fitness[k] += d_fixc[j];
			}
		}

		for(int i = 0 ; i < n_cli ;  i++){

			//Mini-fI
			c = 0;
			Wi = DBL_MAX;
			for( int j = 0 ; j < n_fac ; j++){
				//if(xStar[k*n_fac + j] > 1e-5){
				if( opened[c] == j ){
					Cij = d_costs[i*n_fac + j];
					if(Cij  < Wi ){
						Wi = Cij;
					}
					c++;
				}
			}
			fitness[k] += Wi;
		}
	}
}

//////////////////////////////////NOT WORKING////////////////////////////////////////
//__global__
//void compute_fitness(double* fitness, double* xStar, int n_ind, double* d_fixc	,
//																double* d_costs,
//																int n_fac		,
//																int n_cli ){
//
//	__shared__ double fitnesses[32];
//
//	int k = blockDim.x*blockIdx.x + threadIdx.x;
//	fitness[k] = 0.0;
//
//	if(k < n_ind){
//
//		double Wi;
//		double Cij;
//
//		for(int j = 0 ; j < n_fac ; j++){
//
//			if(xStar[k*n_fac + j] > 1e-5){
//				//int fi = d_fixc[j];
//				fitnesses[threadIdx.x] += d_fixc[j]*xStar[k*n_fac + j];	//fitness[k] += d_fixc[j];
//			}
//		}
//
//		__syncthreads();
//
//		for(int i = 0 ; i < n_cli ;  i++){
//
//			//Mini-fI
//			Wi = DBL_MAX;
//			for( int j = 0 ; j < n_fac ; j++){
//				if(xStar[k*n_fac + j] > 1e-5){
//					Cij = d_costs[i*n_fac + j];//*xStar[k*n_fac + j];
//					if( Cij < Wi ){
//						Wi = Cij;
//					}
//				}
//				//__syncthreads();
//			}
//
//			__syncthreads();
//			fitnesses[threadIdx.x] += Wi;
//		}
//	}
//	__syncthreads();
//
//	if(k < n_ind){
//		fitness[k] = fitnesses[threadIdx.x];
//	}
//}

void execute_compute_fitness_good(double* fitness, double* xStars, int n_ind,	double* d_fixc, double* d_costs,
																			int n_fac,		int n_cli ,
																			int blocks, 	int threads){

	int tpb = 32;
	//cout <<"n_ind: " << n_ind << endl;
	//double* t_xStars;
//	compute_fitness <<< 16 , n_ind/16 >>> (fitness, xStars, n_ind, d_fixc, d_costs, n_fac, n_cli);

	compute_fitness <<< n_ind/tpb, tpb >>> (fitness, xStars, n_ind, d_fixc, d_costs, n_fac, n_cli);

	cudaThreadSynchronize();
}

void execute_compute_fitness(double* fitness, double* xStars, int n_ind,	double* d_fixc, double* d_costs,
																			int n_fac,		int n_cli ,
																			int blocks, 	int threads){

	int tpb = 128;
	int gridsize = floor( (double) n_ind/tpb ) + 1 ;

	if( VERBOSE ){
		cout << "gridsize: " << gridsize << endl;
		cout <<"n_ind: " << n_ind << endl;
	}
	//	compute_fitness <<< 16 , n_ind/16 >>> (fitness, xStars, n_ind, d_fixc, d_costs, n_fac, n_cli);

	compute_fitness <<< gridsize, tpb >>> (fitness, xStars, n_ind, d_fixc, d_costs, n_fac, n_cli);

	cudaThreadSynchronize();
}

void hybridFitness(vector< Solution >& pop, double* fitness, 	double* d_fixc,
																double* d_costs,
																int n_fac,
																int n_cli,
																double* fixed_costs,
																double* costs){


		double n_ind = pop.size();
		double* d_xstar;
		double* d_fitness;

		cudaMalloc((void**)&d_xstar, (n_fac*n_ind) * sizeof(double));
		CUDAErrorCheck();

		cudaMalloc((void**)&d_fitness, 	sizeof(double)*n_ind);
		CUDAErrorCheck();

		vector< Solution > specials;

		int n_blocks;

		//Memcpy every xStar
		for(int k = 0 ; k < n_ind ; k++){



			int offset = k*n_fac;
			double* ystar = (double*)malloc(sizeof(double)*n_fac);
			//cudaMallocHost((void**)&ystar, sizeof(double)*n_fac);
			memcpy(ystar, pop[k].xStar, sizeof(double)*n_fac);
			//cout << k << "): "; prettyPrintXStar(ystar, n_fac);
			//cout << endl;
			//cudaMemcpy(ystar, pop[k].xStar, sizeof(double)*n_fac, cudaMemcpyHostToHost);

			if(pop[k].non_zero_count < 128 && !( k < floor(n_ind/32)*32 ) )
				cudaMemcpy(d_xstar + offset , ystar, sizeof(double)*n_fac, cudaMemcpyHostToDevice );
			else
				specials.push_back(pop[k]);

			CUDAErrorCheck();

			free(ystar);
			cudaDeviceSynchronize();
		}

		execute_compute_fitness(d_fitness, d_xstar, n_ind, d_fixc, d_costs, n_fac, n_cli, 1, 1);
		CUDAErrorCheck();

		//CPUfitness(specials, fixed_costs, costs);

		cudaMemcpy(fitness, d_fitness, sizeof(double)*n_ind, cudaMemcpyDeviceToHost);
		CUDAErrorCheck();

		//Write fitnesses on sons
		for(int k = 0 ; k < pop.size() ; k++){
			pop[k].z_opt = fitness[k];
			//prettyPrintXStar(pop[k].xStar, n_fac); cout << ": " << pop[k].z_opt << endl;
		}

		cudaFree(d_xstar);
		cudaFree(d_fitness);

}

void cudaFitness(vector< Solution >& pop, double* fitness, 	double* d_fixc,
															double* d_costs,
															int n_fac,
															int n_cli){

	int n_ind = pop.size();
	double* d_xstar;
	double* d_fitness;

	cudaMalloc((void**)&d_xstar, (n_fac*n_ind) * sizeof(double));
	CUDAErrorCheck();

	cudaMalloc((void**)&d_fitness, 	sizeof(double)*n_ind);
	CUDAErrorCheck();

	//vector< Solution > specials;

	//Memcpy every xStar
	for(int k = 0 ; k < n_ind ; k++){

		int offset = k*n_fac;
		double* ystar = (double*)malloc(sizeof(double)*n_fac);
		//cudaMallocHost((void**)&ystar, sizeof(double)*n_fac);
		memcpy(ystar, pop[k].xStar, sizeof(double)*n_fac);
		//cout << k << "): "; prettyPrintXStar(ystar, n_fac);
		//cout << endl;
		//cudaMemcpy(ystar, pop[k].xStar, sizeof(double)*n_fac, cudaMemcpyHostToHost);
		//if(pop[k].non_zero_count < 128)
		cudaMemcpy(d_xstar + offset , ystar, sizeof(double)*n_fac, cudaMemcpyHostToDevice );
		//else
		//	specials.push_back(pop[k]);

		CUDAErrorCheck();

		free(ystar);
		cudaDeviceSynchronize();
	}

	execute_compute_fitness(d_fitness, d_xstar, n_ind, d_fixc, d_costs, n_fac, n_cli, 1, 1);
	CUDAErrorCheck();

	//CPUfitness(specials, ins)

	cudaMemcpy(fitness, d_fitness, sizeof(double)*n_ind, cudaMemcpyDeviceToHost);
	CUDAErrorCheck();

	//Write fitnesses on sons
	for(int k = 0 ; k < n_ind ; k++){
		pop[k].z_opt = fitness[k];
		//prettyPrintXStar(pop[k].xStar, n_fac); cout << ": " << pop[k].z_opt << endl;
	}

	cudaFree(d_xstar);
	cudaFree(d_fitness);

}


void cudaFitness2(vector< Solution >& pop, double* fitness, 	double* d_fixc,
																double* d_costs,
																int n_fac,
																int n_cli){
	int n_ind = pop.size();
//	cout << "sons to calc: " << n_ind << endl;
//	if( true ){	//Split if too big ( or else Cuda Kernel goes timeout )
//	double* fitness = (double*)calloc(n_ind, sizeof(double));

	cudaFitness(pop, fitness, d_fixc, d_costs, n_fac, n_cli);//, d_fixc, d_costs, instance.n_facilities, instance.n_clients);

//			for(int k = 0 ; k < pop.size() ; k++ ){
//
//				cout << k << ")" << "GPU fitness: " << fitnesses[k] << " CPU fitness: " << pop[k].z_opt << endl;
//
//			}
	//fitness = fitnesses;

	//free(fitnesses);

//	} else {

//		vector< Solution > sons_c;
//
//		int CHUNCK_SIZE = 128;
//		//int N_CHUNCKS =	n_ind/CHUNCK_SIZE;
//
//		//Create Chuncks
//		for(int c = 0 ; c < ceil(n_ind/CHUNCK_SIZE) ; c++ ){
//
//			vector< Solution > chunck (CHUNCK_SIZE);
//			double* fitnesses = (double*)malloc(sizeof(double)*CHUNCK_SIZE);
//			copy(pop.begin() + CHUNCK_SIZE*c, pop.begin() + CHUNCK_SIZE*(c+1), chunck.begin() );
//			cudaFitness(chunck, fitnesses, d_fixc, d_costs, n_fac, n_cli);
//			free(fitnesses);
//			sons_c.insert(sons_c.end(), chunck.begin(), chunck.end());
//
//		}
//		pop = sons_c;

//	}

	//cout << "after: " << pop.size();

}



//void cudaSort(Instance* inst, double* ystar, YStar** y_ord){
//
//	//Initializing and loading device vectors
//	thrust::device_vector<int> 		d_indexes( inst->n_facilities);//	((int) inst->n_facilities);
//	thrust::device_vector<double> 	d_lpValues(ystar, ystar + inst->n_facilities) ;//	((int) inst->n_facilities);
//
////	thrust::copy(ystar, ystar + ( (uint)inst->n_facilities*sizeof(double) ), d_lpValues.begin());
//
//	//Filling index vector from 0 to (#Facilities-1)
//	thrust::sequence(d_indexes.begin(), d_indexes.end());
//
///*
//	for (int k = 0; k < inst->n_facilities; k++) {
//		d_lpValues[k] = ystar[k];
//	}
//
//	//thrust::copy(ystar, ystar + inst->n_facilities, d_lpValues);
//	//cudaMemcpy(&d_lpValues, ystar, sizeof(double) * inst->n_facilities, cudaMemcpyHostToDevice);
//
//	//d_lpValues = ystar;
//
//	//printf("sfkjghfelkgjfhgkjf\n");
//*/
//
//	//SORT!
//	thrust::stable_sort_by_key(d_lpValues.begin(), d_lpValues.end() , d_indexes.begin(), thrust::greater<double>());
//
//	//Device2Host copy
//	thrust::host_vector<int> index = d_indexes;
//	thrust::host_vector<double> lp_Values = d_lpValues;
//
//	//Re-Create array of structures
//	for (int k = 0; k < inst->n_facilities; k++) {
//		y_ord[k] = (YStar *) malloc(sizeof(YStar));
//		y_ord[k]->index = index[k];
//		y_ord[k]->lpValue = lp_Values[k];
//	}
//}

//void cudaFitness(Instance& inst, vector< Solution >& population, double* fitness){
//
//	thrust::device_vector<double*> 	d_xstars(population.size());
//	thrust::device_vector<double >	d_f_values(population.size());
//
//	//Upload xstars
//	for (int k = 0; k < population.size(); k++) {
//		d_xstars[k] = population[k].xStar;
//	}
//	thrust::transform(d_xstars.begin(), d_xstars.end(), d_f_values.begin() , calc_fitness<double*>(d_xstars, d_f_values ));
//
//}

//void cudaSort(Instance* inst, double* ystar, Map** y_ord){
//
//	//YStar **y_ord = (YStar **) malloc(inst->n_facilities * sizeof(YStar*));
//
//	//Initializing and loading device vectors
//	thrust::device_vector<int> 		d_indexes(inst->n_facilities);//	((int) inst->n_facilities);
//	thrust::device_vector<double> 	d_lpValues(ystar, ystar + inst->n_facilities) ;//	((int) inst->n_facilities);
//
//	//Filling index vector from 0 to (#Facilities-1)
//	thrust::sequence(d_indexes.begin(), d_indexes.end());
//
//	//SORT!
//	thrust::stable_sort_by_key(d_lpValues.begin(), d_lpValues.end() , d_indexes.begin(), thrust::greater<double>());
//
//	//Device2Host copy
//	thrust::host_vector<int> index = d_indexes;
//	thrust::host_vector<double> lp_Values = d_lpValues;
//
//	//Re-Create array of structures
//	for (int k = 0; k < inst->n_facilities; k++) {
//		y_ord[k] = (Map *) malloc(sizeof(Map));
//		y_ord[k]->index = index[k];
//		y_ord[k]->value = lp_Values[k];
//	}
//}
