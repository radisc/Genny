#include <cuda.h>
#include <cuda_runtime_api.h>


void CUDAErrorCheck();

//YStar** cudaSort(Instance* inst, double* ystar);
void cudaSort(Instance* inst, double* ystar, YStar** y_ord);
void cudaSort(Instance* inst, double* ystar, Map** y_ord);

double cudaFunctionI(Instance& instance, double *xstar, int i);

//double cudaFunctionI(const double *xstar, int i,	const thrust::device_vector< double >&	d_fixcs	,
//													const thrust::device_vector< double >& 	d_costs	,
//													const int			d_nfac	,
//													const int			d_ncli	);

void loadInstance(Instance& instance,	double& d_fixcs	,
										double& d_costs,
										double& d_t_costs,
										int		d_nfac	,
										int 	d_ncli	);

void loadInstance(Instance& instance,	double* d_fixcs	,
										double* d_costs,
										double* d_t_costs,
										int		d_nfac	,
										int 	d_ncli	);

extern void execute_compute_fitness(double* fitness, double* xStar, int n_ind,	double* d_fixc, double* d_costs,
																				int n_fac,		int n_cli ,
																				int blocks, int threads);

//thrust::device_ptr<Instance> loadInstance(Instance& instance);

void cudaFitness2(vector< Solution >& pop, double* fitness, 	double* d_fixc,
																double* d_costs,
																int n_fac,
																int n_cli);

void cudaFitness(vector< Solution >& population, double* fitness);

void cudaFitness(vector< Solution >& pop, double* fitness, 	double* d_fixc,
															double* d_costs,
															int n_fac,
															int n_cli);
