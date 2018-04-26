/*
 * genetics.cpp
 *
 *  Created on: 19/mag/2015
 *      Author: nicola
 */

#include "vector"
#include "math.h"
#include "iostream"
#include "limits"
#include "algorithm"
#include "genetics.h"
#include "ufl.h"

#include <sys/time.h>

/*
 * Generate random individuals ( with parent_nzcnt opened facilities for the first iteration
 *  or if the population doesn't have minimum number at beginning of iteration
 */
Solution generate_ramdom_solution(Instance& inst, int parent_nzcnt ,double sparseness, bool integer){

    int ncols = inst.n_facilities + inst.n_clients;

    Solution random;
    random.xStar = (double*)malloc(sizeof(double)*(ncols));
    memset(random.xStar, 0.0f, sizeof(double)*(ncols));
    random.z_opt = 0.0f;
    random.non_zero_count = 0;

    int variability = 0;

    //int nzcnt_max =  (rand()%(int(inst.n_facilities)))/sparseness ;//*sparseness; //+ ((rand()%(variability/2)) - variability);

    int nzcnt; //= (parent_nzcnt == 0 ) ? (rand()%inst.n_facilities + 1) : parent_nzcnt; //or parent_nzcnt

    if( parent_nzcnt < 2 ){
    	nzcnt = rand()%3 + 1;
    } else
    	nzcnt = parent_nzcnt;

    //cout << "nzcnt: " << nzcnt << endl;

    //Extract random genes

    while(nzcnt>0){

        int r_idx = rand()%inst.n_facilities;

        if(random.xStar[ypos(r_idx, &inst)] != 1.0f){

            random.xStar[ypos(r_idx, &inst)] = 1.0f;
            random.z_opt += inst.fixed_costs[r_idx];
            random.non_zero_count++;
            nzcnt--;
            //printf("Newcomers cameth\n");
        }
    }

/*
	for(int j = 0 ; j < inst.n_facilities ; j++){

        if(rand() < RAND_MAX*sparseness){
            random.xStar[ypos(j, &inst)] = 1.0f;
            random.z_opt += inst.fixed_costs[j];
            nzcnt++;
            //cout <<  inst.fixed_costs[j] << endl;
        } else {
            random.xStar[j] = 0.0f;
        }
    }
*/

//	Open at least one
//  if(nzcnt == 0){
//      //cout << "WDHFGBCWGINORWGTITBRVETVO" << endl;
//      int idx = rand()%inst.n_facilities;
//      random.xStar[ypos(idx, &inst)] = 1.0;
//      random.z_opt += inst.fixed_costs[idx];
//  }
//	cout << "New individual zopt beef: " << random.z_opt << endl;

	#pragma omp parallel for
    for(int i = 0 ; i < inst.n_clients ; i++){

       	double Wi = intFunctionI(inst, random.xStar, i);
       	//cout << "W(" << i << "):" <<  Wi << endl;
       	random.xStar[wpos(i, &inst)] = Wi;
       	random.z_opt += Wi;
    }

    //cout << "New individual zopt afft: " << random.z_opt << endl;

    return random;
}

/*
 *Harshness ratio of surviving individuals
 */
void selection(vector<Solution>& population, int size, int elite){




	//population.resize(floor(population.size()*harshness));
	while(population.size() > size ){

		//Deadpool
		int id = elite + rand()%(population.size() - elite);

		free(population[id].xStar);
		population.erase(population.begin() + id);
		//delete (population[population.size() - 1]);

	}

}

int mutate(int index, int gene){

//    if(index == 1925){
//        cout << "E' uscito!" << endl;
//        //exit(0);
//    }

    if (gene == 0.0){
        //cout << "mutation of " << index << " to " << 1 <<  endl;
        return 1.0;
    } else{
        //cout << "mutation of " << index << " to " << 0 <<  endl;
        return 0.0;
    }

}

int mutate(int index, int gene, Solution& gene_pool){

	gene_pool.xStar[index]++;

    if (gene == 0.0){
        //cout << "mutation of " << index << " to " << 1 <<  endl;
        return 1.0;
    } else{
        //cout << "mutation of " << index << " to " << 0 <<  endl;
        return 0.0;
    }

}

Solution reverse(Instance& inst, Solution& solution, bool compute_fitness){

	Solution reversed;
	reversed.non_zero_count = 0;
	reversed.xStar = (double*)malloc( sizeof(double)*(inst.n_facilities+inst.n_clients) );
	memset(reversed.xStar, 0.0f, sizeof(double)*(inst.n_facilities+inst.n_clients));

	//cout << "fkunhtregvypwwctrugybvwliutvmtvb"<< endl;

	for(int j = 0 ; j < inst.n_facilities ; j++){

		if(solution.xStar[j] == 0.0)
			reversed.xStar[j] = 1.0;
		else
			reversed.xStar[j] = 0.0;

		if(reversed.xStar[j] > 1e-5){
			reversed.non_zero_count++;
		}
	}

	//cout << "fdkuthsiugvntoiuthvno"<< endl;

	if(compute_fitness)
		reversed.z_opt = calcZ(reversed.xStar, inst, inst.n_facilities, inst.n_clients, inst.fixed_costs, inst.costs);

	return reversed;
}

Solution* kCutsCrossover(Instance& inst, Solution a, Solution b, double mutation, Solution& gene_pool, int n_cut_points, bool compute_fitness){

	int ncols = inst.n_facilities + inst.n_clients;

	Solution* children = new Solution[2];
	children[0].xStar = (double*) calloc(ncols, sizeof(double));        //new double[inst.n_facilities + inst.n_clients];
	children[1].xStar = (double*) calloc(ncols, sizeof(double));
	children[0].z_opt = 0.0;
	children[1].z_opt = 0.0;
	children[0].non_zero_count = 0;
	children[1].non_zero_count = 0;

	//cout << "n_cut_points: " << n_cut_points << endl;

	int* cut_points = new int[n_cut_points + 1];

	if(false){
		prettyPrint(a, inst);
		cout << " and " << endl;
		prettyPrint(b, inst);
	}


	for(int k = 0 ; k < n_cut_points ; k++ ){
		cut_points[k] = rand()%inst.n_facilities;
	}
	cut_points[n_cut_points] = inst.n_facilities;

	sort(cut_points, cut_points + n_cut_points + 1);

//	int second_cut_point = rand()%inst.n_facilities;          //Matheuristics ( ??? )
//	int mutated_gene_a = rand()%inst.n_facilities;              //extract one/two from pop0's and one/two from others
//	int mutated_gene_b = rand()%inst.n_facilities;
//	int gene_idx;
//    for(int j = 0 ; j < cut_point; j++){
//
//        //if(j == mutated_gene){                                                            //if this on...
//        if(rand() < RAND_MAX/inst.n_facilities){//mutation){                             //...this off
//            children[0].xStar[ypos(j, &inst)] = memcpy(a + cutpoint);                           //mutate(j, b.xStar[ypos(j, &inst)]);
//            children[1].xStar[ypos(j, &inst)] = mutate(j, a.xStar[ypos(j, &inst)]);
//        } else {
//            children[0].xStar[ypos(j, &inst)] = b.xStar[ypos(j, &inst)];
//            children[1].xStar[ypos(j, &inst)] = a.xStar[ypos(j, &inst)];
//        }
//    }

	//First chunk
	memcpy(children[0].xStar, a.xStar, sizeof(double)*(cut_points[0]));
	memcpy(children[1].xStar, b.xStar, sizeof(double)*(cut_points[0]));

	for(int k = 0 ; k < n_cut_points ; k++){

		if(k%2 == 0){
			//crossover odd chunk
			memcpy(&children[0].xStar[cut_points[k]], &a.xStar[cut_points[k]] , sizeof(double) * (cut_points[k+1] - cut_points[k]));
			memcpy(&children[1].xStar[cut_points[k]], &b.xStar[cut_points[k]] , sizeof(double) * (cut_points[k+1] - cut_points[k]));

		} else {
			//Crossover even chunk
			memcpy(&children[0].xStar[cut_points[k]], &b.xStar[cut_points[k]] , sizeof(double) * (cut_points[k+1] - cut_points[k]));
			memcpy(&children[1].xStar[cut_points[k]], &a.xStar[cut_points[k]] , sizeof(double) * (cut_points[k+1] - cut_points[k]));

		}
	}


	//Last chunk

	if(false){
		cout << " made: " << endl;
		prettyPrint(children[0], inst);
		cout << " and " << endl;
		prettyPrint(children[1], inst);
		cout << "Cut points("<< n_cut_points <<"): ";
		for(int k = 0 ; k < (n_cut_points + 1) ; k++ ){
			cout << cut_points[k] << " " ;
		}
		cout << endl;
		cout << "*********" << endl;
	}

	//Apply mutation with 85% probability
	//if( rand() < RAND_MAX*0.70 ){
		for(int j = 0; j < inst.n_facilities ; j++){

			if(rand() < RAND_MAX/inst.n_facilities*mutation)	//2, 4			3: Ok per ga250a-5
				children[0].xStar[ypos(j, &inst)] = mutate(j, b.xStar[ypos(j, &inst)], gene_pool);//, inst, children[0]);

			if(rand() < RAND_MAX/inst.n_facilities*mutation)	//2, 4
				children[1].xStar[ypos(j, &inst)] = mutate(j, a.xStar[ypos(j, &inst)], gene_pool);//, inst, children[1]);

		}
	//}

	if(compute_fitness){

		//struct timeval start, end;

		//gettimeofday(&start, NULL);
		//Fitness computing
		for(int j = 0 ; j < inst.n_facilities ; j++){
			if(children[0].xStar[ypos(j, &inst)] == 1.0){
				children[0].z_opt += inst.fixed_costs[j];
				children[0].non_zero_count++;
			}
			if(children[1].xStar[ypos(j, &inst)] == 1.0){
				children[1].z_opt += inst.fixed_costs[j];
				children[1].non_zero_count++;
			}
		}

		for(int i = 0 ; i < inst.n_clients ; i++){

			double Wi_0 = intFunctionI(inst, children[0].xStar, i);
			double Wi_1 = intFunctionI(inst, children[1].xStar, i);
			//cout << "W(" << i << "):" <<  Wi_0;
			//cout << "W(" << i << "):" <<  Wi_1;
			children[0].xStar[wpos(i, &inst)] = Wi_0;
			children[0].z_opt += Wi_0;
			children[1].xStar[wpos(i, &inst)] = Wi_1;
			children[1].z_opt += Wi_1;
		}

		//gettimeofday(&end, NULL);
		//duration += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	} else {

		//Non zero count
		for(int j = 0 ; j < inst.n_facilities ; j++){
			if(children[0].xStar[ypos(j, &inst)] == 1.0){
				//children[0].z_opt += inst.fixed_costs[j];
				children[0].non_zero_count++;
			}
			if(children[1].xStar[ypos(j, &inst)] == 1.0){
				//children[1].z_opt += inst.fixed_costs[j];
				children[1].non_zero_count++;
			}
		}

//		for(int i = 0 ; i < inst.n_clients ; i++){
//
//			double Wi_0 = intFunctionI(inst, children[0].xStar, i);
//			double Wi_1 = intFunctionI(inst, children[1].xStar, i);
//			cout << "W(" << i << "):" <<  Wi_0;
//			cout << "W(" << i << "):" <<  Wi_1;
//			children[0].xStar[wpos(i, &inst)] = Wi_0;
//			children[0].z_opt += Wi_0;
//			children[1].xStar[wpos(i, &inst)] = Wi_1;
//			children[1].z_opt += Wi_1;
//		}

	}
	return children;
}

/*
 * k cuts crossover
 */
//Solution* kCutsCrossover(Instance& inst, Solution a, Solution b, double mutation, Solution& gene_pool, int n_cut_points, bool compute_fitness){
//
//	int ncols = inst.n_facilities + inst.n_clients;
//
//	Solution* children = new Solution[2];
//	children[0].xStar = (double*) calloc(ncols, sizeof(double));        //new double[inst.n_facilities + inst.n_clients];
//	children[1].xStar = (double*) calloc(ncols, sizeof(double));
//	children[0].z_opt = 0.0;
//	children[1].z_opt = 0.0;
//	children[0].non_zero_count = 0;
//	children[1].non_zero_count = 0;
//
//	//cout << "n_cut_points: " << n_cut_points << endl;
//
//	int* cut_points = new int[n_cut_points + 1];
//
//	if(false){
//		prettyPrint(a, inst);
//		cout << " and " << endl;
//		prettyPrint(b, inst);
//	}
//
//
//	for(int k = 0 ; k < n_cut_points ; k++ ){
//		cut_points[k] = rand()%inst.n_facilities;
//	}
//	cut_points[n_cut_points] = inst.n_facilities;
//
//	sort(cut_points, cut_points + n_cut_points + 1);
//
////	int second_cut_point = rand()%inst.n_facilities;          //Matheuristics ( ??? )
////	int mutated_gene_a = rand()%inst.n_facilities;              //extract one/two from pop0's and one/two from others
////	int mutated_gene_b = rand()%inst.n_facilities;
////	int gene_idx;
//
////    for(int j = 0 ; j < cut_point; j++){
////
////        //if(j == mutated_gene){                                                            //if this on...
////        if(rand() < RAND_MAX/inst.n_facilities){//mutation){                             //...this off
////            children[0].xStar[ypos(j, &inst)] = memcpy(a + cutpoint);                           //mutate(j, b.xStar[ypos(j, &inst)]);
////            children[1].xStar[ypos(j, &inst)] = mutate(j, a.xStar[ypos(j, &inst)]);
////        } else {
////            children[0].xStar[ypos(j, &inst)] = b.xStar[ypos(j, &inst)];
////            children[1].xStar[ypos(j, &inst)] = a.xStar[ypos(j, &inst)];
////        }
////    }
//
//
//	//First chunk
//	memcpy(children[0].xStar, a.xStar, sizeof(double)*(cut_points[0]));
//	memcpy(children[1].xStar, b.xStar, sizeof(double)*(cut_points[0]));
//
//	for(int k = 0 ; k < n_cut_points ; k++){
//
//		if(k%2 == 0){
//			//crossover odd chunk
//			memcpy(&children[0].xStar[cut_points[k]], &a.xStar[cut_points[k]] , sizeof(double) * (cut_points[k+1] - cut_points[k]));
//			memcpy(&children[1].xStar[cut_points[k]], &b.xStar[cut_points[k]] , sizeof(double) * (cut_points[k+1] - cut_points[k]));
//
//		} else {
//			//Crossover even chunk
//			memcpy(&children[0].xStar[cut_points[k]], &b.xStar[cut_points[k]] , sizeof(double) * (cut_points[k+1] - cut_points[k]));
//			memcpy(&children[1].xStar[cut_points[k]], &a.xStar[cut_points[k]] , sizeof(double) * (cut_points[k+1] - cut_points[k]));
//
//		}
//	}
//
//
//	//Last chunk
//
//	if(false){
//		cout << " made: " << endl;
//		prettyPrint(children[0], inst);
//		cout << " and " << endl;
//		prettyPrint(children[1], inst);
//		cout << "Cut points("<< n_cut_points <<"): ";
//		for(int k = 0 ; k < (n_cut_points + 1) ; k++ ){
//			cout << cut_points[k] << " " ;
//		}
//		cout << endl;
//		cout << "*********" << endl;
//	}
//
//	//Apply mutation with 85% probability
//	//if( rand() < RAND_MAX*0.85 ){
//		for(int j = 0; j < inst.n_facilities ; j++){
//
//			if(rand() < RAND_MAX/inst.n_facilities*mutation)	//2, 4			3: Ok per ga250a-5
//				children[0].xStar[ypos(j, &inst)] = mutate(j, b.xStar[ypos(j, &inst)], gene_pool);//, inst, children[0]);
//
//			if(rand() < RAND_MAX/inst.n_facilities*mutation)	//2, 4
//				children[1].xStar[ypos(j, &inst)] = mutate(j, a.xStar[ypos(j, &inst)], gene_pool);//, inst, children[1]);
//
//		}
//	//}
//
//	if(compute_fitness){
//
//		//Fitness computing
//		for(int j = 0 ; j < inst.n_facilities ; j++){
//			if(children[0].xStar[ypos(j, &inst)] == 1.0){
//				children[0].z_opt += inst.fixed_costs[j];
//				children[0].non_zero_count++;
//			}
//			if(children[1].xStar[ypos(j, &inst)] == 1.0){
//				children[1].z_opt += inst.fixed_costs[j];
//				children[1].non_zero_count++;
//			}
//		}
//
//		for(int i = 0 ; i < inst.n_clients ; i++){
//
//			double Wi_0 = intFunctionI(inst, children[0].xStar, i);
//			double Wi_1 = intFunctionI(inst, children[1].xStar, i);
//			//cout << "W(" << i << "):" <<  Wi_0;
//			//cout << "W(" << i << "):" <<  Wi_1;
//			children[0].xStar[wpos(i, &inst)] = Wi_0;
//			children[0].z_opt += Wi_0;
//			children[1].xStar[wpos(i, &inst)] = Wi_1;
//			children[1].z_opt += Wi_1;
//		}
//
//	} else {
//
//		//Non zero count
//		for(int j = 0 ; j < inst.n_facilities ; j++){
//			if(children[0].xStar[ypos(j, &inst)] == 1.0){
//				//children[0].z_opt += inst.fixed_costs[j];
//				children[0].non_zero_count++;
//			}
//			if(children[1].xStar[ypos(j, &inst)] == 1.0){
//				//children[1].z_opt += inst.fixed_costs[j];
//				children[1].non_zero_count++;
//			}
//		}
//
////		for(int i = 0 ; i < inst.n_clients ; i++){
////
////			double Wi_0 = intFunctionI(inst, children[0].xStar, i);
////			double Wi_1 = intFunctionI(inst, children[1].xStar, i);
////			cout << "W(" << i << "):" <<  Wi_0;
////			cout << "W(" << i << "):" <<  Wi_1;
////			children[0].xStar[wpos(i, &inst)] = Wi_0;
////			children[0].z_opt += Wi_0;
////			children[1].xStar[wpos(i, &inst)] = Wi_1;
////			children[1].z_opt += Wi_1;
////		}
//
//	}
//	return children;
//}


/*
 * Dual cut crossover
 */
Solution* crossover(Instance& inst, Solution a, Solution b, double mutation, Solution& gene_pool, int n_cut_points, bool compute_fitness){

	///////////////////////////////////////
	return kCutsCrossover(inst, a, b, mutation, gene_pool, n_cut_points, compute_fitness);
//	return simpleCrossover(inst, a, b, mutation, gene_pool);
	///////////////////////////////////////
	int ncols = inst.n_facilities + inst.n_clients;

	Solution* children = new Solution[2];
	children[0].xStar = (double*) calloc(ncols, sizeof(double));        //new double[inst.n_facilities + inst.n_clients];
	children[1].xStar = (double*) calloc(ncols, sizeof(double));
	children[0].z_opt = 0.0;
	children[1].z_opt = 0.0;
	children[0].non_zero_count = 0;
	children[1].non_zero_count = 0;
	int first_cut_point = 0;
	int second_cut_point = 0;

	//if(rand() < RAND_MAX)
	//if(rand()*85 < RAND_MAX/100)
	first_cut_point = rand()%inst.n_facilities;

	second_cut_point = rand()%inst.n_facilities;

	//switch in case
	if(second_cut_point < first_cut_point){
		int temp;
		temp = first_cut_point;
		first_cut_point = second_cut_point;
		second_cut_point = temp;
	}


	//int second_cut_point = rand()%inst.n_facilities;          //Matheuristics
	int mutated_gene_a = rand()%inst.n_facilities;              //extract one/two from pop0's and one/two from others
	int mutated_gene_b = rand()%inst.n_facilities;
																//for(  0 to mutations ) *pick for mutation*
	int gene_idx;

//    for(int j = 0 ; j < cut_point; j++){
//
//        //if(j == mutated_gene){                                                            //if this on...
//        if(rand() < RAND_MAX/inst.n_facilities){//mutation){                             //...this off
//            children[0].xStar[ypos(j, &inst)] = memcpy(a + cutpoint);                           //mutate(j, b.xStar[ypos(j, &inst)]);
//            children[1].xStar[ypos(j, &inst)] = mutate(j, a.xStar[ypos(j, &inst)]);
//        } else {
//            children[0].xStar[ypos(j, &inst)] = b.xStar[ypos(j, &inst)];
//            children[1].xStar[ypos(j, &inst)] = a.xStar[ypos(j, &inst)];
//        }
//    }

	//cout << "Memcopying\n";

	//Crossover first part
	memcpy(children[0].xStar, a.xStar, sizeof(double)*(first_cut_point));
	memcpy(children[1].xStar, b.xStar, sizeof(double)*(first_cut_point));
	//copy(a.xStar, a.xStar + cut_point, children[0].xStar);
	//copy(b.xStar, b.xStar + cut_point, children[1].xStar);

   // cout << "First half done\n";

	//Crossover second part
	memcpy(&children[0].xStar[first_cut_point], &b.xStar[first_cut_point] , sizeof(double) * (second_cut_point - first_cut_point));
	memcpy(&children[1].xStar[first_cut_point], &a.xStar[first_cut_point] , sizeof(double) * (second_cut_point - first_cut_point));

	//Crossover third part

	memcpy(&children[0].xStar[second_cut_point], &a.xStar[second_cut_point] , sizeof(double) * (inst.n_facilities - second_cut_point));
	memcpy(&children[1].xStar[second_cut_point], &b.xStar[second_cut_point] , sizeof(double) * (inst.n_facilities - second_cut_point));


//	prettyPrint(a, inst);
//	cout << " and " << endl;
//	prettyPrint(b, inst);
//	cout << " made: " << endl;
//	prettyPrint(children[0], inst);
//	cout << " and " << endl;
//	prettyPrint(children[1], inst);

	//Apply mutation with 85% probability

//	if( rand() < RAND_MAX*(15/100) ){
//		for(int j = 0; j < inst.n_facilities ; j++){
//
//			if(rand() < RAND_MAX/inst.n_facilities)
//				children[0].xStar[ypos(j, &inst)] = mutate(j, b.xStar[ypos(j, &inst)], gene_pool);//, inst, children[0]);
//		}
//	}
//
//	if( rand() < RAND_MAX*(15/100) ){
//		for(int j = 0; j < inst.n_facilities ; j++){
//			if(rand() < RAND_MAX/inst.n_facilities)
//				children[1].xStar[ypos(j, &inst)] = mutate(j, a.xStar[ypos(j, &inst)], gene_pool);//, inst, children[1]);
//		}
//	}

	if( rand() < RAND_MAX*0.90 ){
		for(int j = 0; j < inst.n_facilities ; j++){

			if(children[0].xStar[ypos(j, &inst)] == 0.0){
				if(rand() < RAND_MAX/(inst.n_facilities)*2.5)
					children[0].xStar[ypos(j, &inst)] = mutate(j, b.xStar[ypos(j, &inst)], gene_pool);//, inst, children[0]);
			} else {
				if(rand() < RAND_MAX/(inst.n_facilities)*0.4)
					children[0].xStar[ypos(j, &inst)] = mutate(j, b.xStar[ypos(j, &inst)], gene_pool);//, inst, children[0]);
			}

			if(children[1].xStar[ypos(j, &inst)] == 0.0){
				if(rand() < RAND_MAX/(inst.n_facilities)*2.5)
					children[1].xStar[ypos(j, &inst)] = mutate(j, a.xStar[ypos(j, &inst)], gene_pool);//, inst, children[1]);
			} else {
				if(rand() < RAND_MAX/(a.non_zero_count)*0.4){
					children[1].xStar[ypos(j, &inst)] = mutate(j, a.xStar[ypos(j, &inst)], gene_pool);//, inst, children[1]);
				}

			}
		}
	}

//	//if( true ){//rand() < RAND_MAX*0.95 ){
//		for(int j = 0; j < inst.n_facilities ; j++){
//
//			if(rand()/2 < RAND_MAX/inst.n_facilities)
//				children[0].xStar[ypos(j, &inst)] = mutate(j, b.xStar[ypos(j, &inst)], gene_pool);//, inst, children[0]);
//
//			if(rand()/2 < RAND_MAX/inst.n_facilities)
//				children[1].xStar[ypos(j, &inst)] = mutate(j, a.xStar[ypos(j, &inst)], gene_pool);//, inst, children[1]);
//		}
//	//}

	//}

	if(compute_fitness){			//We are NOT using Cuda
		//Fitness computing
		for(int j = 0 ; j < inst.n_facilities ; j++){
			if(children[0].xStar[ypos(j, &inst)] == 1.0){
				children[0].z_opt += inst.fixed_costs[j];
				children[0].non_zero_count++;
			}
			if(children[1].xStar[ypos(j, &inst)] == 1.0){
				children[1].z_opt += inst.fixed_costs[j];
				children[1].non_zero_count++;
			}
		}

		for(int i = 0 ; i < inst.n_clients ; i++){

			double Wi_0 = intFunctionI(inst, children[0].xStar, i);
			double Wi_1 = intFunctionI(inst, children[1].xStar, i);
			cout << "W(" << i << "):" <<  Wi_0;
			cout << "W(" << i << "):" <<  Wi_1;
			children[0].xStar[wpos(i, &inst)] = Wi_0;
			children[0].z_opt += Wi_0;
			children[1].xStar[wpos(i, &inst)] = Wi_1;
			children[1].z_opt += Wi_1;
		}
	} else {
		for(int j = 0 ; j < inst.n_facilities ; j++){
			if(children[0].xStar[ypos(j, &inst)] == 1.0){
				//children[0].z_opt += inst.fixed_costs[j];
				children[0].non_zero_count++;
			}
			if(children[1].xStar[ypos(j, &inst)] == 1.0){
				//children[1].z_opt += inst.fixed_costs[j];
				children[1].non_zero_count++;
			}
		}
	}

	return children;
}

Solution* simpleCrossover(Instance& inst, Solution a, Solution b, double mutation, Solution& gene_pool){

	int ncols = inst.n_facilities + inst.n_clients;

	Solution* children = new Solution[2];
	children[0].xStar = (double*) calloc(ncols, sizeof(double));        //new double[inst.n_facilities + inst.n_clients];
	children[1].xStar = (double*) calloc(ncols, sizeof(double));
	children[0].z_opt = 0.0;
	children[1].z_opt = 0.0;
	children[0].non_zero_count = 0;
	children[1].non_zero_count = 0;
	int cut_point = 0;

	//if(rand() < RAND_MAX)
	//if(rand()*85 < RAND_MAX/100)
		cut_point = rand()%inst.n_facilities;

	//int second_cut_point = rand()%inst.n_facilities;          //Matheuristics
	int mutated_gene_a = rand()%inst.n_facilities;              //extract one/two from pop0's and one/two from others
	int mutated_gene_b = rand()%inst.n_facilities;
																//for(  0 to mutations ) *pick for mutation*
	int gene_idx;

//    for(int j = 0 ; j < cut_point; j++){
//
//        //if(j == mutated_gene){                                                            //if this on...
//        if(rand() < RAND_MAX/inst.n_facilities){//mutation){                             //...this off
//            children[0].xStar[ypos(j, &inst)] = memcpy(a + cutpoint);                           //mutate(j, b.xStar[ypos(j, &inst)]);
//            children[1].xStar[ypos(j, &inst)] = mutate(j, a.xStar[ypos(j, &inst)]);
//        } else {
//            children[0].xStar[ypos(j, &inst)] = b.xStar[ypos(j, &inst)];
//            children[1].xStar[ypos(j, &inst)] = a.xStar[ypos(j, &inst)];
//        }
//    }

	//cout << "Memcopying\n";

	//Crossover first part
	memcpy(children[0].xStar, a.xStar, sizeof(double)*(cut_point));
	memcpy(children[1].xStar, b.xStar, sizeof(double)*(cut_point));
	//copy(a.xStar, a.xStar + cut_point, children[0].xStar);
	//copy(b.xStar, b.xStar + cut_point, children[1].xStar);

   // cout << "First half done\n";

	//Crossover second part
	memcpy(&children[0].xStar[cut_point], &b.xStar[cut_point] , sizeof(double) * (inst.n_facilities - cut_point));
	memcpy(&children[1].xStar[cut_point], &a.xStar[cut_point] , sizeof(double) * (inst.n_facilities - cut_point));

	//copy(a.xStar + cut_point , a.xStar + inst.n_facilities, children[1].xStar + cut_point );
	//copy(b.xStar + cut_point , b.xStar + inst.n_facilities, children[0].xStar + cut_point );

//	  cout << "Second half done\n";

//    children[0].xStar[mutated_gene_a] = mutate(mutated_gene_a, children[0].xStar[mutated_gene_a]);
//    children[1].xStar[mutated_gene_b] = mutate(mutated_gene_b, children[1].xStar[mutated_gene_b]);


	//Apply mutation
	for(int j = 0; j < inst.n_facilities ; j++){

		if(rand() < RAND_MAX/inst.n_facilities)
			children[0].xStar[ypos(j, &inst)] = mutate(j, b.xStar[ypos(j, &inst)], gene_pool);//, inst, children[0]);
//        }else
//            children[0].xStar[ypos(j, &inst)] = a.xStar[ypos(j, &inst)];

		if(rand() < RAND_MAX/inst.n_facilities)
			children[1].xStar[ypos(j, &inst)] = mutate(j, a.xStar[ypos(j, &inst)], gene_pool);//, inst, children[1]);
//      else
//          children[1].xStar[ypos(j, &inst)] = b.xStar[ypos(j, &inst)];
	}

	//Fitness computing
	for(int j = 0 ; j < inst.n_facilities ; j++){
		if(children[0].xStar[ypos(j, &inst)] == 1.0){
			children[0].z_opt += inst.fixed_costs[j];
			children[0].non_zero_count++;
		}
		if(children[1].xStar[ypos(j, &inst)] == 1.0){
			children[1].z_opt += inst.fixed_costs[j];
			children[1].non_zero_count++;
		}
	}

	for(int i = 0 ; i < inst.n_clients ; i++){

		double Wi_0 = intFunctionI(inst, children[0].xStar, i);
		double Wi_1 = intFunctionI(inst, children[1].xStar, i);
		//cout << "W(" << i << "):" <<  Wi_0;
		//cout << "W(" << i << "):" <<  Wi_1;
		children[0].xStar[wpos(i, &inst)] = Wi_0;
		children[0].z_opt += Wi_0;
		children[1].xStar[wpos(i, &inst)] = Wi_1;
		children[1].z_opt += Wi_1;
	}

	return children;
}


/*
 * mutation: probability of change one gene
 */
Solution* slowCrossover(Instance& inst, Solution a, Solution b, double mutation){       // ,int nzcnt_best          //mutation= quantity  of mutations

    Solution* children = new Solution[2];
    children[0].xStar = new double[inst.n_facilities + inst.n_clients];
    children[1].xStar = new double[inst.n_facilities + inst.n_clients];
    children[0].z_opt = 0.0;
    children[1].z_opt = 0.0;


    int cut_point = rand()%inst.n_facilities;
    int mutated_gene = rand()%inst.n_facilities;
    int gene_idx;

    for(int j = 0 ; j < cut_point; j++){

        //if(j == mutated_gene){                                            //
        if(rand() < RAND_MAX/inst.n_facilities){//mutation){             //
            children[0].xStar[ypos(j, &inst)] = mutate(j, b.xStar[ypos(j, &inst)]);
            children[1].xStar[ypos(j, &inst)] = mutate(j, a.xStar[ypos(j, &inst)]);
        } else {
            children[0].xStar[ypos(j, &inst)] = b.xStar[ypos(j, &inst)];
            children[1].xStar[ypos(j, &inst)] = a.xStar[ypos(j, &inst)];
        }
    }

    for(int j = cut_point; j < inst.n_facilities ; j++){

        //if(j == mutated_gene){
        if(rand() < RAND_MAX/inst.n_facilities){
            children[0].xStar[ypos(j, &inst)] = mutate(j, b.xStar[ypos(j, &inst)]);
            children[1].xStar[ypos(j, &inst)] = mutate(j, a.xStar[ypos(j, &inst)]);
        } else {
            children[0].xStar[ypos(j, &inst)] = a.xStar[ypos(j, &inst)];
            children[1].xStar[ypos(j, &inst)] = b.xStar[ypos(j, &inst)];
        }
    }

    for(int j = 0 ; j < inst.n_facilities ; j++){
        if(children[0].xStar[ypos(j, &inst)] == 1.0){
            children[0].z_opt += inst.fixed_costs[j];
            children[0].non_zero_count++;
        }
        if(children[1].xStar[ypos(j, &inst)] == 1.0){
            children[1].z_opt += inst.fixed_costs[j];
            children[1].non_zero_count++;
        }
    }


    for(int i = 0 ; i < inst.n_clients ; i++){

        double Wi_0 = intFunctionI(inst, children[0].xStar, i);
        double Wi_1 = intFunctionI(inst, children[1].xStar, i);
        //cout << "W(" << i << "):" <<  Wi << endl;
        children[0].xStar[wpos(i, &inst)] = Wi_0;
        children[0].z_opt += Wi_0;
        children[1].xStar[wpos(i, &inst)] = Wi_1;
        children[1].z_opt += Wi_1;
    }

    return children;
}



