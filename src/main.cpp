//#include <ilcplex/cplex.h>
//#include <ilcplex/ilocplex.h>
//#include <ilconcert/ilomodel.h>
#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <istream>

#include "ufl.h"
#include "genetics.h"
#include "cuda_functions.h"
#include <exception>
#include <stdexcept>


#define VERBOSE 					0
#define VERBOSE_OUTPUT				1


const char* FILENAME = "B1.1";
char* PATH;

using namespace std;

/*
void CppApiTest (IloModel model, IloNumVarArray x, IloRangeArray c, int** rows){

     IloEnv env = model.getEnv();

     x.add(IloNumVar(env, 0.0, 40.0));
     x.add(IloNumVar(env));
     x.add(IloNumVar(env));
     model.add(IloMaximize(env, x[0] + 2 * x[1] + 3 * x[2]));

     c.add( - x[0] +     x[1] + x[2] <= 20);
     c.add(   x[0] - 3 * x[1] + x[2] <= 30);
     model.add(c);

}
*/

void parseSimple(const char* PATH, Instance& instance, int** rows){

	printf("Start parsing simple format\n");

	cout << "Parsing: " << PATH << endl;
	ifstream file (PATH);

	string line;
	string first_line;
	string second_line, third_line;

	int n;	//	# of Plants
	int m;	//	# of Clients

	double** rows_t;

	if ( file.is_open() ){

		//Get first line with header and infos
		getline ( file, first_line);

		//Get M, N of the problem from first line
		getline ( file, second_line);
		istringstream is( second_line );
		is >> n >> m;

		printf("#Facilities(N): %d, #Clients(M): %d\n", n, m );

		rows_t = new double*[n];

		//Read it line by line
		int j = 0;
		while ( getline (file, line) ){

			//cout << line << '\n';
			  double* facility = new double[m+2];
			  istringstream fac_is( line );

			  //cout << "cacca" << endl;

			  //Parse line argument, first two are facility number and fixed cost for particular facility
			  for(int i = 0; i < m + 2; i++){

				  fac_is >> facility[i];

				  //cout << "facility[" << i << "]:" << facility[i] <<  endl;
			  }

			  rows_t[j] = facility;
			  j++;
			}

		file.close();
		}

	printf("Done Parsing, creating instance for simple problem\n");

	//Instance inst;
	instance.n_facilities = n;
	instance.n_clients = m;
	instance.fixed_costs = new double[instance.n_facilities];
	instance.costs = new double[instance.n_facilities*instance.n_clients];

	instance.t_costs = new double[instance.n_clients * instance.n_facilities];

	rows = new int*[instance.n_clients];

	for(int j = 0; j < instance.n_clients; j++){
		rows[j] = new int[instance.n_facilities];
	}

	for(int i = 0; i < n; i++){

		//cout << rows_t[i][1] << endl;
		instance.fixed_costs[i] = rows_t[i][1];

		for(int j = 0; j < m; j++){

			if(VERBOSE){
				cout << "(" << i << "," << j << ") v:" << rows_t[i][j+2] << " ";
			}
			rows[j][i] = rows_t[i][j+2];
		}
		if(VERBOSE ) cout << endl;
	}

	printf("Checking\n");

	for(int j = 0; j < instance.n_clients; j++){

		//cout << "riga:" << j << endl;

		for(int i = 0; i < instance.n_facilities; i++){

			instance.costs[j*instance.n_facilities + i] = double(rows[j][i]);
			instance.t_costs[i*instance.n_clients + j] =  double(rows[j][i]);

			if (VERBOSE){
				cout << instance.costs[j*instance.n_facilities + i] << " ";
				//cout << rows[i][j] << " ";
				//cout << rows_t[i][j] << " ";
			}


		}
		if(VERBOSE) cout << endl;
	}

}

void parseOrlLib(const char* PATH, int** rows, Instance& inst){

	printf("Start parsing OrlLib format\n");

	cout << "PArsing: " << PATH << endl;
	ifstream file (PATH);

	string line;
	string first_line;
	//string second_line, third_line;

	int n;	//	# of Plants
	int m;	//	# of Clients

	int** rows_t;

	if ( file.is_open() ){
		//Get first line with header and infos
		//getline ( file, first_line);

		//Get N, M of the problem from first line
		getline ( file, first_line);
		istringstream is( first_line );
		is >> n >> m;

		printf("#Facilities(N): %d, #Clients(M): %d\n", n, m );

		//Instance inst;
		inst.n_facilities = n;
		inst.n_clients = m;
		inst.fixed_costs = new double[inst.n_facilities];
		inst.costs = new double[inst.n_facilities*inst.n_clients];

		inst.t_costs = new double[inst.n_clients * inst.n_facilities];

		printf("Reading lines\n");

		int i = 0;
		while ( getline(file, line) && i < n){

			istringstream iss(line);
			double fc, foo;

		    while(iss >> foo >> fc){

		    	//cout << fc << endl; // typing all the terms
		    	inst.fixed_costs[i] = fc;
		    }
		    i++;
		}

		for(int i = 0; i < n ; i++){

			string line, fooline;

			//Discard odd lines
			getline(file, line);
			getline(file, fooline);

			istringstream fac_line(line);
			for(int j = 0; j < m; j++){
				double Cij;
				fac_line >> Cij;
				inst.costs[i*inst.n_clients + j] = Cij;
				inst.t_costs[j*inst.n_facilities + i] = Cij;
			}
		}

		if( VERBOSE ){
			//Checking
			for(int i = 0; i < n; i++){
				cout << inst.fixed_costs[i] << endl;
			}

			for(int i = 0; i < n; i++){
				for(int j = 0; j < m; j++){
					cout <<  inst.costs[i*m + j] << " ";
				}
				cout << endl;
			}
		}
	} else {
		printf("Couldn't find input file!\n");
		exit(0);
	}
}

int main(int argc, char* argv[]) {

	printf("Hello, UFL\n");

	struct Instance inst;
	struct Solution sol;

	//Use an arbitrary file from the dataset if none is given
	const char* PATH = "data/M/P/MP1";

	if( argc >= 2){
		PATH = argv[1];
	}

	int** rows;
	ifstream file (PATH);
	string first_line;
	if ( file.is_open() ){
		//Get first line with header and infos
		//getline ( file, first_line);

		//Get N, M of the problem from first line
		getline ( file, first_line);
		if(first_line[0] == 'F'){
			parseSimple(PATH, inst, rows);
		} else
			parseOrlLib(PATH, rows, inst);
	}

	if(argc>2){

		printf("Detected multiple cmd-line arguments\n");
		solve_problem(inst, sol);

	} else {
		solve_problem(inst, sol);
	}
	 // Print the solution
	clog << "Optimal solution found: [ " << sol.z_opt << " ]" << endl;
	clog << "Variables:" << endl;

	if(VERBOSE_OUTPUT){

		for (long j = 0; j<inst.n_facilities; j++) {
			if ( sol.xStar[ypos(j, &inst)] == 0 )
				continue; // skip if zero

			cout << "y(" << j+1 << ") = " << sol.xStar[ypos(j, &inst)] << endl;
		}
		if(true){	//BENDERS_ON

			for (long i = 0; i<inst.n_clients; i++) {
				//if ( sol.xStar[inst.n_facilities + i] == 0.0 )
				//	continue; // skip if zero

				//clog << "w(" << i+1 << ") = " << sol.xStar[wpos(i, &inst)] << endl;
			}

		} else {
			for (long i = 0; i<inst.n_clients; i++) {
				for (long j = 0; j<inst.n_facilities; j++) {
					if ( sol.xStar[xpos(i,j,&inst)] == 0 )
						continue; // skip if zero

					cout << "x(" << i << ", " << j+1 << ") = " << sol.xStar[xpos(i,j,&inst)] << endl;
				}
			}
		}
	}

	cout << "The remaining variables are zero." << endl;

	// Print statistics to standard output
	cout << "STAT;" << argv[1] << " ; " << sol.z_opt << endl;

	return 0;
}
