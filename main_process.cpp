#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <omp.h>
#include <random>
#include <unistd.h>
#include "tegbp.hpp"

int main(int argc, char **argv)
{
	// Setup OMP
	int num_thread 			= omp_get_max_threads();
	// int W 			    	= 272;
	// int H 			    	= 208;
	// int B 		        	= 30000;
	// std::string filename 	= "bricks";
	int W 			    	= 368;
	int H 			    	= 288;
	int B 		        	= 100000;
	std::string filename 	= "indoor_flying2";

	if (argc>1){
		num_thread = atoi(argv[1]);
	}
	omp_set_num_threads(num_thread);

	if (argc>2){
		B = atoi(argv[2]);
	}

	// Initialize Global varialble 	by Allocating memory
	mem_pool pool=initialize(B, H, W);

	// Load data
	// load_data_dummy(pool);
	// load_data_brick(pool);
	load_data_txt(pool, filename);

	// init_sae(pool); // should update sae in process batch

	// Run the main image processing function
	double ellapse 	= 0;
	int b_ptr 		= 0;
	int n_itr 		= 10; // shoud be B/WINSIZE
	int n_itr_show 	= 1;
	for (int itr=0; itr<n_itr; itr++){
		double start = omp_get_wtime();
		printf("process_batch %d, %d\n",b_ptr, (b_ptr+WINSIZE));
		// #pragma omp parallel
		// {
		process_batch(pool, b_ptr);
		// }
 		double end = omp_get_wtime();

		b_ptr 	= b_ptr + (WINSIZE*1); 
		ellapse = ellapse + (end-start);

		if (itr%n_itr_show==0){
			save_data(pool, itr, 0); // estimated flow
			save_data(pool, itr, 1); // normal flow
		}
	}
	printf("Work took %f seconds for %d K events (num_thread: %d)\n",ellapse, B/1000, num_thread);

	// Visualize results
	// debug_output(pool);

    return 0;
}