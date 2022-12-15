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
	int num_thread 		= 64;
	int W 			    = 272;
	int H 			    = 208;
	int B 		        = 10000;

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
	load_data_brick(pool);

	init_sae(pool);

	// Run the main image processing function
	double ellapse = 0;
	int b_ptr = 0;
	for (int itr=0;itr<10;itr++){
		double start = omp_get_wtime();
		process_batch(pool, b_ptr);
 		double end = omp_get_wtime();

		b_ptr = b_ptr + (WINSIZE*0); 
		ellapse = ellapse + (end-start);
		save_data(pool, itr, 0);
		save_data(pool, itr, 1);
	}
	printf("Work took %f seconds for %d K events (num_thread: %d)\n",ellapse, B/1000, num_thread);

	// Visualize results
	// debug_output(pool);

    return 0;
}