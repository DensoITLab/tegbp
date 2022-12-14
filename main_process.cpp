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
	int B 		        = 1000;

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
	load_data(pool);

	// Run the main image processing function
	double start = omp_get_wtime();
	process_batch(pool);
	double end = omp_get_wtime();
	printf("Work took %f seconds for %d K events (num_thread: %d)\n", end - start, B/1000, num_thread);

	// Visualize results
	save_data(pool);

    return 0;
}