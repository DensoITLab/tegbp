#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <omp.h>
#include <random>
#include <unistd.h>
#include "tegbp.hpp"

int32 main(int32 argc, char **argv)
{
	// Setup OMP
	int32 num_thread 		= omp_get_max_threads();
	int32 WINSIZE 		    = 1000;
	int32 n_itr_save 		= 1;
	std::string data_name 	= "bricks";

	if (argc>1){
		if (atoi(argv[1])>0){
			num_thread = atoi(argv[1]);
		}
	}
	omp_set_num_threads(num_thread);

	if (argc>2){
		data_name = argv[2];
		printf("Dataset: %s\n", data_name.c_str());
	}

	if (argc>3){
		if (atoi(argv[3])>0){
			WINSIZE = atoi(argv[3]);
		}
	}

	if (argc>4){
		n_itr_save = atoi(argv[4]);
	}

	printf("Start main dataset: %s WINSIZE %d\n", data_name.c_str(), WINSIZE);

	// Initialize Global varialble 	by Allocating memory
	// Load data
	mem_pool pool = load_data(data_name);
	pool.WINSIZE = WINSIZE;

	// init_sae(pool); // should update sae in process batch

	// Run the main image processing function
	double ellapse 	= 0;
	int32 b_ptr 		= 0;
	int32 n_itr 		= pool.B/WINSIZE; // shoud be B/WINSIZE
	int32 increment 	= 1;
	if (pool.WINSIZE>=pool.B){
		// For debug
		n_itr = 100;
		increment=0;
	}

	int32 c_time 		= 0;
	for (int32 itr=0; itr<n_itr; itr++){
		double start = omp_get_wtime();
		printf("process_batch %d - %d\n",b_ptr, (b_ptr+WINSIZE));
		#pragma omp parallel
		{
		process_batch(pool, b_ptr);
		}
		c_time = pool.timestamps[b_ptr+WINSIZE];
 		double end = omp_get_wtime();

		b_ptr 	= b_ptr + (WINSIZE*increment); 
		ellapse = ellapse + (end-start);

		if (n_itr_save>0 & itr%n_itr_save==0){
			save_data(pool, itr, 0, c_time); // estimated flow
			// save_data(pool, itr, 1, c_time); // normal flow
		}
	}
	printf("Work took %f seconds for %d K events (num_thread: %d) %5.3f KEPS\n",ellapse, pool.B/1000, num_thread, (pool.B/1000)/ellapse);

	// Visualize results
	// debug_output(pool);

    return 0;
}