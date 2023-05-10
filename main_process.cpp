#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <omp.h>
#include <random>
#include <unistd.h>
#include "tegbp.hpp"
#include "progressbar.hpp"
#include <sstream>

using namespace std;

int32 main(int32 argc, char **argv)
{
	// Setup OMP
	int32 num_thread 		= omp_get_max_threads();
	int32 WINSIZE 		    = 1000;
	int32 n_itr_save 		= 1;
	std::string data_name 	= "bricks";
	bool is_full 			= true; // true, false

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

	if (argc>5){
		is_full = (bool)atoi(argv[5]);
	}
	printf("Start main dataset: %s WINSIZE %d num_thread %d  is_full %d \n", data_name.c_str(), WINSIZE, num_thread, is_full);
	// Load data and Initialize Global varialble 	by Allocating memory
	data_cfg cfg = load_cfg(data_name);
	mem_pool pool;
	if (is_full){
		pool = initialize_full(cfg);
	}else{
		pool = initialize_nrml(cfg);
	}

	pool = load_data(pool);

	pool.WINSIZE = WINSIZE;

	// Run the main image processing function
	double ellapse 		= 0;
	int32 b_ptr 		= 0;
	int32 n_itr 		= pool.B/WINSIZE-1; // shoud be B/WINSIZE
	int32 increment 	= 1;
	if (pool.WINSIZE>=pool.B){
		// For benchmarking
		n_itr = 100;
		increment=0;
	}

	int32 s_time 		= 0;
	int32 e_time 		= 0;
	double start = omp_get_wtime();
	progressbar bar(n_itr);

	std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

	if (is_full){
		for (int32 itr=0; itr<n_itr; itr++){
			s_time = pool.timestamps[b_ptr];
			e_time = pool.timestamps[b_ptr+WINSIZE];
			// bar.update(format(" %dK - %dK  %5.2f - %5.2f",b_ptr/1000, (b_ptr+WINSIZE)/1000, (double)s_time/1000000.0, (double)e_time/1000000.0));
			// printf("%dK - %dK  %5.2f - %5.2f\n",b_ptr/1000, (b_ptr+WINSIZE)/1000, (double)s_time/1000000.0, (double)e_time/1000000.0);

			#pragma omp parallel
			{
				process_batch_full(pool, b_ptr);
			}

			b_ptr 	= b_ptr + (WINSIZE*increment); 

			if (n_itr_save>0 & itr%n_itr_save==0){
				save_flow(pool, itr, 0, e_time);			// estimated flow
				save_img(pool, itr, 0, b_ptr+WINSIZE); 		// image
				save_img(pool, itr, 1, b_ptr+WINSIZE); 		// image
			}
		}
	}else{
		for (int32 itr=0; itr<n_itr; itr++){
			s_time = pool.timestamps[b_ptr];
			e_time = pool.timestamps[b_ptr+WINSIZE];
			// bar.update(format(" %dK - %dK  %5.2f - %5.2f",b_ptr/1000, (b_ptr+WINSIZE)/1000, (double)s_time/1000000.0, (double)e_time/1000000.0));
			printf("%dK - %dK  %5.2f - %5.2f\n",b_ptr/1000, (b_ptr+WINSIZE)/1000, (double)s_time/1000000.0, (double)e_time/1000000.0);

			#pragma omp parallel
			{
				process_batch_nrml(pool, b_ptr);
			}

			b_ptr 	= b_ptr + (WINSIZE*increment); 

			if (n_itr_save>0 & itr%n_itr_save==0){
				save_flow(pool, itr, -1, e_time); // estimated flow
			}
		}
	}
	double end = omp_get_wtime();
	ellapse = ellapse + (end-start);
	printf("\n%5.3f [sec] for %d [K events] %5.3f [sec] (num_thread: %d) %5.3f KEPS\n",ellapse, pool.B/1000, (double)s_time/1000000.0,  num_thread, (pool.B/1000)/ellapse);

	// Visualize results
	// debug_output(pool);

    return 0;
}