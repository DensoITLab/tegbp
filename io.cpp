#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <random>
#include "tegbp.hpp"

void load_data(mem_pool pool){
    printf("Generating dummy data");
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<uint16> rand_idx(1, pool.W-2);
	std::uniform_real_distribution<> rand_val(0.0, 1.0);
	for(int i=0;i<pool.B;i++){
		pool.indices[2*i+0] = rand_idx(gen);
		pool.indices[2*i+1] = rand_idx(gen);
		pool.v_norms[2*i+0] = rand_val(gen);
		pool.v_norms[2*i+1] = rand_val(gen);
		// printf("([%3d] %2d, %2d,  %0.1f, %0.1f)\n",i, indices[2*i+0], indices[2*i+0], v_norms[2*i+0], v_norms[2*i+1]);
	}
	printf("done..\n");
    return;
}

void save_data(mem_pool pool){
	printf("Saving data");
	for(uint16 y=0;y<pool.H/6;y++){
		printf("\n");
		for(uint16 x=0;x<pool.W/6;x++){
			printf("%04.1f|", pool.node[sub2ind(x, y, 0, 0, pool.W, pool.H)]);
		}
	}
	printf("done..\n");
}