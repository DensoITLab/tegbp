#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <fstream>
#include "tegbp.hpp"
#include <iomanip>
using namespace std;

vector<string> split(string& input, char delimiter)
{
    istringstream stream(input);
    string field;
    vector<string> result;
    while (getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}

void load_data_dummy(mem_pool pool){
    // printf("Generating dummy data");
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int32> rand_idx(32, pool.W-32);
	std::uniform_int_distribution<int32> rand_time(1000, 1000000);
	std::uniform_real_distribution<> rand_val(0.0, 1.0);
	for(int32 i=0;i<pool.B;i++){
		pool.indices[2*i+0] = rand_idx(gen);
		pool.indices[2*i+1] = rand_idx(gen);
		pool.v_norms[2*i+0] = rand_val(gen);
		pool.v_norms[2*i+1] = rand_val(gen);
		pool.timestamps[i] 	= rand_time(gen);
		// printf("([%3d] %2d, %2d,  %0.1f, %0.1f)\n",i, indices[2*i+0], indices[2*i+0], v_norms[2*i+0], v_norms[2*i+1]);
	}
    return;
}


void load_data_txt(mem_pool pool, std::string data_name){
	std::string data_path = "data/" + data_name + ".txt";
	// printf("load_data_txt %s %d\n", data_path.c_str(), pool.B);
	ifstream ifs(data_path);
	string line;

	getline(ifs, line);
	// cout << line << '\n';

	int32 n = 0;
    while (getline(ifs, line) && n<pool.B) {
        vector<string> strvec = split(line, ',');
		pool.indices[2*n+0] = (int32)stoi(strvec.at(0)); // x
		pool.indices[2*n+1] = (int32)stoi(strvec.at(1)); // y
		pool.v_norms[2*n+0] = stod(strvec.at(4)); // vx_perp
		pool.v_norms[2*n+1] = stod(strvec.at(5)); // vy_perp
		pool.timestamps[n] 	= (int32)stoi(strvec.at(2)); // t
		// printf("timestamps %d %d\n", n, pool.timestamps[n] );
		n++;
    }
	// printf("done..\n");
    return;
}


mem_pool load_data(std::string data_name){
	mem_pool pool;
	printf("Loading %s\n", data_name.c_str());

	if (data_name == "bricks"){
		pool.W = 272;
    	pool.H = 208;
    	pool.B = 15000;
	}
	if (data_name == "bricks_1slide"){
		pool.W = 272;
    	pool.H = 208;
    	pool.B = 15000;
	}
	if (data_name == "indoor_flying2"){
		pool.W = 368;
    	pool.H = 288;
    	pool.B = 2921002;
	}
	if (data_name == "dummy"){
		pool.W = 272;
    	pool.H = 208;
    	pool.B = 15000;
	}
	if (data_name == "qr000"){ // substr(0, 2) == "qr"
		pool.W = 1312;
    	pool.H = 752;
    	pool.B = 1418870;
	}
	if (data_name == "qr001"){
		pool.W = 1312;
    	pool.H = 752;
    	pool.B = 684966;
	}
	pool.data_name = data_name;
	pool = initialize(pool);
	// printf("initialize memory pool %d  %d  %d %d\n", pool.B, pool.H, pool.W, pool.timestamps);

	if (data_name == "dummy"){
		load_data_dummy(pool);
	}else{
		load_data_txt(pool, data_name);
	}


	printf("done..\n");
	return pool;
}


void debug_output(mem_pool pool){
	printf("Saving data");
	for(int32 y=0;y<pool.H/6;y++){
		printf("\n");
		for(int32 x=0;x<pool.W/6;x++){
			printf("%04.1f|", pool.node[sub2ind_nod(x, y,  pool.H, pool.W)]);
		}
	}
	printf("done..\n");
}

void save_data(mem_pool pool, int32 seq_id, int32 index, int32 c_time){
	double *fimg	= (double *) malloc(2*pool.W*pool.H*sizeof(double));
	memset(fimg, 0.0, 2*pool.W*pool.H*sizeof(double));

	int32 time;
	#pragma omp parallel for
	for(int32 y=0; y<pool.H; y++)
		for(int32 x=0;x<pool.W;x++){
			time = pool.sae[(pool.W*y + x)];
			if ((c_time-time)<DT_ACT){
				if (index==0){
					fimg[2*(pool.W*y + x)] 		= pool.node[sub2ind_nod(x, y, pool.H, pool.W)+0];
					fimg[2*(pool.W*y + x)+1] 	= pool.node[sub2ind_nod(x, y, pool.H, pool.W)+1];
				}else{
					fimg[2*(pool.W*y + x)] 		= pool.node[sub2ind_nod(x, y, pool.H, pool.W)+0+STS_DIM];
					fimg[2*(pool.W*y + x)+1] 	= pool.node[sub2ind_nod(x, y, pool.H, pool.W)+1+STS_DIM];
				}

			}
	}

	stringstream filename;
	filename << "result/" << pool.data_name << "/bin/flo_"  << index <<  "_" << std::setw(5) << std::setfill('0') << seq_id << ".bin";
	std::ofstream myFile (filename.str(), std::ios::out | std::ios::binary);
	myFile.write ((char *)fimg, 2*pool.W*pool.H*sizeof(double));
	printf("Saving data %s, %d_%04d\n",filename.str().c_str(), index, seq_id);
	// printf("comleted\n");
}