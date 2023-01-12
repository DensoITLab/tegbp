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

// void load_data_dummy(mem_pool pool){
//     // printf("Generating dummy data");
// 	std::random_device rd;
// 	std::mt19937 gen(rd());
// 	std::uniform_int_distribution<int32> rand_idx(32, pool.W-32);
// 	std::uniform_int_distribution<int32> rand_time(1000, 1000000);
// 	std::uniform_real_distribution<> rand_val(0.0, 1.0);
// 	for(int32 i=0;i<pool.B;i++){
// 		pool.indices[2*i+0] = rand_idx(gen);
// 		pool.indices[2*i+1] = rand_idx(gen);
// 		pool.v_norms[2*i+0] = rand_val(gen);
// 		pool.v_norms[2*i+1] = rand_val(gen);
// 		pool.timestamps[i] 	= rand_time(gen);
// 		// printf("([%3d] %2d, %2d,  %0.1f, %0.1f)\n",i, indices[2*i+0], indices[2*i+0], v_norms[2*i+0], v_norms[2*i+1]);
// 	}
//     return;
// }


void load_txt(mem_pool pool){
	std::string data_path = "data/" + pool.data_name + ".txt";
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
		pool.timestamps[n] 	= (int32)stoi(strvec.at(2)); // t
		pool.polarities[n] 	= (int32)stoi(strvec.at(3)); // p
		pool.norms[2*n+0] 	= stod(strvec.at(4)); // vx_perp
		pool.norms[2*n+1] 	= stod(strvec.at(5)); // vy_perp
		pool.timestamps[n] 	= (int32)stoi(strvec.at(2)); // t
		// printf("timestamps %d %d\n", n, pool.timestamps[n] );
		n++;
    }
	// printf("done..\n");
    return;
}


data_cfg load_cfg(std::string data_name){
	data_cfg cfg;
	printf("load config %s\n", data_name.c_str());

	if (data_name == "bricks"){
		cfg.W = 272;
    	cfg.H = 208;
    	cfg.B = 15000;
	}
	if (data_name == "bricks_1slide"){
		cfg.W = 272;
    	cfg.H = 208;
    	cfg.B = 15000;
	}
	if (data_name == "indoor_flying2"){
		cfg.W = 368;
    	cfg.H = 288;
    	cfg.B = 2921002;
	}
	if (data_name == "dummy"){
		cfg.W = 272;
    	cfg.H = 208;
    	cfg.B = 15000;
	}
	if (data_name == "qr000"){ // substr(0, 2) == "qr"
		cfg.W = 1312;
    	cfg.H = 752;
    	cfg.B = 1418870;
	}
	if (data_name == "qr001"){
		cfg.W = 1312;
    	cfg.H = 752;
    	cfg.B = 684966;
	}
	cfg.data_name = data_name;

	return cfg;
}

mem_pool load_data(mem_pool pool){

	// printf("initialize memory pool %d  %d  %d %d\n", pool.B, pool.H, pool.W, pool.timestamps);

	load_txt(pool);

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

void save_flow(mem_pool pool, int32 seq_id, int32 index, int32 c_time){
	double *fimg	= (double *) malloc(2*pool.W*pool.H*sizeof(double));
	memset(fimg, 0.0, 2*pool.W*pool.H*sizeof(double));

	// for(int32 y=0; y<pool.H; y++){
	// 	for(int32 x=0;x<pool.W;x++){
	// 		fimg[2*(pool.W*y + x)+0] 	= 0;
	// 		fimg[2*(pool.W*y + x)+1] 	= 0;
	// 	}
	// }
	int32 time;
	#pragma omp parallel for
	for(int32 y=0; y<pool.H; y++){
		for(int32 x=0;x<pool.W;x++){
			time = pool.sae[(pool.W*y + x)];
			if ((c_time-time)<DT_ACT){
				switch (index){
				case 0:
					fimg[2*(pool.W*y + x)] 		= pool.node[sub2ind_nod(x, y, pool.H, pool.W)+0];
					fimg[2*(pool.W*y + x)+1] 	= pool.node[sub2ind_nod(x, y, pool.H, pool.W)+1];
					break;
				case 1:
					fimg[2*(pool.W*y + x)] 		= pool.node[sub2ind_nod(x, y, pool.H, pool.W)+0+STS_DIM];
					fimg[2*(pool.W*y + x)+1] 	= pool.node[sub2ind_nod(x, y, pool.H, pool.W)+1+STS_DIM];
					break;
				default:
					int32 ind = sub2ind_(x, y, pool.H, pool.W);
					fimg[2*ind] 	= pool.flow_norm[2*ind];
					fimg[2*ind+1] 	= pool.flow_norm[2*ind+1];
					break;
				}
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


void save_img(mem_pool pool, int32 seq_id, int32 index, int32 c_idx){
	int32 *fimg	= (int32 *) malloc(pool.W*pool.H*sizeof(int32));
	for(int32 y=0; y<pool.H; y++){
		for(int32 x=0;x<pool.W;x++){
			fimg[pool.W*y + x] 		= 1;
		}
	}
	// memset(fimg, 1.0, pool.W*pool.H*sizeof(double));

	int32 c_time 	= pool.timestamps[c_idx];
	int32 idx 		= 0;
	for(int32 i=c_idx; i>0; i--){
		if (c_time - pool.timestamps[i] > DT_ACT){
			idx = i+1;
			break;
		}
	}
	// printf("%d, %d\n", idx, c_idx);
	double scale;
	double x_to;
	double y_to;
	for(int32 i=idx; i<=c_idx; i++){
		int32 x 	= pool.indices[2*i+0];
		int32 y 	= pool.indices[2*i+1]; 
        int32 t 	= pool.timestamps[i];	
        int32 p 	= pool.polarities[i];				
		double vx 	= pool.fulls[2*i+0];
		double vy 	= pool.fulls[2*i+1];
		switch (index){
			case 0:
				fimg[pool.W*y + x] 		= p*2;
				break;
			case 1:
				if ((vx==0.0) && (vy==0.0)) break;
				scale 	= (double)(c_time - t) / (double)DT;
				x_to  	= (double)x + scale * vx;
				y_to  	= (double)y + scale * vy;
				if (x_to >= pool.W) x_to = pool.W-1;
				if (y_to >= pool.H) y_to = pool.H-1;
				if (x_to < 0) x_to = 0;
				if (y_to < 0) y_to = 0;
				// printf("%.1f, %.1f, %d, %d, %f\n", vx, vy, t, c_time, scale);
				// printf("%d, %d\n", (int32)x, (int32)y);
				// printf("%d, %d\n", (int32)x_to, (int32)y_to);
				fimg[pool.W*(int32)y_to + (int32)x_to] = p*2;
				break;
			default:
				break;
		}
	}

	stringstream filename;
	filename << "result/" << pool.data_name << "/bin/img_"  << index <<  "_" << std::setw(5) << std::setfill('0') << seq_id << ".bin";
	std::ofstream myFile (filename.str(), std::ios::out | std::ios::binary);
	myFile.write ((char *)fimg, pool.W*pool.H*sizeof(int32));
	printf("Saving data %s, %d_%04d\n",filename.str().c_str(), index, seq_id);
	// printf("comleted\n");
}