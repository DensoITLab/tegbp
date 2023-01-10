#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <omp.h>
#include <random>
#include <unistd.h>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include "plane_fitting.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

int32 dircPF[2][N_PATCH] = {-1, 0, +1, -1, 0, +1, -1, 0, +1, -1, -1, -1,  0, 0, 0, +1, +1, +1};
V9D Bx = (V9D() << -1.0, 0.0, 1.0, -1.0, 0.0, 1.0, -1.0, 0.0, 1.0).finished();
V9D By = (V9D() << -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0).finished();
V9D B1 = (V9D() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0).finished();
Matrix<double, 9, 3> base(N_PATCH, BASE);

// #pragma omp declare simd
inline int32 sub2ind_(int32 x, int32 y, int32 H, int32 W){
	return (x + W*y);
}

M3D* get_baseinv(double *BtBinv_lut, int32 ind)
{
	return (M3D*)&BtBinv_lut[ind];
}

void set_vperp(double *flow_norm, int32 ind, V2D* data)
{
	V2D* ptr = (V2D*)&flow_norm[ind];
	*ptr = *data;
	// memcpy(&node[ind], data, sizeof(V6D));
	return;
}

/* binary to index */
void bit2ind(bool valid[N_PATCH], int32 *ind)
{
	*ind = 0;
	for(int32 dir=0; dir<N_PATCH; dir++){
		if (valid[dir]){
			*ind = *ind + (int32)pow(2, dir);
		}
	}
}

/* temporal window for plane fitting */
// #pragma omp declare simd
bool isValidPF(int32 ind, int32 t, int32* sae)
{
	return ((t - sae[ind]) < WIN_T);
}

/* compute valid event in a patch */
void precompute_idx_pf(mem_pool pool, XYTPI* xytpi, int32 xyis[][N_PATCH], bool valid[N_PATCH], int32 *n_valid){
	// https://qiita.com/hareku/items/3d08511eab56a481c7db
	*n_valid = 0;
	int32 ind_pol;
	// #pragma simd
	for(int32 dir=0; dir<N_PATCH; dir++){
		xyis[0][dir] 	= (*xytpi).at(0) + dircPF[0][dir];
		xyis[1][dir] 	= (*xytpi).at(1) + dircPF[1][dir];
		xyis[2][dir] 	= sub2ind_(xyis[0][dir], xyis[1][dir], pool.H, pool.W);
		// polarity
		valid[dir] 		= isValidPF(2*xyis[2][dir]+(*xytpi).at(3), (*xytpi).at(2), pool.sae_pol);
		if (valid[dir]){
			*n_valid 	= *n_valid + 1;;
		}
	}
}

void plane_fitting(mem_pool pool,  XYTPI* xytpi) // plane fitting
{
	int32 xyi_pf[4][N_PATCH];
	bool valid[N_PATCH];
	int32 n_valid;
	int32 ind_BtBinv;

	/* iterativeにやるならvalidを変えながら繰り返し */
	/* calculate valid event in a patch */
	precompute_idx_pf(pool, xytpi, xyi_pf, valid, &n_valid);
	for(int32 dir=0; dir<N_PATCH; dir++){
		printf("%d", valid[dir]);
	}
	printf("\n"); 						// 000000000
	// printf("%d\n", n_valid);

	bit2ind(valid, &ind_BtBinv);		// binary pattern to index
	printf("%d\n", ind_BtBinv); 		// 0~2^9

	if (n_valid < N_TH_PF){				// require a sufficient number of valid events
		return;
	}
	// look up the inverse of base for valid pixels.
	M3D BtBinv = *get_baseinv(pool.BtBinv_lut, ind_BtBinv*N_PATCH);
	// printf("%f,%f,%f\n", BtBinv(0, 0), BtBinv(0, 1), BtBinv(0, 2));
	// printf("%f,%f,%f\n", BtBinv(1, 0), BtBinv(1, 1), BtBinv(1, 2));
	// printf("%f,%f,%f\n", BtBinv(2, 0), BtBinv(2, 1), BtBinv(2, 2));
	// printf("============\n");

	// calculate Btf for valid pixels in a patch
	int f0 = 0;
	int f1 = 0;
	int f2 = 0;
	for(int32 dir=0; dir<N_PATCH; dir++){
		if (valid[dir]){
			// 0, 1, -1だから掛け算要らないので書き直す
			f0 += Bx(dir) * pool.sae[xyi_pf[2][dir]];
			f1 += By(dir) * pool.sae[xyi_pf[2][dir]];
			f2 += pool.sae_pol[2*xyi_pf[2][dir]+(*xytpi).at(3)];
		}
	}

	V3D Btf 	= (V3D() << (double)f0, (double)f1, (double)f2).finished();
	V3D abc 	= BtBinv * Btf; 				// r = (BtB)^-1 Btf
	// printf("%d, %d, %d\n", f0, f1, f2);
	// printf("%f, %f, %f\n", abc(0), abc(1), abc(2));
	// printf("============\n");

	double a2b2 	= pow(abc(0), 2) + pow(abc(1), 2);
	if (a2b2 == 0){ // nan flow
		return;
	}
	V2D v_perp 		= (V2D() << abc(0) / a2b2 * DT, abc(1) / a2b2 * DT).finished(); // pix/DT
	printf("%d, %d, %d, %d, %d\n", (*xytpi).at(0), (*xytpi).at(1), (*xytpi).at(2), (*xytpi).at(3), (*xytpi).at(4));
	printf("%f, %f\n", v_perp(0), v_perp(1));
	set_vperp(pool.flow_norm, (*xytpi).at(4)*2, &v_perp);

	return;
}

mem_pool initialize(mem_pool pool)
{
	pool.BtBinv_lut = (double *) malloc((1<<N_PATCH)*sizeof(M3D)); // 2^9 x (3 x 3)
	double cond = 0.00001;
	M3D C = (M3D() << cond, 0.0, 0.0, 0.0, cond, 0.0, 0.0, 0.0, cond).finished();
	base << Bx, By, B1;
	for (int bit=0; bit<(1<<N_PATCH); bit++){
		V9D Wv;
		for (int i=0; i<N_PATCH; i++){ // binary for each pix
			if (bit & (1<<i)){ // i-th element
				Wv(i, 0) = 1.0;
			}else{
				Wv(i, 0) = 0.0;
			}
		}
		M9D Wm 			= Wv.array().matrix().asDiagonal();
		M3D BtB 		= base.transpose() * Wm * base;
		M3D BtBinv 	= (BtB+C).inverse();
		int32 ind 	= bit*N_PATCH; // index of each pattern
		M3D *ptr 		= (M3D*)&(pool.BtBinv_lut[ind]);
		*ptr 				= BtBinv;

		// M3D BtBinv_ = *get_baseinv(pool.BtBinv_lut, ind);
		// cout << BtBinv_ << "\n";
		// cout << "---------" << "\n";
		// nanやinfを0にしておけば、勝手にv_perpが0になってくれるかも
	}

	pool.flow_norm  = (double *) malloc(2*pool.W*pool.H*sizeof(double));
	pool.sae  	    = (int32 *) malloc(1*pool.W*pool.H*sizeof(int32));
	memset(pool.sae, -10*DT_ACT, 1*pool.W*pool.H*sizeof(int32));
	pool.sae_pol    = (int32 *) malloc(2*pool.W*pool.H*sizeof(int32));
	memset(pool.sae_pol, -10*DT_ACT, 2*pool.W*pool.H*sizeof(int32));
	pool.v_norms 	= (double *) malloc(2*pool.B*sizeof(double));
	pool.indices    = (int32 *) malloc(2*pool.B*sizeof(int32));
	pool.timestamps = (int32 *) malloc(1*pool.B*sizeof(int32));
	pool.polarities = (int32 *) malloc(1*pool.B*sizeof(int32));
	printf("initialize memory pool %d  %d  %d\n", pool.B, pool.H, pool.W);
  return pool;
}

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

void load_data_event(mem_pool pool, std::string data_name){
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
		pool.v_norms[2*n+0] = 0.0; // vx_perp
		pool.v_norms[2*n+1] = 0.0; // vy_perp
		pool.timestamps[n] 	= (int32)stoi(strvec.at(2)); // t
		pool.polarities[n] 	= (int32)stoi(strvec.at(3)); // p
		// printf("timestamps %d %d\n", n, pool.timestamps[n] );
		n++;
    }
	// printf("done..\n");
    return;
}

void save_data_pf(mem_pool pool, int32 seq_id, int32 index, int32 c_time){
	double *fimg	= (double *) malloc(2*pool.W*pool.H*sizeof(double));
	memset(fimg, 0.0, 2*pool.W*pool.H*sizeof(double));

	int32 time;
	#pragma omp parallel for
	for(int32 y=0; y<pool.H; y++){
		for(int32 x=0;x<pool.W;x++){
			time = pool.sae[(pool.W*y + x)];
			if ((c_time-time)<DT_ACT){
				int32 ind = sub2ind_(x, y, pool.H, pool.W);
				fimg[2*ind] 	= pool.flow_norm[2*ind];
				fimg[2*ind+1] 	= pool.flow_norm[2*ind+1];
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

mem_pool load_data_pf(std::string data_name){
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

	load_data_event(pool, data_name);

	printf("done..\n");
	return pool;
}

void process_batch_pf(mem_pool pool, int32 b_ptr)
{
  for(int32 i=b_ptr; i<(b_ptr+pool.WINSIZE); i++){
		int32 x = pool.indices[2*i+0];
		int32 y = pool.indices[2*i+1];
        int32 t = pool.timestamps[i];
        int32 p = pool.polarities[i];
		int32 ind = sub2ind_(x, y, pool.H, pool.W);

		XYTPI xytpi = {x, y, t, p, ind};
		// printf("([%3d] %3.2f, %3.2f, %2d, %2d, %3d)  thread: %d\n",i, v_perp(0), v_perp(1), x, y, t, omp_get_thread_num());

		pool.sae[ind] 			= t;
		pool.sae_pol[2*ind+p] 	= t; // 極性ごとのSAE

		// Core of message passing
		plane_fitting(pool, &xytpi);
	}
	return;
}