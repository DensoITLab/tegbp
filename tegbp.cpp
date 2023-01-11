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
#include "tegbp.hpp"

using namespace std;

// OMP magic
// #pragma omp critical
// #pragma omp barrier
// #pragma omp barrier
// #pragma omp for nowait

int32 dirc[2][N_EDGE];
int32 dirc_idx[N_EDGE];

// self:0  obs:1  from up:2 from down:3  from left:4: from down:5
#pragma omp declare simd
int32 sub2ind_(int32 x, int32 y, int32 H, int32 W){
	return (x + W*y);
}

#pragma omp declare simd
int32 sub2ind_nod(int32 x, int32 y, int32 H, int32 W){
	return NOD_DIM*sub2ind_(x, y, H, W);
}

#pragma omp declare simd
bool isActive(int32 ind, int32 t, int32* sae)
{
	return ((t - sae[ind]) < DT_ACT);
}

// Geter and Setter
#pragma omp declare simd
V6D* get_state(double * node, int32 ind)
{
	return (V6D*)&node[ind];
}
#pragma omp declare simd
V2D* get_mu(double * node, int32 ind)
{
	return (V2D*)&node[ind];
}
#pragma omp declare simd
void set_state(double * node, int32 ind, V6D* data)
{
	// V6D* ptr = (V6D*)&node[ind];
	// *ptr = *data;
	memcpy(&node[ind], data, sizeof(V6D));
	return;
}
#pragma omp declare simd
void set_mu(double * node, int32 ind, V2D* data)
{
	// V2D* ptr = (V2D*)&node[ind];
	// *ptr = *data;
	memcpy(&node[ind], data, sizeof(V2D));
	return;
}
V2D belief_vec_to_mu(V6D belief)
{
	constexpr double cond = 0.00001;
	M2D C = (M2D() << cond, 0.0, 0.0,  cond).finished();
	V2D eta = belief.head(2);
	Map<M2D,Eigen::Aligned> Lam_0(belief.tail(4).data());
	return (Lam_0+C).inverse() * eta;
}

V6D update_state(double * node, bool* active, int32 ind_slf)
{
	V6D belief = *get_state(node, ind_slf+STS_DIM);
	for(int32 dir=0; dir<N_EDGE; dir++){
		if (active[dir]){
			belief += *get_state(node, ind_slf + (dir+2)*STS_DIM);
		}
	}

	V2D mu = belief_vec_to_mu(belief);

    // Set updated state
    set_mu(node, ind_slf, &mu);
	return belief;
}

#pragma omp declare simd
double huber_scale_(double M2){
	constexpr double huber	= 1.0;
	constexpr double huber2 = huber*huber;
	constexpr double scale	= 1.0;
	if (huber > 0){
		if (M2 < huber2){
			double scale   = 1.0;
		}else{
			double scale   = 2.0*huber/sqrt(M2) - huber2/M2;
		}
	}
	return scale;
}

double huber_scale_v2d(V2D *v_perp, M2D* Lam){
	double M2      = ((*v_perp).transpose() * (*Lam) * (*v_perp))(0);
	return huber_scale_(M2);
}

double huber_scale_v4d(V4D *state, M4D* Lam){
	double M2      = ((*state).transpose() * (*Lam) * (*state))(0);
	return huber_scale_(M2);
}

void set_observation(mem_pool* pool, V2D* v_perp, int32 ind_slf)
{
	constexpr double inv_sigma2_obs_r	= 1.0/9.0;
	constexpr double inv_sigma2_obs_t 	= 1.0/100.0;
	double theta   		   				= atan2((*v_perp)(1),  (*v_perp)(0));

	M2D R = (M2D() << cos(theta), -sin(theta), sin(theta), cos(theta)).finished();
	M2D Lam_0 = (M2D() << inv_sigma2_obs_r, 0.0, 0.0,  inv_sigma2_obs_t).finished();
	
	M2D Lam_R = R * Lam_0 * R.transpose();
	V2D eta     = Lam_R * (*v_perp);
	V6D obs_msg = (V6D() << eta(0),  eta(1), Lam_R(0,0),  Lam_R(0,1), Lam_R(1,0), Lam_R(1,1)).finished();

	obs_msg = obs_msg * huber_scale_v2d(v_perp, &Lam_R);
	set_state((*pool).node, ind_slf +STS_DIM, &obs_msg); //IDX_OBS =1
	return;
}

V6D smoothness_factor(V6D* msg_v, V4D* state)
{
	constexpr double inv_sigma2_prior = 1.0 / 1.0;
	M4D Lam = (M4D() << inv_sigma2_prior, 0.0, -inv_sigma2_prior, 0.0, 0.0, inv_sigma2_prior, 0.0, -inv_sigma2_prior, -inv_sigma2_prior, 0.0, inv_sigma2_prior, 0.0, 0.0, -inv_sigma2_prior, 0.0, inv_sigma2_prior).finished(); // [4 x 4]
    V4D eta = (V4D() << 0.0, 0.0, 0.0, 0.0).finished();
	
	//  huber
    double scale   = huber_scale_v4d(state, &Lam);

	Lam = Lam * scale;
	eta = eta * scale;

	M2D Lam_aa = (M2D() << Lam(0, 0), Lam(0, 1), Lam(1, 0), Lam(1, 1)).finished();
	M2D Lam_ab = (M2D() << Lam(0, 2), Lam(0, 3), Lam(1, 2), Lam(1, 3)).finished();
	M2D Lam_bb = (M2D() << Lam(2, 2) + (*msg_v)(2), Lam(2, 3) + (*msg_v)(3), Lam(3, 2) + (*msg_v)(4), Lam(3, 3) + (*msg_v)(5)).finished();
    V2D eta_a  = eta.head(2);
    V2D eta_b  = (V2D() << eta(2, 0) + (*msg_v)(0), eta(3, 0) + (*msg_v)(1)).finished();
	
   	V2D eta_ = eta_a - Lam_ab * Lam_bb.inverse() * eta_b;
   	M2D Lam_ = Lam_aa - Lam_ab * Lam_bb.inverse() * Lam_ab.transpose();
	V6D msg = (V6D() << eta_.head(2), Lam_(0, 0),  Lam_(0, 1), Lam_(1, 0), Lam_(1, 1)).finished();
	return msg;
}


void precompute_idx(mem_pool* pool, XYTI* xyti, int32 xyi_n[4][N_EDGE], bool active[N_EDGE]){
	// #pragma simd
	// for(int32 dir=0; dir<N_EDGE; dir++){
	// 	xs[dir] = x + dirc[0][dir];
	// }
	// for(int32 dir=0; dir<N_EDGE; dir++){
	// 	ys[dir] = y + dirc[1][dir];
	// }
	// #pragma simd
	// for(int32 dir=0; dir<N_EDGE; dir++){
	// 	inds[dir] = sub2ind_(xs[dir], ys[dir], pool.H, pool.W);
	// }
	// #pragma simd
	// for(int32 dir=0; dir<N_EDGE; dir++){
	// 	active[dir] = isActive(inds[dir], t, pool.sae);
	// }

	#pragma omp simd
	for(int32 dir=0; dir<N_EDGE; dir++){
		xyi_n[0][dir] = (*xyti)[0] + dirc[0][dir];
		xyi_n[1][dir] = (*xyti)[1] + dirc[1][dir];
		xyi_n[2][dir] = sub2ind_(xyi_n[0][dir], xyi_n[1][dir], (*pool).H, (*pool).W);
		active[dir]   = isActive(xyi_n[2][dir], (*xyti).at(2), (*pool).sae);
	}
}
void send_message_Nconnect(mem_pool* pool,  XYTI* xyti)
{
    V4D state;
    V6D msg_v_n[N_EDGE];
	int32 xyi_n[4][N_EDGE];
	bool act_n[N_EDGE];

	// For SIMD 
	precompute_idx(pool, xyti, xyi_n, act_n);

	#pragma omp simd
	for(int32 dir=0; dir<N_EDGE; dir++){ 	
		 xyi_n[3][dir] = NOD_DIM*xyi_n[2][dir]  + STS_DIM;
	}

	V6D belief = update_state((*pool).node, act_n, (*xyti)[3]);

    // get state of self node
    V2D *mu = get_mu((*pool).node, (*xyti)[3]);
	state.head(2) <<*mu;

	// #pragma omp simd
    for(int32 dir=0; dir<N_EDGE; dir++){
        if (act_n[dir]){ 	// 送る方向から来るmassageを確認  activeなときはそれを引いておかなければいけない
            msg_v_n[dir] = belief - *get_state((*pool).node, (dir+2)*STS_DIM+(*xyti)[3]);
		}else{
            msg_v_n[dir] = belief;
        }
	}
	// #pragma omp simd
	for(int32 dir=0; dir<N_EDGE; dir++){
		V2D *mu_ = get_mu((*pool).node,  xyi_n[3][dir]);  //dst state
		state.tail(2) << *mu_;
        V6D msg_p       = smoothness_factor(&msg_v_n[dir], &state); // prior factor message
        set_state((*pool).node,  xyi_n[3][dir] + (dirc_idx[dir]+1)*STS_DIM, &msg_p);
    }
}

void message_passing_event(mem_pool* pool,  XYTI* xyti) // 多分再帰でかける？ k=1 hopだからいいか
{
    V6D belief;
	int32 xyi_n[4][N_EDGE];
	bool act_n[N_EDGE];

	// For SIMD 
	precompute_idx(pool, xyti, xyi_n, act_n);

	// self node
    send_message_Nconnect(pool, xyti);
	
	#pragma omp simd
    for(int32 dir=0; dir<N_EDGE; dir++){
		if (act_n[dir]){
			XYTI xyti_ = {xyi_n[0][dir], xyi_n[1][dir], (*xyti).at(2), xyi_n[2][dir]*NOD_DIM};
        	send_message_Nconnect(pool, &xyti_);
        }
    }

	// self node
	update_state((*pool).node, act_n, (*xyti)[3]);
	return;
}


void process_batch_full(mem_pool pool, int32 b_ptr)
{
	// #pragma omp for nowait
	#pragma omp for nowait schedule(dynamic)
	// #pragma omp for simd nowait schedule(dynamic)
    for(int32 i=b_ptr; i<(b_ptr+pool.WINSIZE); i++){
		V2D v_perp = (V2D() << pool.v_norms[2*i+0], pool.v_norms[2*i+1]).finished();
		int32 x = pool.indices[2*i+0]; int32 y = pool.indices[2*i+1]; 
        int32 t = pool.timestamps[i];
		int32 ind = sub2ind_(x, y,  pool.H, pool.W);

		XYTI xyti = {x, y, t, ind*NOD_DIM};
		// printf("([%3d] %3.2f, %3.2f, %2d, %2d, %3d)  thread: %d\n",i, v_perp(0), v_perp(1), x, y, t, omp_get_thread_num());
		
		pool.sae[ind] = t;
		// Compute data factor and Set the observation
		set_observation(&pool, &v_perp, xyti[3]);

		// Core of message passing
		message_passing_event(&pool, &xyti);
	}
return;
}


mem_pool initialize_full(data_cfg cfg){
	mem_pool pool;
	pool.data_name = cfg.data_name;
	pool.H = cfg.H;
	pool.W = cfg.W;
	pool.B = cfg.B;

	// memset(dirc,0, 2*(N_EDGE));
	// memset(dirc_idx,0, 1*(N_EDGE));
	constexpr int32 base_dirc[2][8] 	= {0,  0, -1, +1, -1, +1, -1, +1, -1, +1,  0,  0, -1, +1, +1, -1};  
	constexpr int32 base_dirc_idx[8]	= {1,  0,  3,  2,  5,  4,  7,  6};  

	for (int32 s_=0; s_<N_SCALE; s_++){
		int32 s = N_SCALE - s_ - 1;
		int32 offset = s*NEIGHBOR;
		for (int32 col =0; col<NEIGHBOR; col++){
			for (int32 row =0; row<2; row++){
				dirc[row][offset+col] = (int32)pow(2, s) * base_dirc[row][col];
			}
			dirc_idx[offset+col] = offset+base_dirc_idx[col]; 	// [from left message] of [right node]
		}
	}

	// printf("\n\n");
	// for (int32 row =0; row<2; row++){
	// 	for (int32 col =0; col<(N_EDGE); col++){
	// 		printf("%d,", dirc[row][col]);
	// 	}
	// }
	// printf("\n\n");

	// for (int32 col =0; col<(N_EDGE); col++){
	// 	printf("%d,", dirc_idx[col]);
	// }
	// printf("\n");

    // mem_pool pool;
    pool.node 		= (double *) malloc(NOD_DIM*pool.W*pool.H*sizeof(double));
	memset(pool.node, 0.0, NOD_DIM*pool.W*pool.H*sizeof(double));
	pool.sae  	    = (int32 *) malloc(1*pool.W*pool.H*sizeof(int32));
	memset(pool.sae, -10*DT_ACT, 1*pool.W*pool.H*sizeof(int32));
	pool.v_norms 	= (double *) malloc(2*pool.B*sizeof(double));
	pool.indices    = (int32 *) malloc(2*pool.B*sizeof(int32));
	pool.timestamps = (int32 *) malloc(1*pool.B*sizeof(int32));
	printf("initialize memory pool %d  %d  %d\n", pool.B, pool.H, pool.W);
    return pool;
}