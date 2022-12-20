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

// OMP magic
// #pragma omp critical
// #pragma omp barrier
// #pragma omp barrier
// #pragma omp for nowait

int32 dirc[2][N_EDGE];
int32 dirc_idx[N_EDGE];

// self:0  obs:1  from up:2 from down:3  from left:4: from down:5
#pragma omp declare simd
int32 sub2ind(int32 x, int32 y, int32 H, int32 W){
	return NOD_DIM*(x + W*y);
}

#pragma omp declare simd
bool isActive(int32 x, int32 y, int32 dir, int32 H, int32 W, int32 t, int32* sae)
{
	return ((t - sae[(x+dirc[0][dir] + W*(y+dirc[1][dir]))]) < DT_ACT);
}

// Geter and Setter
void get_state(double * node, int32 ind, V6D* data)
{
	memcpy(data, &node[ind], sizeof(V6D));
	return;
}
V6D* get_state_ptr(double * node, int32 ind)
{
	return (V6D*)&node[ind];
}

void set_state(double * node, int32 ind, V6D* data)
{
	memcpy(&node[ind], data, sizeof(V6D));
	return;
}

V2D belief_vec_to_mu(V6D belief)
{
	double cond = 0.00001;
	M2D C = (M2D() << cond, 0.0, 0.0,  cond).finished();
	V2D eta = belief.head(2);
	Map<M2D,Eigen::Aligned> Lam_0(belief.tail(4).data());
	return (Lam_0+C).inverse() * eta;
}

V6D update_state(double * node, bool* active, int32 ind_slf,  int32 ind_obs,  int32 ind_nod)
{
	V6D belief;
	V6D *msg_from;
	get_state(node, ind_obs, &belief); // ind of obs node

	for(int32 dir=0; dir<N_EDGE; dir++){
		if (active[dir]){
			msg_from = get_state_ptr(node, ind_nod + dir*STS_DIM);
			belief = belief + *msg_from;
		}
	}

	V2D mu = belief_vec_to_mu(belief);

    // Set updated state
    V6D *self;
    self = get_state_ptr(node, ind_slf);
	(*self).head(2)=mu;
	set_state(node, ind_slf, self);
	return belief;
}

double huber_scale_(double M2){
	double huber	= 1.0;
	double huber2   = huber*huber;
	double scale	= 1.0;
	if (huber > 0){
		if (M2 < huber2){
			double scale   = 1.0;
		}else{
			double scale   = 2.0*huber/std::sqrt(M2) - huber2/M2;
		}
	}
	return scale;
}

double huber_scale_v2d(V2D v_perp, M2D Lam){
	double M2      = (v_perp.transpose() * Lam * v_perp)(0);
	return huber_scale_(M2);
}

double huber_scale_v4d(V4D state, M4D Lam){
	double M2      = (state.transpose() * Lam * state)(0);
	return huber_scale_(M2);
}

V6D calc_data_term(V2D v_perp)
{
	double inv_sigma2_obs_r	=  1.0/9.0;
	double inv_sigma2_obs_t =  1.0/100.0;
	double theta   		    = std::atan(v_perp(1) / v_perp(0));

	M2D R = (M2D() << std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta)).finished();
	M2D Lam_0 = (M2D() << inv_sigma2_obs_r, 0.0, 0.0,  inv_sigma2_obs_t).finished();
	
	M2D Lam_R = R * Lam_0 * R.transpose();
	V2D eta     = Lam_R * v_perp;
	V6D obs_msg = (V6D() << eta(0),  eta(1), Lam_R(0,0),  Lam_R(0,1), Lam_R(1,0), Lam_R(1,1)).finished();

	//  huber
	return obs_msg * huber_scale_v2d(v_perp, Lam_R);
}

V6D smoothness_factor(V6D msg_v, V4D state)
{
	double inv_sigma2_prior = 1.0 / 0.25;
	M4D Lam = (M4D() << inv_sigma2_prior, 0.0, -inv_sigma2_prior, 0.0, 0.0, inv_sigma2_prior, 0.0, -inv_sigma2_prior, -inv_sigma2_prior, 0.0, inv_sigma2_prior, 0.0, 0.0, -inv_sigma2_prior, 0.0, inv_sigma2_prior).finished(); // [4 x 4]
    V4D eta = (V4D() << 0.0, 0.0, 0.0, 0.0).finished();
	
	//  huber
    double scale   = huber_scale_v4d(state, Lam);

	Lam = Lam * scale;
	eta = eta * scale;

	M2D Lam_aa = (M2D() << Lam(0, 0), Lam(0, 1), Lam(1, 0), Lam(1, 1)).finished();
	M2D Lam_ab = (M2D() << Lam(0, 2), Lam(0, 3), Lam(1, 2), Lam(1, 3)).finished();
	M2D Lam_bb = (M2D() << Lam(2, 2) + msg_v(2), Lam(2, 3) + msg_v(3), Lam(3, 2) + msg_v(4), Lam(3, 3) + msg_v(5)).finished();
    V2D eta_a  = (V2D() << eta(0, 0), eta(1, 0)).finished();
    V2D eta_b  = (V2D() << eta(2, 0) + msg_v(0), eta(3, 0) + msg_v(1)).finished();
	
   	V2D eta_ = eta_a - Lam_ab * Lam_bb.inverse() * eta_b;
   	M2D Lam_ = Lam_aa - Lam_ab * Lam_bb.inverse() * Lam_ab.transpose();
	V6D msg = (V6D() << eta_(0), eta_(1), Lam_(0, 0),  Lam_(0, 1), Lam_(1, 0), Lam_(1, 1)).finished();
	return msg;
}

void send_message_Nconnect(double* node, int32* sae, int32 x, int32 y, int32 t, int32 H, int32 W)
{
    V4D state;
    V6D *slf, *src, *dst;
    V6D msg_v_all[N_EDGE];
	V6D msg_p;
	int32 ind_obs_all_[N_EDGE];
	bool active[N_EDGE];

	// For SIMD 
	#pragma simd
	for(int32 dir=0; dir<N_EDGE; dir++){
		active[dir] = isActive(x, y, dir, H, W, t, sae);
	}
	int32 ind_slf 	= sub2ind(x, y,  H, W); //IDX_SLF=0
	int32 ind_obs 	= ind_slf+STS_DIM;
	int32 ind_nod 	= ind_obs+STS_DIM;

	V6D belief = update_state(node, active, ind_slf, ind_obs, ind_nod);

    // get state of self node
    slf = get_state_ptr(node, ind_slf);
    state(0) = (*slf)(0);
    state(1) = (*slf)(1);
	
	// For SIMD 
	#pragma simd
	for(int32 dir=0; dir<N_EDGE; dir++){ 	
		ind_obs_all_[dir] = sub2ind(x + dirc[0][dir], y + dirc[1][dir], H, W) + STS_DIM;
	}

    for(int32 dir=0; dir<N_EDGE; dir++){
        if (active[dir]){ 	// 送る方向から来るmassageを確認  activeなときはそれを引いておかなければいけない
            src = get_state_ptr(node, dir*STS_DIM+ind_nod);
            msg_v_all[dir] = belief - *src;
		}else{
            msg_v_all[dir] = belief;
        }
	}
	#pragma simd
	for(int32 dir=0; dir<N_EDGE; dir++){
        dst 		= get_state_ptr(node, ind_obs_all_[dir]); //dst state
        state(2)    = (*dst)(0);
        state(3)    = (*dst)(1);
        msg_p       = smoothness_factor(msg_v_all[dir], state); // prior factor message
        set_state(node, ind_obs_all_[dir]+ STS_DIM + dirc_idx[dir]*STS_DIM, &msg_p);
    }
}


void message_passing_event(double* node, int32* sae, int32 x, int32 y, int32 t, int32 H, int32 W) // 多分再帰でかける？ k=1 hopだからいいか
{
    V6D belief;
	bool active[N_EDGE];

	// For SIMD 
	#pragma simd
	for(int32 dir=0; dir<N_EDGE; dir++){
		active[dir] = isActive(x, y, dir, H, W, t, sae);
	}
	int32 ind_slf = sub2ind(x, y,  H, W); //IDX_SLF=0
	int32 ind_obs = ind_slf+STS_DIM;
	int32 ind_nod = ind_obs+STS_DIM;

	// self node
    send_message_Nconnect(node, sae, x, y, t, H, W);
	
	// neighor node
    for(int32 dir=0; dir<N_EDGE; dir++){
		if (active[dir]){
        	send_message_Nconnect(node, sae, x + dirc[0][dir], y + dirc[1][dir], t, H, W);
        }
    }

	// self node
	update_state(node, active, ind_slf, ind_obs, ind_nod);
	return;
}


void process_batch(mem_pool pool, int32 b_ptr)
{
	// #pragma omp for nowait
	#pragma omp for nowait schedule(dynamic)
	// #pragma omp parallel for schedule(dynamic)
    for(int32 i=b_ptr; i<(b_ptr+pool.WINSIZE); i++){
        // double count    = 0;
        // double tot      = 0;
		V2D v_perp = (V2D() << pool.v_norms[2*i+0], pool.v_norms[2*i+1]).finished();
		int32 x = pool.indices[2*i+0]; int32 y = pool.indices[2*i+1]; 
        int32 t = pool.timestamps[i];

		// printf("([%3d] %3.2f, %3.2f, %2d, %2d, %3d)  thread: %d\n",i, v_perp(0), v_perp(1), x, y, t, omp_get_thread_num());
		pool.sae[(pool.W*y + x)] = t;

		// Compute data factor
		V6D obs_msg = calc_data_term(v_perp);
		// printf("(%3.2f, %3.2f) \n",obs_msg(0), obs_msg(1));

		// Set the observation
		set_state(pool.node, sub2ind(x, y, pool.H, pool.W) +STS_DIM, &obs_msg); //IDX_OBS =1

		// Core of message passing
		message_passing_event(pool.node, pool.sae, x, y, t, pool.H, pool.W);
	}
return;
}


mem_pool initialize(mem_pool pool){
	// memset(dirc,0, 2*(N_EDGE));
	// memset(dirc_idx,0, 1*(N_EDGE));
	int32 base_dirc[2][8] 	= {0,  0, -1, +1, -1, +1, -1, +1, -1, +1,  0,  0, -1, +1, +1, -1};  
	int32 base_dirc_idx[8]	= {1,  0,  3,  2,  5,  4,  7,  6};  

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

	printf("\n\n");
	for (int32 row =0; row<2; row++){
		for (int32 col =0; col<(N_EDGE); col++){
			printf("%d,", dirc[row][col]);
		}
	}
	printf("\n\n");

	for (int32 col =0; col<(N_EDGE); col++){
		printf("%d,", dirc_idx[col]);
	}
	printf("\n");

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