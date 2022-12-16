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

int8 dirc[2][2+N_EDGE]; // 2+NEIGHBOR*N_SCALE
int8 dirc_idx[2+N_EDGE];

// self:0  obs:1  from up:2 from down:3  from left:4: from down:5
int sub2ind_sae(uint16 x, uint16 y, int dir, int H, int W){
	return (x+dirc[0][dir] + W*(y+dirc[1][dir]));
}
// self:0  obs:1  from up:2 from down:3  from left:4: from down:5
int sub2ind(uint16 x, uint16 y, int node_index, int dir, int H, int W){
	return NOD_DIM*(x+dirc[0][dir] + W*(y+dirc[1][dir]))+node_index*STS_DIM;
}

// Geter and Setter
void get_state(double * node, int ind, V6D* data)
{
	memcpy(data, &node[ind], sizeof(V6D));
	return;
}

V6D* get_state_ptr(double * node, int ind)
{
	return (V6D*)&node[ind];
}

void set_state(double * node, int ind, V6D* data)
{
	memcpy(&node[ind], data, sizeof(V6D));
	return;
}

V2D belief_vec_to_mu(V6D belief)
{
	double cond = 0.00001;
	V2D eta = belief.head(2);
	Map<M2D,Eigen::Aligned> Lam_0(belief.tail(4).data());
	Lam_0(0, 0)+=cond; // conditionig
	Lam_0(1, 1)+=cond; // conditionig
	return Lam_0.inverse() * eta;
}

V6D update_state(double * node, int* sae, uint16 x, uint16 y, int t,  int H, int W)
{
	V6D belief;
	V6D msg_from;
	int ind = sub2ind(x, y, IDX_OBS, 0, H, W);
	get_state(node, ind, &belief); // ind of obs node, node_index, dir 

	for(int dir=2; dir<2+N_EDGE; dir++){
		if ((t - sae[sub2ind_sae(x, y, dir, H, W)]) < DT_ACT){
			get_state(node, sub2ind(x, y, dir, 0, H, W), &msg_from);
			belief = belief + msg_from;
		}
	}

	V2D mu = belief_vec_to_mu(belief);

    // Set updated state
    V6D self;
    ind = sub2ind(x, y, IDX_SLF, 0, H, W); // ind of self node, node_index, dir
    get_state(node, ind, &self);
	self.head(2)=mu;
	set_state(node, ind, &self);
	return belief;
}

V6D update_state_fast(double * node, int* sae, uint16 x, uint16 y, int t,  int H, int W)
{
	V6D belief;
	V6D *msg_from;
	int ind = sub2ind(x, y, IDX_OBS, 0, H, W);
	get_state(node, ind, &belief); // ind of obs node, node_index, dir 

	for(int dir=2; dir<2+N_EDGE; dir++){
		if ((t - sae[sub2ind_sae(x, y, dir, H, W)]) < DT_ACT){
			msg_from = get_state_ptr(node, sub2ind(x, y, dir, 0, H, W));
			belief = belief + *msg_from;
		}
	}

	V2D mu = belief_vec_to_mu(belief);

    // Set updated state
    V6D *self;
    ind = sub2ind(x, y, IDX_SLF, 0, H, W); // ind of self node, node_index, dir
    self = get_state_ptr(node, ind);
	(*self).head(2)=mu;
	set_state(node, ind, self);
	return belief;
}

double huber_scale_(double M2){
	double huber	= 1.0;
	double scale	= 1.0;
	if (huber > 0){
		if (M2 < huber*huber){
			double scale   = 1.0;
		}else{
			double scale   = 2.0*huber/std::sqrt(M2) - huber*huber/M2;
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

	M2D R; R << std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta);
	M2D Lam_0; Lam_0 << inv_sigma2_obs_r, 0.0, 0.0,  inv_sigma2_obs_t;
	M2D Lam_R = R * Lam_0 * R.transpose();
	V2D eta     = Lam_R * v_perp;
	V6D obs_msg; obs_msg << eta(0),  eta(1), Lam_R(0,0),  Lam_R(0,1), Lam_R(1,0), Lam_R(1,1);

	//  huber
	return obs_msg * huber_scale_v2d(v_perp, Lam_R);
}

V6D smoothness_factor(V6D msg_v, V4D state)
{
	double huber    		= 1.0;
	double huber2    		= huber*huber;
	double inv_sigma2_prior = 1.0 / 0.25;
    M4D Lam; Lam << inv_sigma2_prior, 0.0, -inv_sigma2_prior, 0.0, 0.0, inv_sigma2_prior, 0.0, -inv_sigma2_prior, -inv_sigma2_prior, 0.0, inv_sigma2_prior, 0.0, 0.0, -inv_sigma2_prior, 0.0, inv_sigma2_prior; // [4 x 4]
    V4D eta; eta << 0.0, 0.0, 0.0, 0.0;
    
	//  huber
    double scale   = huber_scale_v4d(state, Lam);;

	Lam = Lam * scale;
	eta = eta * scale;

    M2D Lam_aa, Lam_ab, Lam_bb;
    V2D eta_a, eta_b;
    Lam_aa << Lam(0, 0), Lam(0, 1), Lam(1, 0), Lam(1, 1);
    Lam_ab << Lam(0, 2), Lam(0, 3), Lam(1, 2), Lam(1, 3);
    Lam_bb << Lam(2, 2) + msg_v(2), Lam(2, 3) + msg_v(3), Lam(3, 2) + msg_v(4), Lam(3, 3) + msg_v(5);
    eta_a << eta(0, 0), eta(1, 0);
    eta_b << eta(2, 0) + msg_v(0), eta(3, 0) + msg_v(1);

   	V2D eta_ = eta_a - Lam_ab * Lam_bb.inverse() * eta_b;
   	M2D Lam_ = Lam_aa - Lam_ab * Lam_bb.inverse() * Lam_ab.transpose();
    V6D msg;
    msg << eta_(0), eta_(1), Lam_(0, 0),  Lam_(0, 1), Lam_(1, 0), Lam_(1, 1);
	return msg;
}

void send_message_Nconnect_fast(double* node, int* sae, uint16 x, uint16 y, int t, int H, int W, V6D belief)
{
    V4D state;
    V6D *self, *come, *node_to;
    V6D msg_v, msg_p;
    int ind, ind_to, ind_sae;

    // get state of self node
    int ind_self = sub2ind(x, y, 0, 0, H, W); // node_index, dir
    self = get_state_ptr(node, ind_self);
    state(0) = (*self)(0);
    state(1) = (*self)(1);

    // variable message
    for(int dir=2; dir<2+N_EDGE; dir++){ 			// connected edge
        ind_sae = sub2ind_sae(x, y, dir, H, W);     // 送る方向から来るmassageを確認
        if ((t - sae[ind_sae]) < DT_ACT){           // activeなときはそれを引いておかなければいけない
            ind = sub2ind(x, y, dir, 0, H, W);      // [dir, 0]
            come = get_state_ptr(node, ind);
            msg_v = belief - *come;
        }else{
            msg_v = belief;
        }
        // get state of destination node
        ind_to      = sub2ind(x, y, IDX_OBS, dir, H, W); 		// [0, dir]
        node_to = get_state_ptr(node, ind_to);
        state(2)    = (*node_to)(0);
        state(3)    = (*node_to)(1);
        // prior factor message
        msg_p       = smoothness_factor(msg_v, state);
        ind         = sub2ind(x, y, dirc_idx[dir], dir, H, W); 	// [?, ?]
        set_state(node, ind, &msg_p);
    }
}
void send_message_Nconnect(double* node, int* sae, uint16 x, uint16 y, int t, int H, int W, V6D belief)
{
    V4D state;
    // V6D self, come, node_to;
	V6D self, come, node_to;
    V6D msg_v, msg_p;
    int ind, ind_to, ind_sae;

    // get state of self node
    int ind_self = sub2ind(x, y, 0, 0, H, W); // node_index, dir
    get_state(node, ind_self, &self);
    state(0) = self(0);
    state(1) = self(1);

    // variable message
    for(int dir=2; dir<2+N_EDGE; dir++){ // 4 direction
        ind_sae = sub2ind_sae(x, y, dir, H, W);     // 送る方向から来るmassageを確認
        if ((t - sae[ind_sae]) < DT_ACT){           // activeなときはそれを引いておかなければいけない
            ind = sub2ind(x, y, dir, 0, H, W);      // [dir, 0]
            get_state(node, ind, &come);
            msg_v = belief - come;
        }else{
            msg_v = belief;
        }
        // get state of destination node
        ind_to      = sub2ind(x, y, IDX_OBS, dir, H, W); // [0, dir]
        get_state(node, ind_to, &node_to);
        state(2)    = node_to(0);
        state(3)    = node_to(1);
        // prior factor message
        msg_p       = smoothness_factor(msg_v, state);
        ind         = sub2ind(x, y, dirc_idx[dir], dir, H, W); // [?, ?]
        set_state(node, ind, &msg_p);
    }
}

void message_passing_event(double* node, int* sae, uint16 x, uint16 y, int t, int H, int W) // 多分再帰でかける？ k=1 hopだからいいか
{
    int ind, ind_sae;
    int x_, y_;
    V6D belief;

    // at self node
    belief = update_state_fast(node, sae, x, y, t, H, W);
    send_message_Nconnect_fast(node, sae, x, y, t, H, W, belief);

    for(int dir=2; dir<2+N_EDGE; dir++){ 			// 4 direction
        ind_sae = sub2ind_sae(x, y, dir, H, W);    // active確認
        if ((t - sae[ind_sae]) < DT_ACT){                   // message passing at neighbor node if active
            x_ = x + dirc[0][dir];
            y_ = y + dirc[1][dir];
            // at self node
            belief = update_state_fast(node, sae, x_, y_, t, H, W);
            send_message_Nconnect_fast(node, sae, x_, y_, t, H, W, belief);
        }
    }
	update_state_fast(node, sae, x, y, t, H, W);
	return;
}


void process_batch(mem_pool pool, int b_ptr)
{
	#pragma omp parallel for
	// #pragma omp for nowait
    for(int i=b_ptr; i<(b_ptr+pool.WINSIZE); i++){
        // double count    = 0;
        // double tot      = 0;

		V2D v_perp;
		v_perp << pool.v_norms[2*i+0], pool.v_norms[2*i+1];
		uint16 x = pool.indices[2*i+0]; uint16 y = pool.indices[2*i+1]; 
        int t = pool.timestamps[i];

		// printf("([%3d] %3.2f, %3.2f, %2d, %2d, %3d)  thread: %d\n",i, v_perp(0), v_perp(1), x, y, t, omp_get_thread_num());
		
		// update SAE
		pool.sae[(pool.W*y + x)] = t;

		// Compute data factor
		V6D obs_msg = calc_data_term(v_perp);
		// printf("(%3.2f, %3.2f) \n",obs_msg(0), obs_msg(1));

		// Set the observation
		set_state(pool.node, sub2ind(x, y, IDX_OBS, 0, pool.H, pool.W), &obs_msg);

		// Core of message passing
		message_passing_event(pool.node, pool.sae, x, y, t, pool.H, pool.W);
	}
return;
}

// void init_sae(mem_pool pool){
// 	// Init SAE
//     for(int i=0; i<pool.B; i++){
// 		uint16 x = pool.indices[2*i+0]; uint16 y = pool.indices[2*i+1]; 
// 		pool.sae[(pool.W*y + x)] = pool.timestamps[i];
// 		// printf("t: %d..\n",  pool.timestamps[i]);
// 	}
// }


mem_pool initialize(mem_pool pool){

	memset(dirc,0, 2*(2+N_EDGE));
	memset(dirc_idx,0, 1*(2+N_EDGE));

	int8 base_dirc[2][8] 	= {0,  0, -1, +1, -1, +1, -1, +1, -1, +1,  0,  0, -1, +1, +1, -1};  
	int8 base_dirc_idx[8]	= {1,0,3,2,5,4,7,6};  

	for (int s_=0; s_<N_SCALE; s_++){
		int s = N_SCALE - s_ - 1;
		int offset = 2+s*NEIGHBOR;
		for (int col =0; col<NEIGHBOR; col++){
			for (int row =0; row<2; row++){
				dirc[row][offset+col] = (int)pow(2, s) * base_dirc[row][col];
			}
			dirc_idx[offset+col] 	= offset+base_dirc_idx[col]; 	// [from left message] of [right node]
		}
	}

	// for (int s_=0; s_<N_SCALE; s_++){
	// 	for (int col =0; col<NEIGHBOR; col++){
	// 		int offset = 2+s_*NEIGHBOR;
	// 		printf("[%d, %d]\n", dirc[0][offset+col], dirc[1][offset+col]);
	// 	}
	// }

    // mem_pool pool;
    pool.node 		= (double *) malloc(NOD_DIM*pool.W*pool.H*sizeof(double));
	memset(pool.node, 0.0, NOD_DIM*pool.W*pool.H*sizeof(double));
	pool.sae  	    = (int *) malloc(1*pool.W*pool.H*sizeof(int));
	memset(pool.sae, -10*DT_ACT, 1*pool.W*pool.H*sizeof(int));
	pool.v_norms 	= (double *) malloc(2*pool.B*sizeof(double));
	pool.indices    = (uint16 *) malloc(2*pool.B*sizeof(uint16));
	pool.timestamps = (int *) malloc(1*pool.B*sizeof(int));
	printf("initialize memory pool %d  %d  %d\n", pool.B, pool.H, pool.W);


    return pool;
}