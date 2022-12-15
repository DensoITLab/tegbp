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

int dirc[2][2+N_EDGE]; // 2+NEIGHBOR*N_SCALE
int map_set_dirc[2+N_EDGE];

// self:0  obs:1  from up:2 from down:3  from left:4: from down:5
int sub2ind_sae(uint16 x, uint16 y, int dir, int W, int H){
	return (x+dirc[0][dir] + W*(y+dirc[1][dir]));
}
// self:0  obs:1  from up:2 from down:3  from left:4: from down:5
int sub2ind(uint16 x, uint16 y, int node_index, int dir, int W, int H){
	return NOD_DIM*(x+dirc[0][dir] + W*(y+dirc[1][dir]))+node_index*STS_DIM;
}

// Geter and Setter
void get_state(double * node, int ind, V6D* data)
{
	memcpy(data, &node[ind], sizeof(V6D));
	return;
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

V6D update_state(double * node, int* sae, uint16 x, uint16 y, int t,  int W, int H)
{
	V6D belief;
	V6D msg_from;
	int ind = sub2ind(x, y, 1, 0, W, H);
	get_state(node, ind, &belief); // ind of obs node, node_index, dir 

	for(int dir=2; dir<2+N_EDGE; dir++){
		if ((t - sae[sub2ind_sae(x, y, dir, W, H)]) < DT_ACT){
			get_state(node, sub2ind(x, y, dir+1, 0, W, H), &msg_from);
			belief = belief + msg_from;
		}
	}

	V2D mu = belief_vec_to_mu(belief);

    // ここ最適化希望
    // beliefはメッセージ送信に使う
    V6D self;
    ind = sub2ind(x, y, 0, 0, W, H); // ind of self node, node_index, dir
    get_state(node, ind, &self); //ここは値をセットする？　そのままだと値が反映されてない
	self.head(2)=mu;
	set_state(node, ind, &self);
	return belief;
}

V6D calc_data_term(V2D v_perp)
{
	double huber 		    =  1.0;
	double inv_sigma_obs_r	=  1/9.0;
	double inv_sigma_obs_t 	=  1/100.0;
	double theta   		    = std::atan(v_perp(1) / v_perp(0));

	M2D R; R << std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta);
	M2D Lam_0; Lam_0 << inv_sigma_obs_r, 0.0, 0.0,  inv_sigma_obs_t;
	M2D Lam_R = R * Lam_0 * R.transpose();
	V2D eta     = Lam_R * v_perp;
	V6D obs_msg; obs_msg << eta(0),  eta(1), Lam_R(0,0),  Lam_R(0,1), Lam_R(1,0), Lam_R(1,1);

	//  huber
	double scale   = 1.0;
	if (huber > 0){
		double M2      = (v_perp.transpose() * Lam_0 * v_perp)(0);
		double M       = std::sqrt(M2);
		if (M2 < huber*huber){
			double scale   = 1.0;
		}else{
			double scale   = 2.0*huber/M - huber*huber/M2;
		}
	}

	return obs_msg * scale;
}

V6D smoothness_factor(V6D msg_v, V4D state)
{
	double huber    		= 1.0;
    M4D Lam; Lam << 100.0, 0.0, -100.0, 0.0, 0.0, 100.0, 0.0, -100.0, -100.0, 0.0, 100.0, 0.0, 0.0, -100.0, 0.0, 100.0; // [4 x 4]
    V4D eta; eta << 0.0, 0.0, 0.0, 0.0;
    
	//  huber
    double scale   = 1.0;
	if (huber > 0){
		double M2      = (state.transpose() * Lam * state)(0);
		double M       = std::sqrt(M2);
		if (M2 < huber*huber){
			scale   = 1.0;
		}else{
			scale   = 2.0*huber/M - huber*huber/M2;
		}
	}

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

void send_message_4connect(double* node, int* sae, uint16 x, uint16 y, int t, int W, int H, V6D belief)
{
    V4D state;
    V6D self, come, node_to;
    V6D msg_v, msg_p;
    int ind, ind_to, ind_sae;

    // get state of self node
    int ind_self = sub2ind(x, y, 0, 0, W, H); // node_index, dir
    get_state(node, ind_self, &self);
    state(0) = self(0);
    state(1) = self(1);

    // variable message
    for(int dir=2; dir<2+N_EDGE; dir++){ // 4 direction
        ind_sae = sub2ind_sae(x, y, dir, W, H);     // 送る方向から来るmassageを確認
        if ((t - sae[ind_sae]) < DT_ACT){           // activeなときはそれを引いておかなければいけない
            ind = sub2ind(x, y, dir, 0, W, H);      // [dir, 0]
            get_state(node, ind, &come);
            msg_v = belief - come;
        }else{
            msg_v = belief;
        }
        // get state of destination node
        ind_to      = sub2ind(x, y, 0, dir, W, H); // [0, dir]
        get_state(node, ind_to, &node_to);
        state(2)    = node_to(0);
        state(3)    = node_to(1);
        // prior factor message
        msg_p       = smoothness_factor(msg_v, state);
        ind         = sub2ind(x, y, map_set_dirc[dir], dir, W, H); // [?, ?]
        set_state(node, ind, &msg_p);
    }
}

void message_passing_event(double* node, int* sae, uint16 x, uint16 y, int t, int W, int H) // 多分再帰でかける？ k=1 hopだからいいか
{
    int ind, ind_sae;
    int x_, y_;
    V6D belief;

    // at self node
    belief = update_state(node, sae, x, y, t, W, H);
    send_message_4connect(node, sae, x, y, t, W, H, belief);

    for(int dir=2; dir<2+N_EDGE; dir++){ 			// 4 direction
        ind_sae = sub2ind_sae(x, y, dir, W, H);    // active確認
        if ((t - sae[ind_sae]) < DT_ACT){                   // message passing at neighbor node if active
            x_ = x + dirc[0][dir];
            y_ = y + dirc[1][dir];
            // at self node
            belief = update_state(node, sae, x_, y_, t, W, H);
            send_message_4connect(node, sae, x_, y_, t, W, H, belief);
        }
    }
	update_state(node, sae, x, y, t, W, H);
	return;
}

void process_batch(mem_pool pool)
{
	#pragma omp parallel for
    for(int i=0; i<pool.B; i++){
        // double count    = 0;
        // double tot      = 0;

		V2D v_perp;
		v_perp <<  pool.v_norms[2*i+0],  pool.v_norms[2*i+1];
		uint16 x = pool.indices[2*i+0]; uint16 y = pool.indices[2*i+1]; 
        int t = pool.timestamps[i];

		// printf("([%3d] %3.2f, %3.2f, %2d, %2d, %3d)  thread: %d\n",i, v_perp(0), v_perp(1), x, y, t, omp_get_thread_num());

		// Compute data factor
		V6D obs_msg = calc_data_term(v_perp);
		// printf("(%3.2f, %3.2f) \n",obs_msg(0), obs_msg(1));

		// Set the observation
		set_state(pool.node, sub2ind(x, y, 1, 0, pool.W, pool.H), &obs_msg);

		// Core of message passing
		message_passing_event(pool.node, pool.sae, x, y, t, pool.W, pool.H);
    }
return;
}


mem_pool initialize(int B, int H, int W){
	dirc[0][0] = 0;
	dirc[1][0] = 0;
	dirc[0][1] = 0;
	dirc[1][1] = 0;
	for (int s=0; s<N_SCALE; s++){
		int offset = 2+s*NEIGHBOR;
		dirc[0][offset+0] 	= (int)pow(2, s) * 0;     // left
		dirc[1][offset+0] 	= (int)pow(2, s) * -1;
		dirc[0][offset+1] 	= (int)pow(2, s) * 0;     // right
		dirc[1][offset+1] 	= (int)pow(2, s) * 1;
		dirc[0][offset+2] 	= (int)pow(2, s) * -1;    // up
		dirc[1][offset+2] 	= (int)pow(2, s) * 0;
		dirc[0][offset+3] 	= (int)pow(2, s) * 1;     // down
		dirc[1][offset+3] 	= (int)pow(2, s) * 0;
		map_set_dirc[offset+0] 	= offset+1; 	// [from left message] of [right node]
		map_set_dirc[offset+1] 	= offset+0; 	// [from right message] of [left node]
		map_set_dirc[offset+2] 	= offset+3; 	// [from up message] of [down node]
		map_set_dirc[offset+3] 	= offset+2; 	// [from down message] of [up node]
		if (NEIGHBOR==8){
			dirc[0][offset+4] = (int)pow(2, s) * -1;     // upper-left
			dirc[1][offset+4] = (int)pow(2, s) * -1;
			dirc[0][offset+5] = (int)pow(2, s) * 1;     	// lower-right
			dirc[1][offset+5] = (int)pow(2, s) * 1;
			dirc[0][offset+6] = (int)pow(2, s) * -1;     // upper-right
			dirc[1][offset+6] = (int)pow(2, s) * 1;
			dirc[0][offset+7] = (int)pow(2, s) * 1;      // lower-left
			dirc[1][offset+7] = (int)pow(2, s) * -1;
			map_set_dirc[offset+4] = offset+5; // [from upper-left message] of [lower-right node]
			map_set_dirc[offset+5] = offset+4; // [from lower-right message] of [upper-left node]
			map_set_dirc[offset+6] = offset+7; // [from upper-right message] of [lower-left node]
			map_set_dirc[offset+7] = offset+6; // [from lower-left message] of [upper-right node]
		}
		printf("%d\n", offset);
		printf("[%d, %d]~[%d, %d]\n", dirc[0][offset+0], dirc[1][offset+0], dirc[0][offset+7], dirc[1][offset+7]);
	}

    mem_pool pool;
    pool.node 		= (double *) malloc(NOD_DIM*W*H*sizeof(double));
	memset(pool.node, 0.0, NOD_DIM*W*H*sizeof(double));
	pool.sae  	    = (int *) malloc(1*W*H*sizeof(int));
	memset(pool.sae, -10*DT_ACT, 1*W*H*sizeof(int));
	pool.v_norms 	= (double *) malloc(2*B*sizeof(double));
	pool.indices    = (uint16 *) malloc(2*B*sizeof(uint16));
	pool.timestamps = (int *) malloc(1*B*sizeof(int));
    pool.B = B;
    pool.H = H;
    pool.W = W;
    return pool;
}