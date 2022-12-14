// https://ichigopack.net/openmp/omp_base_3.html
// https://msu-cmse-courses.github.io/cmse401-S21-student/assignments/0225-HW2_OMP.html
// g++ -o process -lm main_process.cpp  -fopenmp
// ./process 64

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

// const int DT_ACT = 10;
// OMP magic
// #pragma omp critical
// #pragma omp barrier
// #pragma omp barrier

int dirc[2][6];
int map_set_dirc[6];

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
	double cond = 0.01;
	V2D eta = belief.head(2);
	Map<M2D,Eigen::Aligned> Lam_0(belief.tail(4).data());
	Lam_0(0, 0)+=cond; // conditionig
	Lam_0(1, 1)+=cond; // conditionig
	return Lam_0.inverse() * eta;
}

V6D update_state(double * node, int* sae, uint16 x, uint16 y, int t,  int W, int H)
{
	V6D belief;
	V6D belief_;
	int ind = sub2ind(x, y, 1, 0, W, H);
	get_state(node, ind, &belief); // ind of obs node, node_index, dir

	for(int dir=2;dir<6;dir++){
		if ((t - sae[sub2ind_sae(x, y, dir, W, H)]) < DT_ACT){
			get_state(node, sub2ind(x, y, dir+1, 0, W, H), &belief_);
			belief = belief + belief_;
		}
	}

	V2D mu = belief_vec_to_mu(belief);

    // ここ最適化希望
    // beliefはメッセージ送信に使う
    // V6D self;
    ind = sub2ind(x, y, 0, 0, W, H); // ind of self node, node_index, dir
    // get_state(node, ind, &self); //ここは値をセットする？　そのままだと値が反映されてない
	belief.head(2)=mu;
    // self(0) = mu(0);
    // self(1) = mu(1);

	set_state(node, ind, &belief);
	return belief;
}

V6D calc_data_term(V2D v_perp)
{
	double huber 		=  1.0;
	double sigma_obs 	=  1.0;
	double theta   		= std::atan(v_perp(1) / v_perp(0));

	M2D R; R << std::cos(theta), std::sin(theta), -std::sin(theta), std::cos(theta);
	M2D Lam_0; Lam_0 << sigma_obs, 0.0, 0.0,  sigma_obs;
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
	double huber    = 1.0;
    M4D Lam; Lam << 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0; // [4 x 4]
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
    for(int dir=2; dir<6; dir++){ // 4 direction
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

    for(int dir=2; dir<6; dir++){ // 4 direction
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
#pragma omp parallel
{
	#pragma omp for
    for(int i=0; i<pool.B; i++){
        // double count    = 0;
        // double tot      = 0;

		V2D v_perp;
		v_perp <<  pool.v_norms[2*i+0],  pool.v_norms[2*i+1];
		uint16 x = pool.indices[2*i+0]; uint16 y = pool.indices[2*i+1]; 
        int t = pool.timestamps[i] + 100;

		// printf("([%3d] %3.2f, %3.2f, %2d, %2d, %3d)  thread: %d\n",i, v_perp(0), v_perp(1), x, y, t, omp_get_thread_num());

		// Compute data factor
		V6D obs_msg = calc_data_term(v_perp);

		// Set the observation
		set_state(pool.node, sub2ind(x, y, 0, 0, pool.W, pool.H), &obs_msg);

		// Core of message passing
		message_passing_event(pool.node, pool.sae, x, y, t, pool.W, pool.H);
    }
}
return;
}


mem_pool initialize(int B, int H, int W){
	dirc[0][0] = 0;
	dirc[1][0] = 0;
	dirc[0][1] = 0;
	dirc[1][1] = 0;
	dirc[0][2] = 0;     // left
	dirc[1][2] = -1;
	dirc[0][3] = 0;     // right
	dirc[1][3] = 1;
	dirc[0][4] = -1;    // up
	dirc[1][4] = 0;
	dirc[0][5] = 1;     // down
	dirc[1][5] = 0;
	map_set_dirc[2] = 3; // [from left message] of [right node]
	map_set_dirc[3] = 2; // [from right message] of [left node]
	map_set_dirc[4] = 5; // [from up message] of [down node]
	map_set_dirc[5] = 4; // [from down message] of [up node]

    mem_pool pool;
    pool.node 		= (double *) malloc(NOD_DIM*W*H*sizeof(double));
	pool.sae  	    	= (int *) malloc(1*W*H*sizeof(int));
	pool.v_norms 	= (double *) malloc(2*B*sizeof(double));
	pool.indices     = (uint16 *) malloc(2*B*sizeof(uint16));
	pool.timestamps  	= (int *) malloc(1*B*sizeof(int));
    pool.B = B;
    pool.H = H;
    pool.W = W;
    return pool;
}