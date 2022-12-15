#ifndef __TEGBP__
#define __TEGBP__

#include <eigen3/Eigen/Dense>
using namespace Eigen;
#define N_SCALE 5
#define NEIGHBOR 8
#define N_EDGE NEIGHBOR*N_SCALE
#define STS_DIM 6
#define NOD_DIM STS_DIM*(2+N_EDGE)
#define DT_ACT 50000
#define WINSIZE 10000

// #define N_ITER 1

#define IDX_SLF 0
#define IDX_OBS 1

#define FAST1
#define FAST2


typedef Matrix<double, 2, 2> M2D;
typedef Matrix<double, 2, 1> V2D;
typedef Matrix<double, 4, 4> M4D;
typedef Matrix<double, 4, 1> V4D;
typedef Matrix<double, STS_DIM, 1> V6D;
typedef Matrix<double, NOD_DIM, 1> VND;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef int8_t int8;

struct mem_pool{
    double * node;
    int * sae ;
    double * v_norms;
    uint16 *indices;
    int *timestamps;
    int H;
    int W;
    int B;
};

void load_data_brick(mem_pool pool);
void load_data_dummy(mem_pool pool);
void debug_output(mem_pool pool);
void save_data(mem_pool pool, int seq_id, int index);
void process_batch(mem_pool pool, int b_ptr);
mem_pool initialize(int B, int H, int W);
int sub2ind(uint16 x, uint16 y, int node_index, int dir, int H, int W);
void init_sae(mem_pool pool);

#endif