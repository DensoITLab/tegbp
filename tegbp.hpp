#ifndef __TEGBP__
#define __TEGBP__

#include <eigen3/Eigen/Dense>
using namespace Eigen;
#define IDX_SLF 0
#define IDX_OBS 1
#define IDX_NOD 2

#define N_SCALE 5
#define NEIGHBOR 8
#define N_EDGE NEIGHBOR*N_SCALE
#define STS_DIM 6
#define NOD_DIM STS_DIM*(IDX_NOD+N_EDGE)
#define DT_ACT 50000
// #define WINSIZE 5000
// #define N_ITER 1

typedef Matrix<double, 2, 2> M2D;
typedef Matrix<double, 2, 1> V2D;
typedef Matrix<double, 4, 4> M4D;
typedef Matrix<double, 4, 1> V4D;
typedef Matrix<double, STS_DIM, 1> V6D;
typedef Matrix<double, NOD_DIM, 1> VND;
typedef unsigned short uint16;
typedef short int16;
typedef unsigned int uint32;
typedef int8_t int8;
typedef int int32;

struct mem_pool{
    double * node;
    int * sae ;
    double * v_norms;
    int16 *indices;
    int *timestamps;
    int16 H;
    int16 W;
    int B;
    std::string data_name;
    int WINSIZE;
};

mem_pool load_data(std::string data_name);
void debug_output(mem_pool pool);
void save_data(mem_pool pool, int seq_id, int index, int c_time);
void process_batch(mem_pool pool, int b_ptr);
mem_pool initialize(mem_pool pool);
int sub2ind(int16 x, int16 y, int node_index, int dir, int16 H, int16 W);
int sub2ind_obs(int16 x, int16 y, int16 H, int16 W);
int sub2ind_slf(int16 x, int16 y, int16 H, int16 W);
// void init_sae(mem_pool pool);
#endif