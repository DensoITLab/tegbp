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
    double *node;
    double *v_norms;
    int32 *sae;
    int32 *indices;
    int32 *timestamps;
    int32 H;
    int32 W;
    int32 B;
    std::string data_name;
    int32 WINSIZE;
};

mem_pool load_data(std::string data_name);
void debug_output(mem_pool pool);
void save_data(mem_pool pool, int32 seq_id, int32 index, int32 c_time);
void process_batch(mem_pool pool, int32 b_ptr);
mem_pool initialize(mem_pool pool);
// int32 sub2ind_obs(int32 x, int32 y, int32 H, int32 W);
int32 sub2ind(int32 x, int32 y, int32 H, int32 W);
#endif