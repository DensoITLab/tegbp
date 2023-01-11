#ifndef __TEGBP__
#define __TEGBP__

#include <eigen3/Eigen/Dense>
#define DEBUG 0


using namespace Eigen;
#define N_SCALE 5
#define NEIGHBOR 8
#define N_EDGE NEIGHBOR*N_SCALE
#define STS_DIM 6
#define NOD_DIM STS_DIM*(2+N_EDGE)
#define DT_ACT 50000

#define BASE 3
#define WIN_X 3
#define WIN_T 40000
#define N_PATCH WIN_X*WIN_X
#define DT_ACT 50000
#define DT 5000
#define N_TH_PF 3

typedef Matrix<double, 2, 2> M2D;
typedef Matrix<double, 2, 1> V2D;
typedef Matrix<double, 3, 3> M3D;
typedef Matrix<double, 3, 1> V3D;
typedef Matrix<double, 4, 4> M4D;
typedef Matrix<double, 4, 1> V4D;
typedef Matrix<double, 9, 9> M9D;
typedef Matrix<double, 9, 1> V9D;
typedef Matrix<double, STS_DIM, 1> V6D;
typedef Matrix<double, NOD_DIM, 1> VND;
typedef unsigned short uint16;
typedef short int16;
typedef unsigned int uint32;
typedef int int32;
typedef int8_t int8;
typedef uint8_t uint8;

typedef std::array<int32, 3> XYT;
typedef std::array<int32, 4> XYTI;
typedef std::array<int32, 5> XYTPI;

struct mem_pool{
    double *BtBinv_lut; // look up table for plane fitting
    double *flow_norm;  // for checking normal flow
    double *node;
    double *v_norms;
    double *v_fulls;
    int32 *sae;
    int32 *sae_pol;
    int32 *indices;
    int32 *timestamps;
    int32 *polarities;    // int32も要らない
    int32 H;
    int32 W;
    int32 B;
    std::string data_name;
    int32 WINSIZE;
};

struct data_cfg{
    int32 H;
    int32 W;
    int32 B;
    std::string data_name;
    int32 WINSIZE;
};

data_cfg load_cfg(std::string data_name);
mem_pool load_data(mem_pool pool);
void debug_output(mem_pool pool);
void save_data(mem_pool pool, int32 seq_id, int32 index, int32 c_time);
void process_batch_full(mem_pool pool, int32 b_ptr);
void process_batch_nrml(mem_pool pool, int32 b_ptr);
mem_pool initialize_full(data_cfg cfg);
mem_pool initialize_nrml(data_cfg cfg);
int32 sub2ind_nod(int32 x, int32 y, int32 H, int32 W);
int32 sub2ind_(int32 x, int32 y, int32 H, int32 W);
#endif