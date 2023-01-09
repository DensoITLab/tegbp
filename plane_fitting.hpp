#ifndef __PLANE_FITTING__
#define __PLANE_FITTING__

#include <eigen3/Eigen/Dense>
using namespace Eigen;
#define BASE 3
#define WIN_X 3
#define WIN_T 40000
#define N_PATCH WIN_X*WIN_X
#define DT_ACT 50000
#define DT 5000
#define N_TH_PF 3

typedef Matrix<double, 9, 1> V9D;
typedef Matrix<double, 3, 1> V3D;
typedef Matrix<double, 2, 1> V2D;
typedef Matrix<double, 3, 3> M3D;
typedef Matrix<double, 9, 9> M9D;
typedef unsigned short uint16;
typedef short int16;
typedef unsigned int uint32;
typedef int int32;
typedef int8_t int8;
typedef uint8_t uint8;

typedef std::array<int32, 3> XYT;
typedef std::array<int32, 5> XYTPI;

struct mem_pool{
    double *BtBinv_lut; // look up table for plane fitting
    double *flow_norm;  // for checking normal flow
    double *v_norms;
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

mem_pool load_data_pf(std::string data_name);
void debug_output(mem_pool pool);
void save_data_pf(mem_pool pool, int32 seq_id, int32 index, int32 c_time);
void process_batch_pf(mem_pool pool, int32 b_ptr);
mem_pool initialize(mem_pool pool);
#endif