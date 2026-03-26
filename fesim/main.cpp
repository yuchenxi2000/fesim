//
//  main.cpp
//  BTO
//
//  BaTiO3 体系铁电性仿真程序
//  现在能仿真 Nx >= 20 的体系，不用担心内存问题
//
//  Created by yuchenxi2000 on 2021/1/30.
//  Copyright © 2021年 ChenXi Yu. All rights reserved.
//

#include <iostream>
#include <string>
#include <fstream>

// openmp
#include <omp.h>

// ini
#include "INIReader.h"

// cblas
#ifdef __APPLE__
    // ycx's mac
    #include <Accelerate/Accelerate.h>
#elif defined(__linux__) && defined(__INTEL_COMPILER)
    // PKU HPC
    #include <cmath>
    #include <mkl.h>
#elif defined(__WIN32)
    // Windows
    #include <cmath>
    #include <cblas.h>
#else
    #error "not supported"
#endif

// util
#define SQUARE(x) ((x)*(x))
#define CUBIC(x) ((x)*(x)*(x))
#ifndef MIN
    #define MIN(a, b) ((a)<(b)?(a):(b))
#endif
#ifndef MAX
    #define MAX(a, b) ((a)>(b)?(a):(b))
#endif


// ---------- begin variable def ----------

// random number generator
// no rand_r in minGW, use rand
//
struct rand_gen {
#if defined(__WIN32)
    rand_gen() {
        srand((unsigned int)time(0));
    }
#else
    unsigned int seed;
    rand_gen() {
        seed = (unsigned int)time(0);
    }
#endif
    int operator () () {
#if defined(__WIN32)
        return rand();
#else
        return rand_r(&seed);
#endif
    }
    // [-1, +1)
    double rand_double() {
        return 2.0 * (double)(*this)() / (double)RAND_MAX - 1.0;
    }
    // [0, n) && n in Z
    int rand_int(int n) {
        return (*this)() % n;
    }
    // [0, 1)
    double rand_01() {
        return (double)(*this)() / (double)RAND_MAX;
    }
};

// timer
// if openmp is enabled, time from clock() is not correct
// use omp_get_wtime() instead
//
struct timer {
#ifdef _OPENMP
    double start;
#else
    clock_t start;
#endif
    
    const char * func_str;
    timer(const char * func_str) {
#ifdef _OPENMP
        start = omp_get_wtime();
#else
        start = clock();
#endif
        this->func_str = func_str;
    }
    ~timer() {
#ifdef _OPENMP
        double s = omp_get_wtime() - start;
#else
        double s = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
#endif
        printf("[timer] %s: %lf s\n", func_str, s);
    }
};

// params
// from https://link.aps.org/doi/10.1103/PhysRevB.52.6301 (1995, Vanderbilt)
//
const struct Param {
    const double ax = 7.456;
    const double ay = 7.456;
    const double az = 7.456;
    
    const double z_star = 9.956;
    const double epsilon = 5.24;
    
    const double keppa2 = 0.0568;
    const double alpha = 0.320;
    const double gamma = -0.473;
    
    const double j1 = -0.02734;
    const double j2 = 0.04020;
    const double j3 = 0.00927;
    const double j4 = -0.00815;
    const double j5 = 0.00580;
    const double j6 = 0.00370;
    const double j7 = 0.00185;
    
    const double B11 = 4.64;
    const double B12 = 1.65;
    const double B44 = 1.85;
    
    const double B1xx = -2.18;
    const double B1yy = -0.20;
    const double B4yz = -0.08;
} param;

// physics constants
const struct PhysicsConst {
    // unit conversion: https://en.wikipedia.org/wiki/Hartree_atomic_units
    const double pressure_unit = 2.942101569713e+13;
    const double energy_unit = 4.359744722207185e-18;
    const double kB = 1.380649e-23;
} physics_const;

// system dimension
// N = Nx * Ny * Nz
int Nx, Ny, Nz, N;

// local mode amplitude
// matrix shape = [Nx, Ny, Nz]
double * ux_r;
double * uy_r;
double * uz_r;

// interaction of (0, 0, 0) and (i, j, k), including harmonic term of self energy
//
// harmonic part of local mode energy
// + short range interaction
// + long range interaction (dipole interaction)
//
// 放一起的目的是提高效率，而且二次项放一起处理更合理
// 这篇文章提出的仿真程序就是这样处理的 https://journals.aps.org/prb/abstract/10.1103/PhysRevB.82.134106
// （别看文章，文章里没有，读那个程序的源码）
//
// J0只是相互作用矩阵的第一行：固定第一个点为(0, 0, 0)
// 原因：节省内存（Nx = 20 需要 2.8G）；访存优化，不然时间都消耗在内存取数据上
//
// matrix shape = [Nx, Ny, Nz]
//
double * J0xx_r;
double * J0yy_r;
double * J0zz_r;
double * J0yz_r;
double * J0zx_r;
double * J0xy_r;

// inhomogeneous strain
double * vx_r;
double * vy_r;
double * vz_r;

// homogeneous strain
double eta_h[6];

// future use (not implemented)
// electric field
double E[3];

// configure (read from sim.ini)
struct Config {
    // simulation count
    int sim_cnt;
    
    double u_step, v_step, eta_h_step;
    
    bool monitor_eta_h, monitor_e;
    int monitor_eta_h_step, monitor_e_step;
    std::string monitor_eta_h_path, monitor_e_path;
    
    bool out_u, out_v, out_eta_h;
    std::string out_u_path, out_v_path, out_eta_h_path;
} config;

// temperature
double T;

// pressure
double pressure;

class Monitor {
public:
    bool on;
    int out_step;
    FILE * out_file;
public:
    Monitor(bool on, int out_step, const char * out_file_path) : on(on), out_step(out_step) {
        if (on) {
            out_file = fopen(out_file_path, "wb");
        }
    }
    ~Monitor() {
        if (on) {
            fclose(out_file);
        }
    }
    void monitor(int step) {}
};

double e_global();

class EnergyMonitor: public Monitor {
public:
    EnergyMonitor(bool on, int out_step, const char * out_file_path) : Monitor(on, out_step, out_file_path) {}
    void monitor(int step) {
        if (on && step % out_step == 0) {
            double step_energy = e_global();
            fwrite(&step_energy, sizeof(double), 1, this->out_file);
        }
    }
};

class EtahMonitor : public Monitor {
public:
    EtahMonitor(bool on, int out_step, const char * out_file_path) : Monitor(on, out_step, out_file_path) {}
    void monitor(int step) {
        if (on && step % out_step == 0) {
            fwrite(eta_h, sizeof(double), 6, this->out_file);
        }
    }
};

// ---------- end variable def ----------


// alloc memory for u, v, J0
void alloc() {
    ux_r = new double[N];
    uy_r = new double[N];
    uz_r = new double[N];
    vx_r = new double[N];
    vy_r = new double[N];
    vz_r = new double[N];
    J0xx_r = new double[N];
    J0yy_r = new double[N];
    J0zz_r = new double[N];
    J0yz_r = new double[N];
    J0zx_r = new double[N];
    J0xy_r = new double[N];
}

// free memory of u, v, J0
void dealloc() {
    delete [] ux_r;
    delete [] uy_r;
    delete [] uz_r;
    delete [] vx_r;
    delete [] vy_r;
    delete [] vz_r;
    delete [] J0xx_r;
    delete [] J0yy_r;
    delete [] J0zz_r;
    delete [] J0yz_r;
    delete [] J0zx_r;
    delete [] J0xy_r;
}

// beta = 1 / kB / T
double T2beta(double T) {
    return physics_const.energy_unit / physics_const.kB / T;
}

// a[i][j][k] <-> a[index]
inline int array_index(int i, int j, int k, int Nx, int Ny, int Nz) {
    return Nz * (Ny * i + j) + k;
}

// safe version of array_index
int array_index_s(int i, int j, int k, int Nx, int Ny, int Nz) {
    if (i == -1) i = Nx-1;
    if (j == -1) j = Ny-1;
    if (k == -1) k = Nz-1;
    if (i == Nx) i = 0;
    if (j == Ny) j = 0;
    if (k == Nz) k = 0;
    return array_index(i, j, k, Nx, Ny, Nz);
}

// unsafe
inline double & ux(int i, int j, int k) {
    return ux_r[array_index(i, j, k, Nx, Ny, Nz)];
}
inline double & uy(int i, int j, int k) {
    return uy_r[array_index(i, j, k, Nx, Ny, Nz)];
}
inline double & uz(int i, int j, int k) {
    return uz_r[array_index(i, j, k, Nx, Ny, Nz)];
}

// unsafe
inline double & J0xx(int i, int j, int k) {
    return J0xx_r[array_index(i, j, k, Nx, Ny, Nz)];
}
inline double & J0yy(int i, int j, int k) {
    return J0yy_r[array_index(i, j, k, Nx, Ny, Nz)];
}
inline double & J0zz(int i, int j, int k) {
    return J0zz_r[array_index(i, j, k, Nx, Ny, Nz)];
}
inline double & J0yz(int i, int j, int k) {
    return J0yz_r[array_index(i, j, k, Nx, Ny, Nz)];
}
inline double & J0zx(int i, int j, int k) {
    return J0zx_r[array_index(i, j, k, Nx, Ny, Nz)];
}
inline double & J0xy(int i, int j, int k) {
    return J0xy_r[array_index(i, j, k, Nx, Ny, Nz)];
}

// safe
inline double & J0xx_s(int i, int j, int k) {
    return J0xx_r[array_index_s(i, j, k, Nx, Ny, Nz)];
}
inline double & J0yy_s(int i, int j, int k) {
    return J0yy_r[array_index_s(i, j, k, Nx, Ny, Nz)];
}
inline double & J0zz_s(int i, int j, int k) {
    return J0zz_r[array_index_s(i, j, k, Nx, Ny, Nz)];
}
inline double & J0yz_s(int i, int j, int k) {
    return J0yz_r[array_index_s(i, j, k, Nx, Ny, Nz)];
}
inline double & J0zx_s(int i, int j, int k) {
    return J0zx_r[array_index_s(i, j, k, Nx, Ny, Nz)];
}
inline double & J0xy_s(int i, int j, int k) {
    return J0xy_r[array_index_s(i, j, k, Nx, Ny, Nz)];
}

// unsafe
inline double & vx(int i, int j, int k) {
    return vx_r[array_index(i, j, k, Nx, Ny, Nz)];
}
inline double & vy(int i, int j, int k) {
    return vy_r[array_index(i, j, k, Nx, Ny, Nz)];
}
inline double & vz(int i, int j, int k) {
    return vz_r[array_index(i, j, k, Nx, Ny, Nz)];
}
// safe
inline double & vx_s(int i, int j, int k) {
    return vx_r[array_index_s(i, j, k, Nx, Ny, Nz)];
}
inline double & vy_s(int i, int j, int k) {
    return vy_r[array_index_s(i, j, k, Nx, Ny, Nz)];
}
inline double & vz_s(int i, int j, int k) {
    return vz_r[array_index_s(i, j, k, Nx, Ny, Nz)];
}

// displacement field representation to strain tensor representation
// take average of displacements within a cell
void v2eta(int i, int j, int k, double eta[6]) {
    double dvxx, dvyy, dvzz, dvxy, dvyx, dvxz, dvzx, dvyz, dvzy;
    
    // dv_x
    dvxx = 0.25 * (
                   vx_s(i+1, j, k) - vx(i, j, k)
                   + vx_s(i+1, j+1, k) - vx_s(i, j+1, k)
                   + vx_s(i+1, j, k+1) - vx_s(i, j, k+1)
                   + vx_s(i+1, j+1, k+1) - vx_s(i, j+1, k+1)
                   );
    dvyx = 0.25 * (
                   vy_s(i+1, j, k) - vy(i, j, k)
                   + vy_s(i+1, j+1, k) - vy_s(i, j+1, k)
                   + vy_s(i+1, j, k+1) - vy_s(i, j, k+1)
                   + vy_s(i+1, j+1, k+1) - vy_s(i, j+1, k+1)
                   );
    dvzx = 0.25 * (
                   vz_s(i+1, j, k) - vz(i, j, k)
                   + vz_s(i+1, j+1, k) - vz_s(i, j+1, k)
                   + vz_s(i+1, j, k+1) - vz_s(i, j, k+1)
                   + vz_s(i+1, j+1, k+1) - vz_s(i, j+1, k+1)
                   );
    
    // dv_y
    dvyy = 0.25 * (
                   vy_s(i, j+1, k) - vy(i, j, k)
                   + vy_s(i+1, j+1, k) - vy_s(i+1, j, k)
                   + vy_s(i, j+1, k+1) - vy_s(i, j, k+1)
                   + vy_s(i+1, j+1, k+1) - vy_s(i+1, j, k+1)
                   );
    dvxy = 0.25 * (
                   vx_s(i, j+1, k) - vx(i, j, k)
                   + vx_s(i+1, j+1, k) - vx_s(i+1, j, k)
                   + vx_s(i, j+1, k+1) - vx_s(i, j, k+1)
                   + vx_s(i+1, j+1, k+1) - vx_s(i+1, j, k+1)
                   );
    dvzy = 0.25 * (
                   vz_s(i, j+1, k) - vz(i, j, k)
                   + vz_s(i+1, j+1, k) - vz_s(i+1, j, k)
                   + vz_s(i, j+1, k+1) - vz_s(i, j, k+1)
                   + vz_s(i+1, j+1, k+1) - vz_s(i+1, j, k+1)
                   );
    
    // dv_z
    dvzz = 0.25 * (
                   vz_s(i, j, k+1) - vz(i, j, k)
                   + vz_s(i+1, j, k+1) - vz_s(i+1, j, k)
                   + vz_s(i, j+1, k+1) - vz_s(i, j+1, k)
                   + vz_s(i+1, j+1, k+1) - vz_s(i+1, j+1, k)
                   );
    dvxz = 0.25 * (
                   vx_s(i, j, k+1) - vx(i, j, k)
                   + vx_s(i+1, j, k+1) - vx_s(i+1, j, k)
                   + vx_s(i, j+1, k+1) - vx_s(i, j+1, k)
                   + vx_s(i+1, j+1, k+1) - vx_s(i+1, j+1, k)
                   );
    dvyz = 0.25 * (
                   vy_s(i, j, k+1) - vy(i, j, k)
                   + vy_s(i+1, j, k+1) - vy_s(i+1, j, k)
                   + vy_s(i, j+1, k+1) - vy_s(i, j+1, k)
                   + vy_s(i+1, j+1, k+1) - vy_s(i+1, j+1, k)
                   );
    
    eta[0] = dvxx;
    eta[1] = dvyy;
    eta[2] = dvzz;
    eta[3] = dvyz + dvzy;
    eta[4] = dvzx + dvxz;
    eta[5] = dvxy + dvyx;
}

double parallel_ddot(int N, double * a, double * b) {
    double dot = 0.0;
#pragma omp parallel reduction(+:dot)
    {
        int threads = omp_get_num_threads();
        int stride = N / threads;
        int num = omp_get_thread_num();
        if (num == threads - 1) {
            dot += cblas_ddot(N - stride * num, a + num * stride, 1, b + num * stride, 1);
        } else {
            dot += cblas_ddot(stride, a + num * stride, 1, b + num * stride, 1);
        }
    }
    return dot;
}

// cblas in mkl is very fast, no need to parallel (threading cost may be larger!)
// 经测试发现mkl的cblas很快，并行的线程开销比单线程算大。。（至少对于Nx=10，更大的体系可能要并行）
double sum_ux2() {
    return cblas_ddot(N, ux_r, 1, ux_r, 1);
//    return parallel_ddot(N, ux_r, ux_r);
}
double sum_uy2() {
    return cblas_ddot(N, uy_r, 1, uy_r, 1);
//    return parallel_ddot(N, uy_r, uy_r);
}
double sum_uz2() {
    return cblas_ddot(N, uz_r, 1, uz_r, 1);
//    return parallel_ddot(N, uz_r, uz_r);
}
double sum_uyuz() {
    return cblas_ddot(N, uy_r, 1, uz_r, 1);
//    return parallel_ddot(N, uy_r, uz_r);
}
double sum_uzux() {
    return cblas_ddot(N, uz_r, 1, ux_r, 1);
//    return parallel_ddot(N, uz_r, ux_r);
}
double sum_uxuy() {
    return cblas_ddot(N, ux_r, 1, uy_r, 1);
//    return parallel_ddot(N, ux_r, uy_r);
}

void convert_q_to_Q0(double * q, double * Q0, int Nx, int Ny, int Nz, double sign_change[3]) {
    // q[i, j, k] = dipole interaction between (0, 0, 0) and (i, j, k)
    // 0 <= i < (Nx+2)/2, 0 <= j < (Ny+2)/2, 0 <= k < (Nz+2)/2
    // shape of q = [Nx_h, Ny_h, Nz_h]
    //
    // Q0[i, j, k] = dipole interaction between (0, 0, 0) and (i, j, k)
    // 0 <= i < Nx, 0 <= j < Ny, 0 <= k < Nz
    // shape of Q0 = [Nx, Ny, Nz]
    //
    int Nx_h = (Nx+2)/2;
    int Ny_h = (Ny+2)/2;
    int Nz_h = (Nz+2)/2;
    // 构造Q第一行，得到(0, 0, 0)处点和其他点的相互作用
    // 反射操作可能变号，变号规则：qxx沿任意轴反射不变号；qxy沿z轴反射不变号，沿x/y轴反射变号
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                int i2, j2, k2;
                double sign = 1.0;
                if (i < Nx_h) {
                    i2 = i;
                } else {
                    i2 = Nx - i;
                    sign *= sign_change[0];
                }
                if (j < Ny_h) {
                    j2 = j;
                } else {
                    j2 = Ny - j;
                    sign *= sign_change[1];
                }
                if (k < Nz_h) {
                    k2 = k;
                } else {
                    k2 = Nz - k;
                    sign *= sign_change[2];
                }
                int tq = array_index(i2, j2, k2, Nx_h, Ny_h, Nz_h);
                int tQ = array_index(i, j, k, Nx, Ny, Nz);
                Q0[tQ] = q[tq] * sign;
            }
        }
    }
}

// cache_q file format: version 2
//
// int version
// int header_length
// cache_q_header_params_v2 header
// matrix qxx, qyy, qzz, qyz, qzx, qxy
//
struct cache_q_header_params_v2 {
    int Nx, Ny, Nz, _align_4;  // _align_4 for alignment
    double ax, ay, az, lambda, q_tol;
    bool is_compatible() {
        return this->Nx == Nx
        && this->Ny == Ny
        && this->Nz == Nz
        && abs((this->ax - this->ay) / this->ax) <= 1e-7
        && abs((this->ax - this->az) / this->ax) <= 1e-7;
    }
};

bool init_Q0_from_cache_q(const char * file_path) {
    FILE * fin = fopen(file_path, "rb");
    if (fin == 0) {
        printf("[cache_q] failed to open q cache\n");
        return false;
    }
    int version;
    fread(&version, sizeof(int), 1, fin);
    // drop support for version 1
    if (version == 1) {
        printf("[cache_q] version 1 is deprecated. please use version 2 or above\n");
        fclose(fin);
        return false;
    }
    else if (version == 2) {
        // version 2 及以后版本都会给出文件头长度。这样只关心q时直接跳过文件头即可。
        int header_length;
        fread(&header_length, sizeof(int), 1, fin);
        // file header
        cache_q_header_params_v2 header;
        fread(&header, sizeof(cache_q_header_params_v2), 1, fin);
        if (!header.is_compatible()) {
            printf("[cache_q] incompatible q cache\n");
            fclose(fin);
            return false;
        }
        int Nx_h = (Nx + 2) / 2;
        int Ny_h = (Ny + 2) / 2;
        int Nz_h = (Nz + 2) / 2;
        int N_h = Nx_h * Ny_h * Nz_h;
        // alloc q
        double * qxx = new double[N_h];
        double * qyy = new double[N_h];
        double * qzz = new double[N_h];
        double * qyz = new double[N_h];
        double * qzx = new double[N_h];
        double * qxy = new double[N_h];
        // read q
        fread(qxx, sizeof(double), N_h, fin);
        fread(qyy, sizeof(double), N_h, fin);
        fread(qzz, sizeof(double), N_h, fin);
        fread(qyz, sizeof(double), N_h, fin);
        fread(qzx, sizeof(double), N_h, fin);
        fread(qxy, sizeof(double), N_h, fin);
        // 版本2计算的q没有乘z_star^2/epsilon（为了通用性），这里乘上
        // ax可以是任意值，只要满足ax=ay=az。使用时需要使ax=7.456
        double c1 = SQUARE(param.z_star) / param.epsilon / CUBIC(param.ax) * CUBIC(header.ax);
        cblas_dscal(N_h, c1, qxx, 1);
        cblas_dscal(N_h, c1, qyy, 1);
        cblas_dscal(N_h, c1, qzz, 1);
        cblas_dscal(N_h, c1, qyz, 1);
        cblas_dscal(N_h, c1, qzx, 1);
        cblas_dscal(N_h, c1, qxy, 1);
        // q to Q0
        convert_q_to_Q0(qxx, J0xx_r, Nx, Ny, Nz, (double [3]){1.0, 1.0, 1.0});
        convert_q_to_Q0(qyy, J0yy_r, Nx, Ny, Nz, (double [3]){1.0, 1.0, 1.0});
        convert_q_to_Q0(qzz, J0zz_r, Nx, Ny, Nz, (double [3]){1.0, 1.0, 1.0});
        convert_q_to_Q0(qyz, J0yz_r, Nx, Ny, Nz, (double [3]){1.0, -1.0, -1.0});
        convert_q_to_Q0(qzx, J0zx_r, Nx, Ny, Nz, (double [3]){-1.0, 1.0, -1.0});
        convert_q_to_Q0(qxy, J0xy_r, Nx, Ny, Nz, (double [3]){-1.0, -1.0, 1.0});
        // free q
        delete [] qxx;
        delete [] qyy;
        delete [] qzz;
        delete [] qyz;
        delete [] qzx;
        delete [] qxy;
        // close file
        fclose(fin);
    } else {
        printf("[cache_q] unknown file version: %d\n", version);
        fclose(fin);
        return false;
    }
    
    return true;
}

void init_J0(const char * cache_q_path) {
    // init dipole (long range interaction)
    if (!init_Q0_from_cache_q(cache_q_path)) {
        exit(EXIT_FAILURE);
    }
    
    // dipole x2 (There's 1/2 before J0)
    cblas_dscal(N, 2.0, J0xx_r, 1);
    cblas_dscal(N, 2.0, J0yy_r, 1);
    cblas_dscal(N, 2.0, J0zz_r, 1);
    cblas_dscal(N, 2.0, J0yz_r, 1);
    cblas_dscal(N, 2.0, J0zx_r, 1);
    cblas_dscal(N, 2.0, J0xy_r, 1);
    
    // self energy
    J0xx(0, 0, 0) += 2 * param.keppa2;
    J0yy(0, 0, 0) += 2 * param.keppa2;
    J0zz(0, 0, 0) += 2 * param.keppa2;
    
    // short range interaction
    // j1
    for (int dy = -1; dy <= 1; dy += 2)
        J0xx_s(0, dy, 0) += param.j1;
    
    for (int dz = -1; dz <= 1; dz += 2)
        J0xx_s(0, 0, dz) += param.j1;
    
    for (int dx = -1; dx <= 1; dx += 2)
        J0yy_s(dx, 0, 0) += param.j1;
    
    for (int dz = -1; dz <= 1; dz += 2)
        J0yy_s(0, 0, dz) += param.j1;
    
    for (int dx = -1; dx <= 1; dx += 2)
        J0zz_s(dx, 0, 0) += param.j1;

    for (int dy = -1; dy <= 1; dy += 2)
        J0zz_s(0, dy, 0) += param.j1;
    
    // j2
    for (int dx = -1; dx <= 1; dx += 2)
        J0xx_s(dx, 0, 0) += param.j2;
    
    for (int dy = -1; dy <= 1; dy += 2)
        J0yy_s(0, dy, 0) += param.j2;
    
    for (int dz = -1; dz <= 1; dz += 2)
        J0zz_s(0, 0, dz) += param.j2;
    
    // j3
    for (int dx = -1; dx <= 1; dx += 2)
        for (int dz = -1; dz <= 1; dz += 2)
            J0xx_s(dx, 0, dz) += param.j3;
    
    for (int dx = -1; dx <= 1; dx += 2)
        for (int dy = -1; dy <= 1; dy += 2)
            J0xx_s(dx, dy, 0) += param.j3;
    
    for (int dx = -1; dx <= 1; dx += 2)
        for (int dy = -1; dy <= 1; dy += 2)
            J0yy_s(dx, dy, 0) += param.j3;
    
    for (int dy = -1; dy <= 1; dy += 2)
        for (int dz = -1; dz <= 1; dz += 2)
            J0yy_s(0, dy, dz) += param.j3;
    
    for (int dx = -1; dx <= 1; dx += 2)
        for (int dz = -1; dz <= 1; dz += 2)
            J0zz_s(dx, 0, dz) += param.j3;
    
    for (int dy = -1; dy <= 1; dy += 2)
        for (int dz = -1; dz <= 1; dz += 2)
            J0zz_s(0, dy, dz) += param.j3;
    
    // j4
    for (int dy = -1; dy <= 1; dy += 2)
        for (int dz = -1; dz <= 1; dz += 2)
            J0xx_s(0, dy, dz) += param.j4;
    
    for (int dx = -1; dx <= 1; dx += 2)
        for (int dz = -1; dz <= 1; dz += 2)
            J0yy_s(dx, 0, dz) += param.j4;
    
    for (int dx = -1; dx <= 1; dx += 2)
        for (int dy = -1; dy <= 1; dy += 2)
            J0zz_s(dx, dy, 0) += param.j4;
    
    // j5
    for (int dy = -1; dy <= 1; dy += 2)
        for (int dz = -1; dz <= 1; dz += 2)
            J0yz_s(0, dy, dz) += param.j5 * dy * dz;
    
    for (int dx = -1; dx <= 1; dx += 2)
        for (int dz = -1; dz <= 1; dz += 2)
            J0zx_s(dx, 0, dz) += param.j5 * dx * dz;
    
    for (int dx = -1; dx <= 1; dx += 2)
        for (int dy = -1; dy <= 1; dy += 2)
            J0xy_s(dx, dy, 0) += param.j5 * dx * dy;
    
    // j6
    for (int dx = -1; dx <= 1; dx += 2)
        for (int dy = -1; dy <= 1; dy += 2)
            for (int dz = -1; dz <= 1; dz += 2) {
                J0xx_s(dx, dy, dz) += param.j6;
                J0yy_s(dx, dy, dz) += param.j6;
                J0zz_s(dx, dy, dz) += param.j6;
            }
    
    // j7
    for (int dx = -1; dx <= 1; dx += 2)
        for (int dy = -1; dy <= 1; dy += 2)
            for (int dz = -1; dz <= 1; dz += 2) {
                J0yz_s(dx, dy, dz) += param.j7 * dy * dz;
                J0zx_s(dx, dy, dz) += param.j7 * dx * dz;
                J0xy_s(dx, dy, dz) += param.j7 * dx * dy;
            }
}

// local mode energy, self energy + interaction (local)
double e_u_local(int i, int j, int k) {
    // non-harmonic term
    double p = SQUARE(ux(i, j, k)) + SQUARE(uy(i, j, k)) + SQUARE(uz(i, j, k));
    double q =
    SQUARE(uy(i, j, k) * uz(i, j, k))
    + SQUARE(uz(i, j, k) * ux(i, j, k))
    + SQUARE(ux(i, j, k) * uy(i, j, k));
    double e_u = param.alpha * SQUARE(p) + param.gamma * q;
    
    double e_u_h[3] = {0.0, 0.0, 0.0};
    
    // harmonic term
    // 默认情况，openmp只会并行外层循环（不会并行嵌套循环）
    // collapse的作用是让openmp并行内部循环
    // *但是*，经测试，加collapse后效率反而降低，因为线程开销增大
    // 所以并行外层循环就够了
//#pragma omp parallel for collapse(3) reduction(+:e_u)  // 效率低
#pragma omp parallel for reduction(+:e_u_h)
    for (int i1 = 0; i1 < Nx; i1++) {
        for (int j1 = 0; j1 < Ny; j1++) {
            for (int k1 = 0; k1 < Nz; k1++) {
                int di = i1 - i;
                int dj = j1 - j;
                int dk = k1 - k;
                if (di < 0) di += Nx;
                if (dj < 0) dj += Ny;
                if (dk < 0) dk += Nz;
                e_u_h[0] += J0xx(di, dj, dk) * ux(i1, j1, k1);
                e_u_h[1] += J0yy(di, dj, dk) * uy(i1, j1, k1);
                e_u_h[2] += J0zz(di, dj, dk) * uz(i1, j1, k1);
                
                e_u_h[1] += J0yz(di, dj, dk) * uz(i1, j1, k1);
                e_u_h[2] += J0yz(di, dj, dk) * uy(i1, j1, k1);
                
                e_u_h[2] += J0zx(di, dj, dk) * ux(i1, j1, k1);
                e_u_h[0] += J0zx(di, dj, dk) * uz(i1, j1, k1);
                
                e_u_h[0] += J0xy(di, dj, dk) * uy(i1, j1, k1);
                e_u_h[1] += J0xy(di, dj, dk) * ux(i1, j1, k1);
            }
        }
    }
    
    // 自身能量重算了，减掉
    e_u_h[0] -= 0.5 * J0xx(0, 0, 0) * ux(i, j, k);
    e_u_h[1] -= 0.5 * J0yy(0, 0, 0) * uy(i, j, k);
    e_u_h[2] -= 0.5 * J0zz(0, 0, 0) * uz(i, j, k);
    e_u_h[2] -= J0yz(0, 0, 0) * uy(i, j, k);
    e_u_h[0] -= J0zx(0, 0, 0) * uz(i, j, k);
    e_u_h[1] -= J0xy(0, 0, 0) * ux(i, j, k);
    
    e_u += ux(i, j, k) * e_u_h[0];
    e_u += uy(i, j, k) * e_u_h[1];
    e_u += uz(i, j, k) * e_u_h[2];
    
    return e_u;
}

// local mode energy, self energy + interaction (global)
// !!! this function is VERY slow !!!
// 避免频繁计算总能量。计算软模总能量复杂度O(N^2)，而每个仿真步复杂度只有O(N)
// 随着体系增大，计算总能量会花去大部分时间。
// 计算总能量一次然后把每一步的dE累加会引入舍入误差，
// 误差会随累加次数增加而积累，经测试累加1000次后精度降到5e-3，而一次sweep至少2000次（对于10x10x10体系），所以行不通。
//
double e_u_global() {
    // non-harmonic term
    double p2 = 0.0;
#pragma omp parallel for reduction(+:p2)
    for (int t = 0; t < N; t++) {
        p2 += SQUARE(SQUARE(ux_r[t]) + SQUARE(uy_r[t]) + SQUARE(uz_r[t]));
    }
    double q = 0.0;
#pragma omp parallel for reduction(+:q)
    for (int t = 0; t < N; t++) {
        q += SQUARE(uy_r[t] * uz_r[t])
        + SQUARE(uz_r[t] * ux_r[t])
        + SQUARE(ux_r[t] * uy_r[t]);
    }
    double e_u = param.alpha * p2 + param.gamma * q;
    
    // harmonic term
#pragma omp parallel for collapse(3) reduction(+:e_u)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++) {
                double tmp[3] = {0.0, 0.0, 0.0};
                for (int i1 = 0; i1 < Nx; i1++)
                    for (int j1 = 0; j1 < Ny; j1++)
                        for (int k1 = 0; k1 < Nz; k1++) {
                            int di = i1 - i;
                            int dj = j1 - j;
                            int dk = k1 - k;
                            if (di < 0) di += Nx;
                            if (dj < 0) dj += Ny;
                            if (dk < 0) dk += Nz;
                            tmp[0] += 0.5 * J0xx(di, dj, dk) * ux(i1, j1, k1);
                            tmp[1] += 0.5 * J0yy(di, dj, dk) * uy(i1, j1, k1);
                            tmp[2] += 0.5 * J0zz(di, dj, dk) * uz(i1, j1, k1);
                            
                            tmp[1] += J0yz(di, dj, dk) * uz(i1, j1, k1);
                            tmp[2] += J0zx(di, dj, dk) * ux(i1, j1, k1);
                            tmp[0] += J0xy(di, dj, dk) * uy(i1, j1, k1);
                        }
                e_u += tmp[0] * ux(i, j, k);
                e_u += tmp[1] * uy(i, j, k);
                e_u += tmp[2] * uz(i, j, k);
            }
    
    return e_u;
}

// elastic energy of inhomogeneous strain (local)
double e_elas_i_local(int i, int j, int k) {
    double e_elas_i = 0.0;
    e_elas_i += 0.25 * param.B11 * (SQUARE(vx(i, j, k) - vx_s(i+1, j, k)) + SQUARE(vx(i, j, k) - vx_s(i-1, j, k)));
    e_elas_i += 0.25 * param.B11 * (SQUARE(vy(i, j, k) - vy_s(i, j+1, k)) + SQUARE(vy(i, j, k) - vy_s(i, j-1, k)));
    e_elas_i += 0.25 * param.B11 * (SQUARE(vz(i, j, k) - vz_s(i, j, k+1)) + SQUARE(vz(i, j, k) - vz_s(i, j, k-1)));
    e_elas_i += 0.125 * param.B12 * (
                          (vx(i, j, k) - vx_s(i+1, j, k))*(vy(i, j, k) - vy_s(i, j+1, k))
                          + (vx(i, j, k) - vx_s(i+1, j, k))*(vy(i, j, k) - vy_s(i, j-1, k))
                          + (vx(i, j, k) - vx_s(i-1, j, k))*(vy(i, j, k) - vy_s(i, j+1, k))
                          + (vx(i, j, k) - vx_s(i-1, j, k))*(vy(i, j, k) - vy_s(i, j-1, k))
                          );
    e_elas_i += 0.125 * param.B12 * (
                          (vy(i, j, k) - vy_s(i, j+1, k))*(vz(i, j, k) - vz_s(i, j, k+1))
                          + (vy(i, j, k) - vy_s(i, j+1, k))*(vz(i, j, k) - vz_s(i, j, k-1))
                          + (vy(i, j, k) - vy_s(i, j-1, k))*(vz(i, j, k) - vz_s(i, j, k+1))
                          + (vy(i, j, k) - vy_s(i, j-1, k))*(vz(i, j, k) - vz_s(i, j, k-1))
                          );
    e_elas_i += 0.125 * param.B12 * (
                          (vz(i, j, k) - vz_s(i, j, k+1))*(vx(i, j, k) - vx_s(i+1, j, k))
                          + (vz(i, j, k) - vz_s(i, j, k+1))*(vx(i, j, k) - vx_s(i-1, j, k))
                          + (vz(i, j, k) - vz_s(i, j, k-1))*(vx(i, j, k) - vx_s(i+1, j, k))
                          + (vz(i, j, k) - vz_s(i, j, k-1))*(vx(i, j, k) - vx_s(i-1, j, k))
                          );
    e_elas_i += 0.125 * param.B44 * (
                          SQUARE(vx(i, j, k) - vx_s(i, j+1, k) + vy(i, j, k) - vy_s(i+1, j, k))
                          + SQUARE(vx(i, j, k) - vx_s(i, j+1, k) + vy(i, j, k) - vy_s(i-1, j, k))
                          + SQUARE(vx(i, j, k) - vx_s(i, j-1, k) + vy(i, j, k) - vy_s(i+1, j, k))
                          + SQUARE(vx(i, j, k) - vx_s(i, j-1, k) + vy(i, j, k) - vy_s(i-1, j, k))
                          );
    e_elas_i += 0.125 * param.B44 * (
                          SQUARE(vy(i, j, k) - vy_s(i, j, k+1) + vz(i, j, k) - vz_s(i, j+1, k))
                          + SQUARE(vy(i, j, k) - vy_s(i, j, k+1) + vz(i, j, k) - vz_s(i, j-1, k))
                          + SQUARE(vy(i, j, k) - vy_s(i, j, k-1) + vz(i, j, k) - vz_s(i, j+1, k))
                          + SQUARE(vy(i, j, k) - vy_s(i, j, k-1) + vz(i, j, k) - vz_s(i, j-1, k))
                          );
    e_elas_i += 0.125 * param.B44 * (
                          SQUARE(vz(i, j, k) - vz_s(i+1, j, k) + vx(i, j, k) - vx_s(i, j, k+1))
                          + SQUARE(vz(i, j, k) - vz_s(i+1, j, k) + vx(i, j, k) - vx_s(i, j, k-1))
                          + SQUARE(vz(i, j, k) - vz_s(i-1, j, k) + vx(i, j, k) - vx_s(i, j, k+1))
                          + SQUARE(vz(i, j, k) - vz_s(i-1, j, k) + vx(i, j, k) - vx_s(i, j, k-1))
                          );
    
    return e_elas_i;
}

// elastic energy of inhomogeneous strain (global)
double e_elas_i_global() {
    double e_elas_i = 0.0;
//#pragma omp parallel for collapse(3) reduction(+:e_elas_i)
#pragma omp parallel for reduction(+:e_elas_i)
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                e_elas_i += e_elas_i_local(i, j, k);
            }
        }
    }
    return e_elas_i;
}

// interaction between local mode and inhomogeneous strain (local)
double e_int_i_local(int i, int j, int k) {
    double eta[6];
    v2eta(i, j, k, eta);
    return 0.5 * param.B1xx * (
                               eta[0] * SQUARE(ux(i, j, k))
                               + eta[1] * SQUARE(uy(i, j, k))
                               + eta[2] * SQUARE(uz(i, j, k))
                               )
    + 0.5 * param.B1yy * (
                          eta[0] * SQUARE(uy(i, j, k)) + eta[0] * SQUARE(uz(i, j, k))
                          + eta[1] * SQUARE(ux(i, j, k)) + eta[1] * SQUARE(uz(i, j, k))
                          + eta[2] * SQUARE(ux(i, j, k)) + eta[2] * SQUARE(uy(i, j, k))
                          )
    + param.B4yz * (
                    eta[3] * uy(i, j, k) * uz(i, j, k)
                    + eta[4] * uz(i, j, k) * ux(i, j, k)
                    + eta[5] * ux(i, j, k) * uy(i, j, k)
                    );
}

// interaction between local mode and inhomogeneous strain (global)
double e_int_i_global() {
    double e_int_i = 0.0;
//#pragma omp parallel for collapse(3) reduction(+:e_int_i)
#pragma omp parallel for reduction(+:e_int_i)
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                e_int_i += e_int_i_local(i, j, k);
            }
        }
    }
    return e_int_i;
}

// interaction between local mode and homogeneous strain (local)
double e_int_h_local(int i, int j, int k) {
    return 0.5 * param.B1xx * (
                               eta_h[0] * SQUARE(ux(i, j, k))
                               + eta_h[1] * SQUARE(uy(i, j, k))
                               + eta_h[2] * SQUARE(uz(i, j, k))
                               )
    + 0.5 * param.B1yy * (
                          eta_h[0] * SQUARE(uy(i, j, k)) + eta_h[0] * SQUARE(uz(i, j, k))
                          + eta_h[1] * SQUARE(ux(i, j, k)) + eta_h[1] * SQUARE(uz(i, j, k))
                          + eta_h[2] * SQUARE(ux(i, j, k)) + eta_h[2] * SQUARE(uy(i, j, k))
                          )
    + param.B4yz * (
                    eta_h[3] * uy(i, j, k) * uz(i, j, k)
                    + eta_h[4] * uz(i, j, k) * ux(i, j, k)
                    + eta_h[5] * ux(i, j, k) * uy(i, j, k)
                    );
}

// interaction between local mode and homogeneous strain (global)
double e_int_h_global() {
    double ux2 = sum_ux2();
    double uy2 = sum_uy2();
    double uz2 = sum_uz2();
    double uyuz = sum_uyuz();
    double uzux = sum_uzux();
    double uxuy = sum_uxuy();
    
    return 0.5 * param.B1xx * (
                               eta_h[0] * ux2
                               + eta_h[1] * uy2
                               + eta_h[2] * uz2
                               )
    + 0.5 * param.B1yy * (
                          eta_h[0] * uy2 + eta_h[0] * uz2
                          + eta_h[1] * ux2 + eta_h[1] * uz2
                          + eta_h[2] * ux2 + eta_h[2] * uy2
                          )
    + param.B4yz * (
                    eta_h[3] * uyuz
                    + eta_h[4] * uzux
                    + eta_h[5] * uxuy
                    );
}

// elastic energy of homogeneous strain (global)
double e_elas_h_global() {
    return N * (
                0.5 * param.B11 * (SQUARE(eta_h[0]) + SQUARE(eta_h[1]) + SQUARE(eta_h[2]))
                + param.B12 * (eta_h[0] * eta_h[1] + eta_h[1] * eta_h[2] + eta_h[2] * eta_h[0])
                + 0.5 * param.B44 * (SQUARE(eta_h[3]) + SQUARE(eta_h[4]) + SQUARE(eta_h[5]))
                );
}

// H = U + pV
double e_p_global() {
    return pressure * (eta_h[0] + eta_h[1] + eta_h[2]) * (N * CUBIC(param.ax) / physics_const.pressure_unit);
}

// u改变时局域能量
double e_u_change(int i, int j, int k) {
    return e_u_local(i, j, k) + e_int_i_local(i, j, k) + e_int_h_local(i, j, k);
}

// v改变时局域能量
double e_v_change(int i, int j, int k) {
    double e_v = 0.0;
    for (int i1 = i-1; i1 <= i+1; i1++) {
        for (int j1 = j-1; j1 <= j+1; j1++) {
            for (int k1 = k-1; k1 <= k+1; k1++) {
                int i2, j2, k2;
                if (i1 == -1) i2 = Nx-1;
                else if (i1 == Nx) i2 = 0;
                else i2 = i1;
                if (j1 == -1) j2 = Ny-1;
                else if (j1 == Ny) j2 = 0;
                else j2 = j1;
                if (k1 == -1) k2 = Nz-1;
                else if (k1 == Nz) k2 = 0;
                else k2 = k1;
                e_v += e_elas_i_local(i2, j2, k2);
            }
        }
    }
    for (int i1 = i - 1; i1 <= i; i1++) {
        for (int j1 = j - 1; j1 <= j; j1++) {
            for (int k1 = k - 1; k1 <= k; k1++) {
                // 只会向下越界
                int i2 = (i1 == -1) ? Nx-1 : i1;
                int j2 = (j1 == -1) ? Ny-1 : j1;
                int k2 = (k1 == -1) ? Nz-1 : k1;
                e_v += e_int_i_local(i2, j2, k2);
            }
        }
    }
    return e_v;
}

// eta_h改变时局域能量
double e_eta_h_change() {
    return e_int_h_global() + e_elas_h_global() + e_p_global();
}

// 总能量
double e_global() {
    return e_u_global() + e_elas_i_global() + e_elas_h_global() + e_int_i_global() + e_int_h_global() + e_p_global();
}

// get file suffix, "/folder/example.txt" -> ".txt"
std::string file_suffix(const std::string & path) {
    return path.substr(path.rfind('.'), std::string::npos);
}

void init_u_0() {
    memset(ux_r, 0, N * sizeof(double));
    memset(uy_r, 0, N * sizeof(double));
    memset(uz_r, 0, N * sizeof(double));
}

void init_u_random(double u0[3], double ud[3]) {
    rand_gen rand_g;
    for (int i = 0; i < N; i++) {
        ux_r[i] = u0[0] + ud[0] * rand_g.rand_double();
        uy_r[i] = u0[1] + ud[1] * rand_g.rand_double();
        uz_r[i] = u0[2] + ud[2] * rand_g.rand_double();
    }
}

void init_u_from_file(const std::string & path) {
    if (file_suffix(path) == std::string(".bin")) {
        FILE * fin = fopen(path.c_str(), "r");
        fread(ux_r, sizeof(double), N, fin);
        fread(uy_r, sizeof(double), N, fin);
        fread(uz_r, sizeof(double), N, fin);
        fclose(fin);
    } else {
        printf("[init] read u from other file format (not .bin) is not implemented\n");
        exit(EXIT_FAILURE);
    }
}

void output_u_to_file(const std::string & path) {
    if (file_suffix(path) == std::string(".bin")) {
        FILE * fin = fopen(path.c_str(), "w");
        fwrite(ux_r, sizeof(double), N, fin);
        fwrite(uy_r, sizeof(double), N, fin);
        fwrite(uz_r, sizeof(double), N, fin);
        fclose(fin);
    } else {
        printf("[output] write u to other file format (not .bin) is not implemented\n");
        exit(EXIT_FAILURE);
    }
}

void init_v_random(double v0[3], double vd[3]) {
    rand_gen rand_g;
    for (int i = 0; i < N; i++) {
        vx_r[i] = v0[0] + vd[0] * rand_g.rand_double();
        vy_r[i] = v0[1] + vd[1] * rand_g.rand_double();
        vz_r[i] = v0[2] + vd[2] * rand_g.rand_double();
    }
}

void init_v_from_file(const std::string & path) {
    if (file_suffix(path) == std::string(".bin")) {
        FILE * fin = fopen(path.c_str(), "r");
        fread(vx_r, sizeof(double), N, fin);
        fread(vy_r, sizeof(double), N, fin);
        fread(vz_r, sizeof(double), N, fin);
        fclose(fin);
    } else {
        printf("[init] read v from other file format (not .bin) is not implemented\n");
        exit(EXIT_FAILURE);
    }
}

void output_v_to_file(const std::string & path) {
    if (file_suffix(path) == std::string(".bin")) {
        FILE * fin = fopen(path.c_str(), "w");
        fwrite(vx_r, sizeof(double), N, fin);
        fwrite(vy_r, sizeof(double), N, fin);
        fwrite(vz_r, sizeof(double), N, fin);
        fclose(fin);
    } else {
        printf("[output] write v to other file format (not .bin) is not implemented\n");
        exit(EXIT_FAILURE);
    }
}

void init_eta_h_0() {
    memset(eta_h, 0, 6 * sizeof(double));
}

void init_eta_h_from_file(const std::string & path) {
    if (file_suffix(path) == std::string(".bin")) {
        FILE * fin = fopen(path.c_str(), "r");
        fread(eta_h, sizeof(double), 6, fin);
        fclose(fin);
    } else {
        printf("[init] read eta_h from other file format (not .bin) is not implemented\n");
        exit(EXIT_FAILURE);
    }
}

void output_eta_h_to_file(const std::string & path) {
    if (file_suffix(path) == std::string(".bin")) {
        FILE * fin = fopen(path.c_str(), "w");
        fwrite(eta_h, sizeof(double), 6, fin);
        fclose(fin);
    } else {
        printf("[output] write eta_h to other file format (not .bin) is not implemented\n");
        exit(EXIT_FAILURE);
    }
}

bool read_config(const std::string & config_file_path) {
    timer timer0(__FUNCTION__);
    
    INIReader reader(config_file_path);
    int error_num = reader.ParseError();
    switch (error_num) {
        case 0:
            // no error
            break;
            
        case -1:
            printf("[ini] failed to open config file\n");
            return false;
            
        case -2:
            printf("[ini] failed to alloc memory\n");
            return false;
            
        default:
            printf("[ini] incorrect file format: find error in line %d\n", error_num);
            return false;
    }
    
    // sys
    Nx = (int)reader.GetInteger("sys", "Nx", 10);
    Ny = (int)reader.GetInteger("sys", "Ny", 10);
    Nz = (int)reader.GetInteger("sys", "Nz", 10);
    N = Nx * Ny * Nz;
    
    // alloc u, v, J0
    alloc();
    
    // init
    std::string cache_q_path = reader.Get("init", "q", "./cache_q");
    init_J0(cache_q_path.c_str());
    // u
    std::string init_u_path = reader.Get("init", "u", "random");
    if (init_u_path == "random") {
        double u0[3];
        double ud[3];
        u0[0] = reader.GetReal("init", "ux0", 0.0);
        u0[1] = reader.GetReal("init", "uy0", 0.0);
        u0[2] = reader.GetReal("init", "uz0", 0.0);
        ud[0] = reader.GetReal("init", "uxd", 0.0);
        ud[1] = reader.GetReal("init", "uyd", 0.0);
        ud[2] = reader.GetReal("init", "uzd", 0.0);
        init_u_random(u0, ud);
    } else {
        init_u_from_file(init_u_path.c_str());
    }
    // v
    std::string init_v_path = reader.Get("init", "v", "random");
    if (init_v_path == "random") {
        double v0[3];
        double vd[3];
        v0[0] = reader.GetReal("init", "vx0", 0.0);
        v0[1] = reader.GetReal("init", "vy0", 0.0);
        v0[2] = reader.GetReal("init", "vz0", 0.0);
        vd[0] = reader.GetReal("init", "vxd", 0.0);
        vd[1] = reader.GetReal("init", "vyd", 0.0);
        vd[2] = reader.GetReal("init", "vzd", 0.0);
        init_v_random(v0, vd);
    } else {
        init_v_from_file(init_v_path.c_str());
    }
    // eta_h
    std::string init_eta_h_path = reader.Get("init", "eta_h", std::string());
    if (init_eta_h_path == std::string()) {
        eta_h[0] = reader.GetReal("init", "eta_h0", 0.012);
        eta_h[1] = reader.GetReal("init", "eta_h1", 0.012);
        eta_h[2] = reader.GetReal("init", "eta_h2", 0.012);
        eta_h[3] = reader.GetReal("init", "eta_h3", 0.0);
        eta_h[4] = reader.GetReal("init", "eta_h4", 0.0);
        eta_h[5] = reader.GetReal("init", "eta_h5", 0.0);
    } else {
        init_eta_h_from_file(init_eta_h_path);
    }
    
    // sim
    T = reader.GetReal("sim", "T", 400);
    config.sim_cnt = reader.GetReal("sim", "steps", N * 100);
    pressure = reader.GetReal("sim", "pressure", 4.8e+9);
    // future use
    E[0] = reader.GetReal("sim", "Ex", 0.0);
    E[1] = reader.GetReal("sim", "Ey", 0.0);
    E[2] = reader.GetReal("sim", "Ez", 0.0);
    // update step
    config.u_step = reader.GetReal("sim", "u_step", 0.05);
    config.v_step = reader.GetReal("sim", "v_step", 0.003);
    config.eta_h_step = reader.GetReal("sim", "eta_h_step", 0.003);
    
    // monitor
    config.monitor_e_step = (int)reader.GetInteger("monitor", "e_step", 0);
    config.monitor_eta_h_step = (int)reader.GetInteger("monitor", "eta_h_step", 0);
    config.monitor_e = config.monitor_e_step > 0;
    config.monitor_eta_h = config.monitor_eta_h_step > 0;
    config.monitor_e_path = reader.Get("monitor", "e", "./monitor_e");
    config.monitor_eta_h_path = reader.Get("monitor", "eta_h", "./monitor_eta_h");
    
    // out
    config.out_u_path = reader.Get("out", "u", std::string());
    config.out_v_path = reader.Get("out", "v", std::string());
    config.out_eta_h_path = reader.Get("out", "eta_h", std::string());
    config.out_u = config.out_u_path != std::string();
    config.out_v = config.out_v_path != std::string();
    config.out_eta_h = config.out_eta_h_path != std::string();
    
    return true;
}

void sim() {
    timer timer0(__FUNCTION__);
    
    // monitor
    EnergyMonitor energy_monitor(config.monitor_e, config.monitor_e_step, config.monitor_e_path.c_str());
    EtahMonitor eta_h_monitor(config.monitor_eta_h, config.monitor_eta_h_step, config.monitor_eta_h_path.c_str());
    
    // stat reject ratio
    int cnt = 0;
    int reject_cnt = 0;
    
    double beta = T2beta(T);
    rand_gen rand_g;
    // 1995, Vanderbilt: update eta_h for 2 * Nx times every sweep
    int eta_h_cnt_per_step = 2 * MAX(MAX(Nx, Ny), Nz);
    double E0, E1, dE;
    // 我一直不明白为什么要每次“sweep”所有的格点，而且从左向右按顺序扫？那为什么不从右向左呢？这样不是破坏了随机性吗？
    // 然而看了几篇论文全部是这样扫的，于是我也只能这样扫了
    for (int s = 0; s < config.sim_cnt; s++) {
        // update u, "typewriter sweep"
        // https://journals.aps.org/prb/abstract/10.1103/PhysRevB.55.6161
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++) {
                    // backup
                    double ux_backup = ux(i, j, k);
                    double uy_backup = uy(i, j, k);
                    double uz_backup = uz(i, j, k);
                    // energy before move
                    E0 = e_u_change(i, j, k);
                    // update
                    ux(i, j, k) += config.u_step * rand_g.rand_double();
                    uy(i, j, k) += config.u_step * rand_g.rand_double();
                    uz(i, j, k) += config.u_step * rand_g.rand_double();
                    // after move
                    E1 = e_u_change(i, j, k);
                    // energy change
                    dE = E1 - E0;
                    if (rand_g.rand_01() > exp(-beta * dE)) {
                        // roll back
                        ux(i, j, k) = ux_backup;
                        uy(i, j, k) = uy_backup;
                        uz(i, j, k) = uz_backup;
                        reject_cnt++;
                        dE = 0.0;
                    }
                    cnt++;
                }
        // update v
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++) {
                    // backup
                    double vx_backup = vx(i, j, k);
                    double vy_backup = vy(i, j, k);
                    double vz_backup = vz(i, j, k);
                    // before change
                    E0 = e_v_change(i, j, k);
                    // update
                    vx(i, j, k) += config.v_step * rand_g.rand_double();
                    vy(i, j, k) += config.v_step * rand_g.rand_double();
                    vz(i, j, k) += config.v_step * rand_g.rand_double();
                    // after change
                    E1 = e_v_change(i, j, k);
                    // change
                    dE = E1 - E0;
                    if (rand_g.rand_01() > exp(-beta * dE)) {
                        // roll back
                        vx(i, j, k) = vx_backup;
                        vy(i, j, k) = vy_backup;
                        vz(i, j, k) = vz_backup;
                        reject_cnt++;
                        dE = 0.0;
                    }
                    cnt++;
                }
        // update eta_h
        for (int t = 0; t < eta_h_cnt_per_step; t++) {
            double E0, E1, dE;
            // backup
            double eta_h_backup[6];
            memcpy(eta_h_backup, eta_h, sizeof(double) * 6);
            // energy before move
            E0 = e_eta_h_change();
            // trial move
            for (int t = 0; t < 6; t++) {
                eta_h[t] += config.eta_h_step * rand_g.rand_double();
            }
            // energy after move
            E1 = e_eta_h_change();
            // change of energy
            dE = E1 - E0;
            if (rand_g.rand_01() > exp(-beta * dE)) {
                // roll back
                memcpy(eta_h, eta_h_backup, sizeof(double) * 6);
                reject_cnt++;
            }
            cnt++;
        }
        
        // monitor
        energy_monitor.monitor(s);
        eta_h_monitor.monitor(s);
    }
    
    // print reject ratio
    printf("[sim] total cnt = %d, reject cnt = %d, reject ratio = %lf\n", cnt, reject_cnt, (double)reject_cnt / (double)cnt);
}

void output() {
    if (config.out_u) {
        output_u_to_file(config.out_u_path);
    }
    if (config.out_v) {
        output_v_to_file(config.out_v_path);
    }
    if (config.out_eta_h) {
        output_eta_h_to_file(config.out_eta_h_path);
    }
}

int main(int argc, const char * argv[]) {
    timer timer0(__FUNCTION__);
    
    // print num of CPUs available
    //
    // if you use slurm, you can write script like this:
    // #SBATCH --nodes=1
    // #SBATCH --ntasks-per-node=1
    // #SBATCH --cpus-per-task=32
    //
    printf("[main] openmp processors = %d\n", omp_get_num_procs());
    
    // init from config
    if (!read_config("./sim.ini")) {
        exit(EXIT_FAILURE);
    }
    
    // run simulation
    sim();
    // output
    output();
    
    // clean up
    dealloc();
    return 0;
}
