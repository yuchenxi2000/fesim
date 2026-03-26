// 计算dipole相互作用矩阵
// 目前理论支持所有长方体原胞
// TODO: 测试lambda默认值的收敛性
// 
// Created by yuchenxi2000 on 2021/1/19.
// Copyright © 2021年 ChenXi Yu. Licensed under MIT License.

#include <iostream>
#include <cmath>
#include <cstring>
#include <string>

// cblas不同平台的实现
#if defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#elif defined(__INTEL_COMPILER) || defined(__ICC)
#include <mkl.h>
#else
#include <cblas.h>
#endif

#define SQUARE(x) ((x)*(x))
#define CUBIC(x) ((x)*(x)*(x))

#if defined(__INTEL_COMPILER)
#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

// mpi
#include <mpi.h>
int world_size;
int world_rank;

// 配置文件
#include "INIReader.h"

// constant
double ax;
double ay;
double az;

double q_sum_tol;
double lambda;
int Nx, Ny, Nz;
int N;

// 由于对称性只需要算大约1/8的点
// 考虑R位移一个晶格矢量A(a, 0, 0)后，因为dot(G, A)=0，所以不变
// 又考虑R反射变成-R后因为会取余弦所以也不变
int Nx_h, Ny_h, Nz_h;
int N_h;

struct timer {
    clock_t start;
    const char * func_str;
    timer(const char * func_str) {
        start = clock();
        this->func_str = func_str;
    }
    ~timer() {
        double s = (double)(clock() - start) / (double)CLOCKS_PER_SEC;
        printf("[proc %d] [timer] %s: %lf s\n", world_rank, func_str, s);
    }
};

double * Rx_r;
double * Ry_r;
double * Rz_r;

inline int to_arr_index(int i, int j, int k) {
    return Nz_h * (Ny_h * i + j) + k;
}

inline double & Rx(int i, int j, int k) {
    return Rx_r[to_arr_index(i, j, k)];
}

inline double & Ry(int i, int j, int k) {
    return Ry_r[to_arr_index(i, j, k)];
}

inline double & Rz(int i, int j, int k) {
    return Rz_r[to_arr_index(i, j, k)];
}

void init_R() {
    for (int i = 0; i < Nx_h; i++) {
        for (int j = 0; j < Ny_h; j++) {
            for (int k = 0; k < Nz_h; k++) {
                Rx(i, j, k) = ax*i;
                Ry(i, j, k) = ay*j;
                Rz(i, j, k) = az*k;
            }
        }
    }
}

#if defined(__APPLE__)
// 函数原型和mkl一致，应对本地mac没有mkl的情况，因为本地只做开发所以就用一个效率低的
void vdCos(int n, double * a, double * y) {
    for (int i = 0; i < n; i++) {
        y[i] = cos(a[i]);
    }
}
#elif defined(__INTEL_COMPILER) || defined(__ICC)
// 用 mkl.h 里的 vdCos
#else
#error "not supported"
#endif

double * qxx_r;
double * qyy_r;
double * qzz_r;
double * qyz_r;
double * qzx_r;
double * qxy_r;

inline double & qxx(int i, int j, int k) {
    return qxx_r[to_arr_index(i, j, k)];
}

inline double & qyy(int i, int j, int k) {
    return qyy_r[to_arr_index(i, j, k)];
}

inline double & qzz(int i, int j, int k) {
    return qzz_r[to_arr_index(i, j, k)];
}

inline double & qyz(int i, int j, int k) {
    return qyz_r[to_arr_index(i, j, k)];
}

inline double & qzx(int i, int j, int k) {
    return qzx_r[to_arr_index(i, j, k)];
}

inline double & qxy(int i, int j, int k) {
    return qxy_r[to_arr_index(i, j, k)];
}

double * tmp1;

double * qxx_r_n;
double * qyy_r_n;
double * qzz_r_n;
double * qyz_r_n;
double * qzx_r_n;
double * qxy_r_n;

// 给定G算求和中的一项
void sum_element_q(double Gx, double Gy, double Gz) {
    double G2 = SQUARE(Gx) + SQUARE(Gy) + SQUARE(Gz);
    memset(tmp1, 0, N_h * sizeof(double));  // 这个虽然是double，但是可以这样置0.0
    cblas_daxpy(N_h, Gx, Rx_r, 1, tmp1, 1);
    cblas_daxpy(N_h, Gy, Ry_r, 1, tmp1, 1);
    cblas_daxpy(N_h, Gz, Rz_r, 1, tmp1, 1);
    vdCos(N_h, tmp1, tmp1);
    double c1 = exp(-G2/SQUARE(2 * lambda))/G2;
    cblas_daxpy(N_h, Gx*Gx*c1, tmp1, 1, qxx_r_n, 1);
    cblas_daxpy(N_h, Gy*Gy*c1, tmp1, 1, qyy_r_n, 1);
    cblas_daxpy(N_h, Gz*Gz*c1, tmp1, 1, qzz_r_n, 1);
    cblas_daxpy(N_h, Gy*Gz*c1, tmp1, 1, qyz_r_n, 1);
    cblas_daxpy(N_h, Gz*Gx*c1, tmp1, 1, qzx_r_n, 1);
    cblas_daxpy(N_h, Gx*Gy*c1, tmp1, 1, qxy_r_n, 1);
}

int file_format_version = 2;

struct cache_q_header_params_v2 {
    int Nx, Ny, Nz, _align_4;
    double ax, ay, az, lambda, q_tol;
    bool is_compatible() {
        return true;
    }
};

void save_q_to_cache(const char * file_path) {
    FILE * fcache = fopen(file_path, "wb");
    // int version
    fwrite(&file_format_version, sizeof(int), 1, fcache);
    // int header_length
    int header_length = sizeof(cache_q_header_params_v2);
    fwrite(&header_length, sizeof(int), 1, fcache);
    // cache_q_header_params_v2 info （文件头）
    cache_q_header_params_v2 info;
    info.Nx = Nx;
    info.Ny = Ny;
    info.Nz = Nz;
    info.ax = ax;
    info.ay = ay;
    info.az = az;
    info.lambda = lambda;
    info.q_tol = q_sum_tol;
    fwrite(&info, sizeof(cache_q_header_params_v2), 1, fcache);
    // q mat
    fwrite(qxx_r, sizeof(double), N_h, fcache);
    fwrite(qyy_r, sizeof(double), N_h, fcache);
    fwrite(qzz_r, sizeof(double), N_h, fcache);
    fwrite(qyz_r, sizeof(double), N_h, fcache);
    fwrite(qzx_r, sizeof(double), N_h, fcache);
    fwrite(qxy_r, sizeof(double), N_h, fcache);
    fclose(fcache);
}

// 第n个“壳层”的和，只需要算一半的点，后面乘二即可
void sum_n(int n) {
    // 置零
    memset(qxx_r_n, 0, N_h * sizeof(double));
    memset(qyy_r_n, 0, N_h * sizeof(double));
    memset(qzz_r_n, 0, N_h * sizeof(double));
    memset(qyz_r_n, 0, N_h * sizeof(double));
    memset(qzx_r_n, 0, N_h * sizeof(double));
    memset(qxy_r_n, 0, N_h * sizeof(double));
    
    int i, j, k;
    double c[3] = {
        2.0*M_PI/Nx/ax,
        2.0*M_PI/Ny/ay,
        2.0*M_PI/Nz/az
    };
    int cnt = 0;
    j = -n;
    for (i = -n; i <= n; i++) {
        for (k = -n+1; k < n; k++) {
            if (cnt % world_size == world_rank)
                sum_element_q(c[0]*i, c[1]*j, c[2]*k);
            cnt++;
        }
    }
    k = -n;
    for (i = -n+1; i < n; i++) {
        for (j = -n; j <= n; j++) {
            if (cnt % world_size == world_rank)
                sum_element_q(c[0]*i, c[1]*j, c[2]*k);
            cnt++;
        }
    }
    i = -n;
    for (j = -n+1; j < n; j++) {
        for (k = -n; k <= n; k++) {
            if (cnt % world_size == world_rank)
                sum_element_q(c[0]*i, c[1]*j, c[2]*k);
            cnt++;
        }
    }
    if (cnt % world_size == world_rank)
        sum_element_q(c[0]*(-n), c[1]*(-n), c[2]*(-n));
    cnt++;
    if (cnt % world_size == world_rank)
        sum_element_q(c[0]*(n), c[1]*(-n), c[2]*(-n));
    cnt++;
    if (cnt % world_size == world_rank)
        sum_element_q(c[0]*(-n), c[1]*(n), c[2]*(-n));
    cnt++;
    if (cnt % world_size == world_rank)
        sum_element_q(c[0]*(-n), c[1]*(-n), c[2]*(n));
//    cnt++;
}

#if defined(__APPLE__)
// 本地开发测试用，垃圾
void vdLinearFrac(int n, double * a, double * b, double scalea, double shifta, double scaleb, double shiftb, double * y) {
    for (int i = 0; i < n; i++) {
        y[i] = (scalea * a[i] + shifta) / (scaleb * b[i] + shiftb);
    }
}
#elif defined(__INTEL_COMPILER)
// 用 mkl.h 里的 vdLinearFrac
#else
#error "not supported"
#endif

double * tmp2;

bool q_converged(double * qn, double * q, double tol) {
    vdLinearFrac(N_h, qn, q, 1.0, 0.0, 1.0, 0.0, tmp2);
    double err = fabs(tmp2[cblas_idamax(N_h, tmp2, 1)]);
    printf("[proc %d] [init_q] err = %lf\n", world_rank, err);
    return err <= tol;
}

void init_q() {
    timer timer0(__FUNCTION__);
    
    // 置零
    memset(qxx_r, 0, N_h * sizeof(double));
    memset(qyy_r, 0, N_h * sizeof(double));
    memset(qzz_r, 0, N_h * sizeof(double));
    memset(qyz_r, 0, N_h * sizeof(double));
    memset(qzx_r, 0, N_h * sizeof(double));
    memset(qxy_r, 0, N_h * sizeof(double));
    
    double * qxx_r_n_g = new double[N_h];
    double * qyy_r_n_g = new double[N_h];
    double * qzz_r_n_g = new double[N_h];
    double * qyz_r_n_g = new double[N_h];
    double * qzx_r_n_g = new double[N_h];
    double * qxy_r_n_g = new double[N_h];
    
    int n = 1;
    int stop = 0;
    while (!stop) {
        sum_n(n);
        MPI_Reduce(qxx_r_n, qxx_r_n_g, N_h, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(qyy_r_n, qyy_r_n_g, N_h, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(qzz_r_n, qzz_r_n_g, N_h, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(qyz_r_n, qyz_r_n_g, N_h, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(qzx_r_n, qzx_r_n_g, N_h, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(qxy_r_n, qxy_r_n_g, N_h, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (world_rank == 0) {
            stop = n > 5;
            stop = stop && q_converged(qxx_r_n_g, qxx_r, q_sum_tol);
            stop = stop && q_converged(qyy_r_n_g, qyy_r, q_sum_tol);
            stop = stop && q_converged(qzz_r_n_g, qzz_r, q_sum_tol);
            stop = stop && q_converged(qyz_r_n_g, qyz_r, q_sum_tol);
            stop = stop && q_converged(qzx_r_n_g, qzx_r, q_sum_tol);
            stop = stop && q_converged(qxy_r_n_g, qxy_r, q_sum_tol);
            // 前面由于空间反射对称，每项只算了一半，所以这里乘二
            cblas_daxpy(N_h, 2.0, qxx_r_n_g, 1, qxx_r, 1);
            cblas_daxpy(N_h, 2.0, qyy_r_n_g, 1, qyy_r, 1);
            cblas_daxpy(N_h, 2.0, qzz_r_n_g, 1, qzz_r, 1);
            cblas_daxpy(N_h, 2.0, qyz_r_n_g, 1, qyz_r, 1);
            cblas_daxpy(N_h, 2.0, qzx_r_n_g, 1, qzx_r, 1);
            cblas_daxpy(N_h, 2.0, qxy_r_n_g, 1, qxy_r, 1);
        }
        MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
        n++;
    }

    if (world_rank == 0) {
        printf("[proc 0] [init_q] sum converged, n = %d\n", n);
        // 现在q是sigma求和后的式子
        // 乘上常数、加上第二项
        double c1 = 2.0*M_PI/(ax*ay*az*N);
        cblas_dscal(N_h, c1, qxx_r, 1);
        cblas_dscal(N_h, c1, qyy_r, 1);
        cblas_dscal(N_h, c1, qzz_r, 1);
        cblas_dscal(N_h, c1, qyz_r, 1);
        cblas_dscal(N_h, c1, qzx_r, 1);
        cblas_dscal(N_h, c1, qxy_r, 1);
        double c2 = 2.0*CUBIC(lambda)/3.0/sqrt(M_PI);
        qxx_r[0] -= c2;
        qyy_r[0] -= c2;
        qzz_r[0] -= c2;
    }
    
    delete [] qxx_r_n_g;
    delete [] qyy_r_n_g;
    delete [] qzz_r_n_g;
    delete [] qyz_r_n_g;
    delete [] qzx_r_n_g;
    delete [] qxy_r_n_g;
    
}

void alloc_init() {
    Rx_r = new double[N_h];
    Ry_r = new double[N_h];
    Rz_r = new double[N_h];
    init_R();
    tmp1 = new double[N_h];
    tmp2 = new double[N_h];
    qxx_r = new double[N_h];
    qyy_r = new double[N_h];
    qzz_r = new double[N_h];
    qyz_r = new double[N_h];
    qzx_r = new double[N_h];
    qxy_r = new double[N_h];
    qxx_r_n = new double[N_h];
    qyy_r_n = new double[N_h];
    qzz_r_n = new double[N_h];
    qyz_r_n = new double[N_h];
    qzx_r_n = new double[N_h];
    qxy_r_n = new double[N_h];
    init_q();
}

void dealloc() {
    delete [] Rx_r;
    delete [] Ry_r;
    delete [] Rz_r;
    delete [] tmp1;
    delete [] tmp2;
    delete [] qxx_r;
    delete [] qyy_r;
    delete [] qzz_r;
    delete [] qyz_r;
    delete [] qzx_r;
    delete [] qxy_r;
    delete [] qxx_r_n;
    delete [] qyy_r_n;
    delete [] qzz_r_n;
    delete [] qyz_r_n;
    delete [] qzx_r_n;
    delete [] qxy_r_n;
}

int main(int argc, const char * argv[]) {
    
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    timer timer0(__FUNCTION__);
    
    std::string save_path;
    
    int config_no_error = 1;
    
    if (world_rank == 0) {
        INIReader reader("./dipole.ini");
        if (reader.ParseError() != 0) {
            printf("[ini] incorrect config file\n");
            config_no_error = 0;
        } else {
            // 超胞
            Nx = (int)reader.GetInteger("sys", "Nx", 5);
            Ny = (int)reader.GetInteger("sys", "Ny", 5);
            Nz = (int)reader.GetInteger("sys", "Nz", 5);
            if (Nx <= 0 || Ny <= 0 || Nz <= 0) {
                printf("[ini] system dim must be positive\n");
                config_no_error = 0;
            } else if (Nx > 80 || Ny > 80 || Nz > 80) {
                printf("[ini] warning: system too large\n");
            }
            N = Nx * Ny * Nz;
            Nx_h = (Nx+2)/2;
            Ny_h = (Ny+2)/2;
            Nz_h = (Nz+2)/2;
            N_h = Nx_h * Ny_h * Nz_h;
            // 晶格常数
            ax = reader.GetReal("lattice", "ax", 7.456);
            ay = reader.GetReal("lattice", "ay", 7.456);
            az = reader.GetReal("lattice", "az", 7.456);
            if (ax <= 0.0 || ay <= 0.0 || az <= 0.0) {
                printf("[ini] lattice constants must be positive\n");
                config_no_error = 0;
            }
            // lambda取值来自 https://gitlab.com/pyseries/PyEwald ，需要测试
            q_sum_tol = reader.GetReal("dipole", "tol", 1e-6);
            if (q_sum_tol <= 0.0) {
                printf("[ini] q_sum_tol must be positive\n");
                config_no_error = 0;
            } else if (q_sum_tol > 1e-3) {
                printf("[ini] q_sum_tol too large\n");
                config_no_error = 0;
            }
            lambda = reader.GetReal("dipole", "lambda", sqrt(-log(q_sum_tol))/MIN(MIN(ax, ay), az));
            save_path = reader.Get("output", "q_file_path", "./cache_q");
        }
    }
    
    MPI_Bcast(&config_no_error, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (!config_no_error) {
        MPI_Finalize();
        return 0;
    }
    
    MPI_Bcast(&Nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Nz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Nx_h, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Ny_h, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Nz_h, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N_h, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ay, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&az, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&q_sum_tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lambda, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    alloc_init();
    if (world_rank == 0)
        save_q_to_cache(save_path.c_str());
    dealloc();
    
    MPI_Finalize();
    
    return 0;
}

