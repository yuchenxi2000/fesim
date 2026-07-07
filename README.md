# FESim

BaTiO3 铁电相变仿真，基于有效哈密顿量理论

[English](https://github.com/yuchenxi2000/fesim/blob/main/README_en.md)

## 编译

本项目使用 CMake 构建。

```bash
# 克隆（含 submodule）
git clone --recurse-submodules https://github.com/yuchenxi2000/fesim.git
cd fesim

# macOS 本地开发（自动使用 Accelerate + Homebrew OpenMP/MPI）
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

# Linux HPC（Intel MKL + OpenMPI）
module load intel openmpi
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

CMake 会自动检测 BLAS 库（MKL / Apple Accelerate / OpenBLAS）。输出文件在 `build/fesim/main`（FESim）和 `build/qmat/qmat`（QMat）。

> 如果自动检测失败，可通过 `-DBLA_VENDOR=Intel10_64lp` 或 `-DBLAS_LIBRARIES=...` 手动指定 BLAS。

以下平台已测试通过：

| 平台 | CPU | 编译器 | BLAS | MPI |
|------|-----|--------|------|-----|
| RHEL 8.5 | Xeon Gold 5320 (×2, 52C/104T) | GCC 8.5.0 / ICC 19.0.4 | MKL 2019.4 | OpenMPI 5.0.6 |
| macOS 26 | Apple M1 Pro (10C) | Apple Clang 21 | Accelerate | OpenMPI 5.0.7 |

## 使用步骤

完整仿真流程分为三步：生成相互作用矩阵 → 蒙特卡洛仿真 → 绘图。

### 1. 计算偶极子相互作用矩阵（QMat）

编辑 `qmat/dipole.ini`，然后运行：

```bash
mpirun -np 4 build/qmat/qmat qmat/dipole.ini
```

在当前目录生成二进制 cache 文件。程序使用 MPI 并行，可按核数调整 `-np`。

> **重要**：给 FESim 使用时，`ax = ay = az = 1.0`，因为 FESim 内部会乘上 BaTiO3 晶格常数（7.456 Bohr）。超胞尺寸 `Nx, Ny, Nz` 必须与 FESim 的 `sim.ini` 一致。

### 2. 蒙特卡洛仿真（FESim）

编写 `sim.ini` 配置文件，然后运行：

```bash
build/fesim/main ./sim.ini
```

不指定配置文件时默认读取 `./sim.ini`。

`sim.ini` 关键参数：

| 节 | 参数 | 说明 |
|----|------|------|
| `[sys]` | `Nx, Ny, Nz` | 超胞尺寸 |
| `[init]` | `u, v` | `"random"` 随机初始化或 .bin 文件路径 |
| `[init]` | `q` | QMat 生成的 cache 文件路径 |
| `[sim]` | `T` | 温度 (K) |
| `[sim]` | `steps` | 蒙特卡洛步数（建议 ≥ 500） |
| `[sim]` | `pressure` | 压强（Pa，负值为拉应力） |
| `[monitor]` | `eta_h, e` | 输出路径和步长间隔 |
| `[out]` | `u, v, eta_h` | 最终态输出路径 |

降温过程示例见 `fesim/run-BTO.sh`（从 350K 降到 55K）。使用前请修改脚本开头的 `CALC_DIR`、`MAIN_BIN`、`CACHE_Q` 三个变量。

> 10×10×10 超胞跑 10000 步约需一天（64 核）。想快一点可以降低 STEPS，但不要小于 500，否则无法得到有意义的结果。

### 3. 绘图

```bash
python fesim/plot.py
```

用 numpy 读取 monitor 输出，matplotlib 绘制应变–温度曲线。

> ini 解析使用 [inih][2] 库（以 git submodule 引入，`external/inih/`）。

## 原理

将有效哈密顿量展开为局域软模（local soft mode）和应变张量的多项式，展开系数由第一性原理 DFT 计算得到。选择声子模式而非直接使用原子位移有两个原因：（1）声子模式对称性高，展开式简洁；（2）根据软模理论，铁电相变主要由软声子模式驱动，可忽略高阶声子。

得到哈密顿量后，用蒙特卡洛方法（原论文 [Zhang et al., 1995][1] 的方法）进行相变仿真。

## 结果展示

![](fesim/fig.png)

10×10×10 超胞降温过程中的应变–温度曲线，展示了 C → T → O → R 的相变过程。由 `fesim/plot.py` 绘制。

![](fesim/fig_prb.png)

作为对照，上图为 [Zhang et al., 1995][1] 的结果。图中 6 个应变分量的顺序疑似与原文定义不一致：原文 A 部分定义为 e<sub>xx</sub>, e<sub>yy</sub>, e<sub>zz</sub>, e<sub>yz</sub>, e<sub>zx</sub>, e<sub>xy</sub>，而图中似乎为 e<sub>xx</sub>, e<sub>yy</sub>, e<sub>zz</sub>, e<sub>xy</sub>, e<sub>yz</sub>, e<sub>zx</sub>。

[1]: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.52.6301 "Effective Hamiltonian theory"
[2]: https://github.com/benhoyt/inih "inih"

