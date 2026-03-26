# FESim

BaTiO3铁电仿真，基于有效哈密顿量理论

[English](https://github.com/yuchenxi2000/fesim/blob/main/README_en.md)

## 使用方法

`fesim`目录下是仿真主程序FESim，实现BaTiO3超胞的蒙特卡洛仿真。体系哈密顿量使用[Zhang et al., 1995][1]中的有效哈密顿量。该程序使用CBLAS库实现矩阵计算，OpenMP并行计算哈密顿量。因为不同机器的CBLAS库不同，加上该程序比较简单（只有两个文件），因此这里暂时不给统一的Makefile，可以参考`make.sh`进行编译（针对自己用的超算）。

FESim的配置文件是ini格式，示例运行脚本见`fesim/run-BTO.sh`。该脚本实现了降温过程中铁电相变的仿真，降温的实现是写一个配置文件`sim.ini`并运行仿真程序，固定温度下跑一定仿真步数后，再写一个不同温度的配置文件，继续这一过程。

> `fesim/run-BTO.sh`需要跑一天左右（在一个64核机器上）。想快一点可以降低STEPS，但不要小于500，否则无法得到有意义的结果。

`qmat`目录下是计算偶极子相互作用矩阵的程序QMat。QMat计算FESim每个蒙特卡洛步计算能量时需要的偶极子相互作用矩阵，输出一个二进制文件，这个文件在FESim运行时被读取。在有效哈密顿量理论中，相互作用矩阵近似不变，因此提前计算能够减小计算量。该程序依赖CBLAS库和OpenMPI运行库，不同机器的配置不同，需要根据具体配置进行编译，因此目录下给了一些参考的编译脚本。

QMat配置文件是ini格式，示例输入文件见`qmat_v2/dipole.ini`。需要注意，如果想生成FESim所需的相互作用矩阵，晶格常数ax, ay, az必须设置成1.0，因为FESim会自动考虑BaTiO3的晶格常数；并且体系大小Nx, Ny, Nz必须和FESim的输入文件一致。

> ini格式读取使用GitHub上一个开源库[inih][2]，只需单个头文件，已经包含在项目文件里

## 原理

这里大概介绍一下仿真BaTiO3铁电相变的原理。首先将哈密顿量近似展开成局域声子模式和应变张量的多项式形式，然后用基于密度泛函理论的第一性原理计算得到其中的系数。声子模式与应变可以认为是晶体内部原子位移的表示，至于为什么选择声子模式与应变而不是简单的原子位移，因为声子模式对称性较高，可以得到简洁的展开式；而且根据软模理论，软声子模式在铁电相变中起主要作用，因此可以只考虑软声子而忽略高阶声子。

得到哈密顿量后，就可以采用蒙特卡洛模拟或者分子动力学对相变进行仿真。原论文[Zhang et al., 1995][1]采用的是蒙特卡洛模拟。

## 结果展示

![](https://github.com/yuchenxi2000/fesim/blob/main/fesim/fig.png)

10x10x10超胞体系使用`fesim/run-BTO.sh`得到的仿真结果，降温过程中应变随温度的变化。应变反映了晶格常数的变化，展示了降温过程中从C，TO到R相的相变过程。该图像用`fesim/plot.py`绘制。

![](https://github.com/yuchenxi2000/fesim/blob/main/fesim/fig_prb.png)

作为对照，上图是[Zhang et al., 1995][1]中的结果。我发现上图中，6个应变分量的顺序貌似和原文定义不一致。原文A部分中定义6个应变分量为exx，eyy，ezz，eyz，ezx，exy。但是图中貌似为exx，eyy，ezz，exy，eyz，ezx。

[1]: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.52.6301 "Effective Hamiltonian theory"
[2]: https://github.com/benhoyt/inih "Ini library"

