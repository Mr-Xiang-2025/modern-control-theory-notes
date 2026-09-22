# 第二章：线性系统运动规律分析

> 来源：《现代控制理论-线性系统-ch2.pptx》，共 40 页。页码均为 PPT 物理页序。本文保留课件 2.1–2.3 的结构、主要术语、定义 2.1–2.2 和结论 2.1–2.16；在此基础上补充公式来源、推导逻辑、物理意义及知识联系。
>
> 文末完整收录例 2.1–2.6，并另列第 10、12 页的两组未编号重根演算示例。课件原题与校正内容明确区分。公式使用标准 LaTeX，按 GitHub 支持的数学语法排版，请使用 **Preview／预览** 阅读。
>
> [第一章笔记](第一章_线性系统状态空间描述.md) · [第二章逐页核对索引](Chapter02_逐页核对索引.md)

## 阅读目录

- [本章主线与符号](#本章主线与符号)
- [2.1 连续时间线性定常系统运动规律分析](#21-连续时间线性定常系统运动规律分析)
- [2.2 线性系统状态转移矩阵](#22-线性系统状态转移矩阵)
- [2.3 连续时间线性时变系统运动规律分析](#23-连续时间线性时变系统运动规律分析)
- [本章小结与学习检查](#本章小结与学习检查)
- [课件校勘与适用条件](#课件校勘与适用条件)
- [附录：全部例题及详细计算](#附录全部例题及详细计算)

## 本章主线与符号

第一章建立模型，回答“系统内部如何变化”；第二章求解模型，回答“从给定初始状态出发，施加给定输入后，各时刻的状态究竟是多少”。

核心区别是：状态方程给出**瞬时变化率**，状态响应给出**经过一段时间后的实际状态**。已知

```math
\dot x(t)=Ax(t)+Bu(t),\qquad x(t_0)=x_0,
```

不能把 $`Ax`$ 当作下一时刻的状态。它是当前状态变化的速度，需要对整个运动过程求解，才能得到 $`x(t)`$。

本章贯穿的公式是

```math
\boxed{
x(t)=\underbrace{\Phi(t,t_0)x_0}_{\text{初始状态的影响}}
+\underbrace{\int_{t_0}^{t}\Phi(t,\tau)B(\tau)u(\tau)\,d\tau}_{\text{各时刻输入影响的累积}}
}.
```

对定常系统，$`\Phi(t,t_0)=e^{A(t-t_0)}`$，$`B(\tau)=B`$。先理解这个特例，再推广到时变系统，就是课件的顺序。

| 课件部分 | 页码 | 解决的问题 |
| --- | --- | --- |
| 2.1 连续时间线性定常系统运动规律分析 | 1–20 | 如何求零输入、零状态及全响应？如何计算矩阵指数？ |
| 2.2 线性系统状态转移矩阵 | 21–31 | 如何用统一的“转移”概念表达运动？基本解阵为什么不唯一？ |
| 2.3 连续时间线性时变系统运动规律分析 | 32–37 | 系统矩阵随时间变化后，哪些公式仍成立，哪些要改变？ |
| 本章小结 | 38–40 | 三部分的公式与计算方法如何对应？ |

符号沿用课件：$`x\in\mathbb R^n`$ 为状态，$`u\in\mathbb R^p`$ 为输入，$`A`$ 为系统矩阵，$`B`$ 为输入矩阵；$`x_{0u}`$ 为零输入响应，$`x_{0x}`$ 为零状态响应；$`P,Q`$ 为坐标变换矩阵，$`\Psi`$ 为基本解阵，$`\Phi`$ 为状态转移矩阵。

本文使用 $`I`$ 表示单位矩阵，$`0`$ 根据上下文表示零标量、零向量或零矩阵。“定常”“时不变”在本章表达同一含义。讨论时变系统时，假定系数在所研究区间连续，或满足通常的分段连续条件，输入可积，从而初值问题有唯一解。

## 2.1 连续时间线性定常系统运动规律分析

### 2.1.1 全响应为什么能拆成两部分

对应课件第 2 页。

零输入响应保留原初始状态，关闭输入：

```math
\dot x_{0u}=Ax_{0u},\qquad x_{0u}(t_0)=x_0.
```

零状态响应保留输入，把初始状态设为零：

```math
\dot x_{0x}=Ax_{0x}+Bu,\qquad x_{0x}(t_0)=0.
```

将两式相加，

```math
\frac{d}{dt}(x_{0u}+x_{0x})
=A(x_{0u}+x_{0x})+Bu,
\qquad (x_{0u}+x_{0x})(t_0)=x_0.
```

它与原系统满足同一微分方程和同一初值，由解的唯一性，

```math
\boxed{x(t)=x_{0u}(t)+x_{0x}(t)}.
```

这就是叠加原理在初值问题中的具体应用。

物理上，电容的初始电压、电感的初始电流可以在没有外部输入时继续引起运动，对应零输入响应；外部电源在初始储能为零的情况下推动系统，对应零状态响应。对一般坐标变换后的状态，“零状态”应理解为所选参考状态为零，不必把每个坐标都直接解释为储能量。

**零状态响应并不等于稳态响应。** 零状态响应通常同时包含建立过程和长期行为；“零输入／零状态”按响应的来源划分，“暂态／稳态”按时间行为划分。

### 2.1.2 矩阵指数函数的定义与本质

对应课件第 3–4 页，结论 2.1【零输入响应】。

标量方程 $`\dot x=ax`$ 的解为 $`e^{at}x_0`$。对矩阵，定义

```math
\boxed{
e^{At}=I+At+\frac{A^2t^2}{2!}+\frac{A^3t^3}{3!}+\cdots
=\sum_{k=0}^{\infty}\frac{A^kt^k}{k!}
}.
```

这里的 $`A^2`$ 是矩阵乘积 $`AA`$，矩阵指数也不是把 $`A`$ 的每个元素分别取指数。特别地，$`e^{A\cdot0}=I`$，不是全 1 矩阵。

**按课件的待定系数法推导。** 先取 $`t_0=0`$，设

```math
x_{0u}(t)=b_0+b_1t+b_2t^2+\cdots,
```

其中 $`b_k`$ 是列向量。代入 $`\dot x=Ax`$，比较 $`t^k`$ 的系数：

```math
(k+1)b_{k+1}=Ab_k,\qquad
b_k=\frac{A^k}{k!}b_0.
```

初值给出 $`b_0=x_0`$，所以

```math
x_{0u}(t)=
\left(I+At+\frac{A^2t^2}{2!}+\cdots\right)x_0
=e^{At}x_0.
```

定常系统允许平移时间原点，因此一般形式为

```math
\boxed{x_{0u}(t)=e^{A(t-t_0)}x_0,\qquad t\ge t_0}.
```

**为什么是矩阵？** 不同初始状态分量会相互影响，不能用一个标量放大倍数描述全部方向。$`e^{A(t-t_0)}`$ 的第 $`j`$ 列，就是把初始状态设为第 $`j`$ 个单位向量时，经过这段时间得到的状态。

例如，若

```math
e^{At}=
\begin{bmatrix}
\phi_{11}(t)&\phi_{12}(t)\\
\phi_{21}(t)&\phi_{22}(t)
\end{bmatrix},
```

则 $`x_1(t)=\phi_{11}(t)x_1(0)+\phi_{12}(t)x_2(0)`$。其中 $`\phi_{12}`$ 描述第二个初始状态如何影响第一个当前状态。

从局部看，

```math
x(t+\Delta t)\approx (I+A\Delta t)x(t).
```

把许多很短的运动连续组合起来，极限就是矩阵指数。于是 $`A`$ 描述局部变化，$`e^{At}`$ 描述有限时间内累积的变化。

### 2.1.3 零输入响应的几何意义、运动形态与稳定性

对应课件第 5 页。

在状态空间中，$`x_{0u}(t)`$ 是从 $`x_0`$ 出发的一条轨迹。它是由初始偏离引起的自由运动，运动规律由 $`A`$ 决定，具体轨迹还取决于 $`x_0`$ 激发了哪些方向。

对原点平衡状态，课件给出的渐近稳定条件为

```math
\boxed{\lim_{t\to\infty}e^{At}=0}.
```

必须理解为：**任意初始状态的自由响应都趋于零**。某一条轨迹趋于零，不足以判断整个系统渐近稳定。

结合第一章特征结构，可以进一步理解为

```math
e^{At}\to0
\quad\Longleftrightarrow\quad
\text{所有特征值的实部均小于零}.
```

实特征值对应指数衰减或增长；共轭复特征值对应带有指数包络的振荡；约当块还会带来 $`t,t^2,\ldots`$ 因子。负实部指数最终压过任何有限次幂，因此有限阶约当块不会改变“全部特征值实部严格为负”的渐近稳定判据。

补充边界：$`\dot x=Ax`$ 的平衡状态满足 $`Ax_e=0`$。原点总是平衡状态；当 $`A`$ 奇异时，可能还有其他平衡状态。课件本处讨论的是原点渐近稳定。

### 2.1.4 矩阵指数函数的性质

对应课件第 6 页。

```math
\begin{aligned}
e^{A0}&=I,\\
(e^{At})^{-1}&=e^{-At},\\
\frac{d}{dt}e^{At}&=Ae^{At}=e^{At}A,\\
\frac{d}{dt}e^{-At}&=-Ae^{-At}=-e^{-At}A,\\
e^{A(t+r)}&=e^{At}e^{Ar}=e^{Ar}e^{At},\\
(e^{At})^m&=e^{Amt},\qquad m=0,1,2,\ldots.
\end{aligned}
```

这些性质来自幂级数以及 $`A`$ 与其各次幂可交换。注意：即使 $`A`$ 不可逆，$`e^{At}`$ 在每个有限时刻仍然可逆。

对两个同维矩阵，课件特别要求可交换：

```math
AF=FA
\quad\Longrightarrow\quad
e^{(A+F)t}=e^{At}e^{Ft}=e^{Ft}e^{At}.
```

不能省略这一条件。一般矩阵乘法不交换，标量指数公式不能不加条件地照搬。

一个很有用的反求关系是

```math
\boxed{A=\left.\frac{d}{dt}e^{At}\right|_{t=0}}.
```

但在使用前，应先检查题给函数确实可能是矩阵指数，至少要满足 $`E(0)=I`$。例 2.2 就涉及这个问题。

### 2.1.5 矩阵指数的四类计算方法

对应课件第 6–15 页，结论 2.2–2.6。

#### 方法一：定义法——结论 2.2

直接计算 $`A,A^2,A^3,\ldots`$ 并代入级数。它说明矩阵指数是什么，但对一般矩阵逐项求和计算量大。

当 $`A`$ 为幂零矩阵，即 $`A^r=0`$ 时，级数恰好截断：

```math
e^{At}=I+At+\cdots+\frac{A^{r-1}t^{r-1}}{(r-1)!}.
```

这里的有限项是精确结果。对一般矩阵只取若干项，则是近似。

#### 方法二：特征值法（1）——结论 2.3

课件先考虑特征值两两互异。取对应右特征向量组成

```math
P=[v_1,\ldots,v_n],\qquad
P^{-1}AP=\Lambda=\mathrm{diag}(\lambda_1,\ldots,\lambda_n).
```

由 $`A=P\Lambda P^{-1}`$，

```math
A^2=P\Lambda^2P^{-1},\qquad
A^k=P\Lambda^kP^{-1}.
```

代入矩阵指数定义，可把左右两端的 $`P,P^{-1}`$ 提取出来：

```math
\boxed{
e^{At}=Pe^{\Lambda t}P^{-1}
=P
\begin{bmatrix}
e^{\lambda_1t}&&0\\
&\ddots&\\
0&&e^{\lambda_nt}
\end{bmatrix}
P^{-1}
}.
```

物理含义按乘法顺序从右向左理解：$`P^{-1}`$ 把原状态换到模态坐标；$`e^{\Lambda t}`$ 让各模态独立演化；$`P`$ 再换回原坐标。

补充：真正需要的条件是 $`A`$ 可对角化。特征值互异是充分条件；有重根的矩阵也可能可对角化，不能看到重根就断定必须出现非平凡约当块。

#### 方法二：特征值法（2）——结论 2.4

普通特征向量不足以组成基时，使用第一章的广义特征向量，构成可逆矩阵 $`Q`$：

```math
Q^{-1}AQ=J,\qquad e^{At}=Qe^{Jt}Q^{-1}.
```

一个 $`r`$ 阶约当块写成

```math
J_r(\lambda)=\lambda I+N,\qquad N^r=0.
```

因 $`\lambda I`$ 与 $`N`$ 可交换，

```math
\boxed{
e^{J_r(\lambda)t}
=e^{\lambda t}\left(I+tN+\frac{t^2}{2!}N^2+\cdots+
\frac{t^{r-1}}{(r-1)!}N^{r-1}\right)
}.
```

例如三阶块：

```math
e^{J_3(\lambda)t}
=e^{\lambda t}
\begin{bmatrix}
1&t&t^2/2\\
0&1&t\\
0&0&1
\end{bmatrix}.
```

这解释了广义特征向量为何在运动分析中有用：约当链使一个模态沿链影响前面的模态，连续积分产生 $`t`$、$`t^2/2`$ 等因子。

与第一章的符号联系：若链按 $`(\lambda I-A)v^{(k)}=v^{(k-1)}`$ 定义，直接以这些向量为列可能得到超对角元为 $`-1`$ 的块。要使用本章超对角元为 $`+1`$ 的标准块，需要按相应的交替符号调整链向量，或直接验证 $`AQ=QJ`$。不可混用两种链约定。

第 10 页的 $`J_3(\lambda_1)\oplus J_2(\lambda_2)`$ 演算见附录 S01。

#### 方法三：有限项展开法／多项式法——结论 2.5

凯莱–哈密尔顿定理给出

```math
p_A(A)=A^n+c_{n-1}A^{n-1}+\cdots+c_0I=0.
```

所以 $`A^n`$ 及更高次幂都能用 $`I,A,\ldots,A^{n-1}`$ 表示。由此可将矩阵指数精确写为

```math
\boxed{
e^{At}=\alpha_0(t)I+\alpha_1(t)A+\cdots+
\alpha_{n-1}(t)A^{n-1}
}.
```

**有限的是矩阵幂的个数，系数仍然是完整的时间函数。** 这不是把泰勒级数截断到 $`n-1`$ 项。

令

```math
q_t(z)=\sum_{k=0}^{n-1}\alpha_k(t)z^k.
```

若特征值两两互异，则要求 $`q_t(\lambda_i)=e^{\lambda_it}`$，得到范德蒙方程：

```math
\begin{bmatrix}
1&\lambda_1&\cdots&\lambda_1^{n-1}\\
1&\lambda_2&\cdots&\lambda_2^{n-1}\\
\vdots&\vdots&&\vdots\\
1&\lambda_n&\cdots&\lambda_n^{n-1}
\end{bmatrix}
\begin{bmatrix}\alpha_0\\\alpha_1\\\vdots\\\alpha_{n-1}\end{bmatrix}
=
\begin{bmatrix}e^{\lambda_1t}\\e^{\lambda_2t}\\\vdots\\e^{\lambda_nt}\end{bmatrix}.
```

对重根只重复相同的取值方程不能提供新信息。约当块还需要匹配导数：

```math
\boxed{
\left.\frac{d^\ell q_t(z)}{dz^\ell}\right|_{z=\lambda_i}
=t^\ell e^{\lambda_it},
\qquad \ell=0,\ldots,r_i-1
}.
```

这里 $`r_i`$ 可取该特征值最大约当块的阶数，即最小多项式中相应因子的次数。若按特征多项式的代数重数施加更多导数条件，也可得到一个有效的插值多项式，但不一定是所需的最低次数。

推导依据是

```math
q_t(\lambda I+N)=
\sum_{\ell=0}^{r_i-1}\frac{q_t^{(\ell)}(\lambda)}{\ell!}N^\ell,
```

与 $`e^{(\lambda I+N)t}`$ 逐项比较。导数是对插值变量 $`z`$ 求的，不能误当作对时间 $`t`$ 求导。

课件第 12 页的三重根、二重根情形及其系数矩阵校正，见附录 S02。

#### 方法四：预解矩阵法——结论 2.6

预解矩阵为 $`(sI-A)^{-1}`$。矩阵指数满足

```math
\boxed{e^{At}=\mathcal L^{-1}\{(sI-A)^{-1}\}},\qquad t\ge0.
```

可以直接由微分方程推导。设 $`E(t)=e^{At}`$，则 $`\dot E=AE`$、$`E(0)=I`$，拉普拉斯变换给出

```math
sE(s)-I=AE(s)
\quad\Longrightarrow\quad
E(s)=(sI-A)^{-1}.
```

课件使用级数证明：在展开收敛的区域，

```math
(sI-A)^{-1}
=\frac{I}{s}+\frac{A}{s^2}+\frac{A^2}{s^3}+\cdots.
```

利用 $`\mathcal L^{-1}\{1/s^{k+1}\}=t^k/k!`$，逐项恢复矩阵指数。实际计算时，通常先求逆，再对每个元素做部分分式展开和逆变换。拉普拉斯变换的收敛域应取在全部特征值实部右侧。

四类方法的关系：

| 方法 | 利用的结构 | 计算重点 |
| --- | --- | --- |
| 定义法 | 幂级数 | 观察高次幂是否截断或有简单规律 |
| 特征值法 | 相似变换与模态／约当块 | 求正确且可逆的特征基或广义特征基 |
| 有限项展开法 | 凯莱–哈密尔顿定理与插值 | 求时间系数，重根时补导数条件 |
| 预解矩阵法 | 时域与拉普拉斯域对应 | 求逆、部分分式、逐元素逆变换 |

它们计算同一个 $`e^{At}`$。无论用哪种方法，最后都能检查

```math
E(0)=I,\qquad \dot E(t)=AE(t).
```

仅满足这两个条件，就由初值问题唯一性确定了正确的矩阵指数。

### 2.1.6 零状态响应及卷积公式的来源

对应课件第 17–19 页，结论 2.7【零状态响应】。

考虑 $`\dot x=Ax+Bu`$、$`x(0)=0`$。为了去掉 $`Ax`$，在左侧引入积分因子 $`e^{-At}`$：

```math
\begin{aligned}
\frac{d}{dt}(e^{-At}x)
&=-Ae^{-At}x+e^{-At}\dot x\\
&=e^{-At}(\dot x-Ax)\\
&=e^{-At}Bu.
\end{aligned}
```

这里用到 $`A`$ 与 $`e^{-At}`$ 可交换，并没有交换 $`B`$ 与指数矩阵。积分后

```math
e^{-At}x(t)=\int_0^t e^{-A\tau}Bu(\tau)\,d\tau.
```

左乘 $`e^{At}`$：

```math
\boxed{x_{0x}(t)=\int_0^t e^{A(t-\tau)}Bu(\tau)\,d\tau}.
```

一般初始时刻：

```math
\boxed{x_{0x}(t)=\int_{t_0}^{t}e^{A(t-\tau)}Bu(\tau)\,d\tau}.
```

**积分里的每一项是什么？**

在时刻 $`\tau`$，输入在短时间 $`d\tau`$ 内向状态注入约 $`Bu(\tau)d\tau`$ 的增量；这个增量继续自由演化到观察时刻 $`t`$，变为

```math
e^{A(t-\tau)}Bu(\tau)d\tau.
```

把所有历史时刻的贡献相加，就是积分。注意时间差是 $`t-\tau`$，而不是固定的 $`t-t_0`$：较早和较晚注入的影响，经历的传播时间不同。

对 $`t_0=0`$，更换积分变量还能写成

```math
x_{0x}(t)
=\int_0^t e^{A\tau}B\,u(t-\tau)\,d\tau.
```

这就是卷积形式。更换变量不等于允许任意交换矩阵乘法顺序。

#### 几何属性、强迫运动与能达性

零状态响应在状态空间中从原点出发，由输入推动形成轨迹。它是强迫运动，依赖输入，也依赖系统本身。

课件写“稳态取决于输入”，应结合系统理解：若稳态存在，它通常同时依赖 $`A,B`$ 和输入。例如在渐近稳定且输入恒为 $`\bar u`$ 时，

```math
0=Ax_{\mathrm{ss}}+B\bar u
\quad\Longrightarrow\quad
x_{\mathrm{ss}}=-A^{-1}B\bar u.
```

若系统不稳定或输入持续变化，则不一定存在常值稳态。

课件第 19 页引入能达性：若对任意目标状态 $`x_f`$，都能选择某个输入及有限终止时刻 $`t_f`$，使从原点出发的响应达到 $`x_f`$，则系统具有完全能达性。这个性质取决于 $`(A,B)`$，不是任何线性系统自动具有的性质。

与后续章节的联系在于：$`B`$ 决定输入首先进入哪些方向，$`A`$ 决定这些方向如何传播、耦合到其他方向。本章的积分公式正是讨论“能否到达某状态”的出发点。

### 2.1.7 系统全响应

对应课件第 20 页，结论 2.8【系统状态运动规律】。

```math
\boxed{
x(t)=e^{At}x_0+\int_0^t e^{A(t-\tau)}Bu(\tau)\,d\tau
},\qquad t_0=0,
```

```math
\boxed{
x(t)=e^{A(t-t_0)}x_0+
\int_{t_0}^{t}e^{A(t-\tau)}Bu(\tau)\,d\tau
},\qquad t\ge t_0.
```

第一项传播初始状态，第二项累计输入影响；两者共同满足初始条件与状态方程。

如果还要求输出，使用第一章的输出方程即可：

```math
y(t)=Cx(t)+Du(t).
```

因此状态响应与输出响应应区分；本章例题没有给出 $`C,D`$ 时，不应擅自把某个状态认定为输出。

## 2.2 线性系统状态转移矩阵

### 2.2.1 为什么引入状态转移矩阵

对应课件第 21–22 页。

$`e^{A(t-t_0)}`$ 的作用是把零输入条件下时刻 $`t_0`$ 的状态变成时刻 $`t`$ 的状态。把这个作用命名为“状态转移”，有两个目的：

1. 直接表达“从哪个时刻传播到哪个时刻”的物理含义。
2. 当 $`A`$ 随时间变化时，仍使用同一种响应结构。

本节沿用课件定常系统的记法 $`\Phi(t-t_0)`$；下一节使用两个时间自变量 $`\Phi(t,t_0)`$。定常系统只是两时间函数仅依赖时间差的特殊情况。

### 2.2.2 定义 2.1：状态转移矩阵

对应课件第 23 页。

给定系统矩阵 $`A`$，满足

```math
\frac{d}{dt}\Phi(t-t_0)=A\Phi(t-t_0),\qquad \Phi(0)=I
```

的 $`n\times n`$ 解阵，称为系统的状态转移矩阵。

其中 $`\Phi(0)=I`$ 的含义是“传播零时长时保持原状态”，不是“绝对时间 $`t=0`$ 时必定这样”。在此处，$`\Phi`$ 的自变量是时间差。

对定常系统：

```math
\boxed{\Phi(t-t_0)=e^{A(t-t_0)}}.
```

状态转移矩阵取决于齐次状态方程，即 $`A`$；$`B`$ 决定输入怎样注入状态，随后这些影响也由同一个 $`\Phi`$ 传播。

### 2.2.3 定义 2.2：基本解阵

对应课件第 23–24 页，结论 2.9–2.11。

取任意非奇异实常阵 $`H`$，考虑

```math
\dot\Psi(t)=A\Psi(t),\qquad \Psi(t_0)=H.
```

其解 $`\Psi(t)`$ 是一个基本解阵。等价地，取齐次系统的 $`n`$ 个线性无关解为列：

```math
\Psi(t)=[\psi_1(t),\ldots,\psi_n(t)].
```

每列都满足 $`\dot\psi_i=A\psi_i`$，所以整个矩阵满足 $`\dot\Psi=A\Psi`$。初始列向量线性无关，使 $`H`$ 可逆；在同一区间内它们保持线性无关。

齐次系统的任意状态解都可以写成

```math
x(t)=\Psi(t)c,
```

其中 $`c`$ 是与时间无关的系数向量。这就是“基本”的含义：它的列提供全部齐次解的一组基。

**结论 2.9：基本解阵不唯一。** 可以选择不同的初始基 $`H`$，得到不同 $`\Psi`$。但如果已经固定同一个 $`H`$，对应初值问题的解是唯一的；“不唯一”来自初始基可选，不是微分方程失去唯一性。

**结论 2.10：基本解阵由线性无关解构成。** 有 $`n`$ 条轨迹还不够，它们的初始状态必须线性无关，否则无法表示任意初始状态。

**结论 2.11：基本解阵的形式。** 指定 $`\Psi(t_0)=H`$ 时，

```math
\boxed{\Psi(t)=e^{A(t-t_0)}H}.
```

课件列出的 $`\Psi(t)=e^{At}`$ 确实是一个基本解阵；它对应的初值为 $`H=e^{At_0}`$，不能把它当作任意指定 $`H`$ 的同一个解。

### 2.2.4 基本解阵与状态转移矩阵的关系

对应课件第 25–26 页，结论 2.12–2.13。

由 $`x(t)=\Psi(t)c`$，在初始时刻

```math
x_0=\Psi(t_0)c
\quad\Longrightarrow\quad
c=\Psi^{-1}(t_0)x_0.
```

代回：

```math
x(t)=\Psi(t)\Psi^{-1}(t_0)x_0.
```

因此

```math
\boxed{\Phi(t-t_0)=\Psi(t)\Psi^{-1}(t_0)}.
```

可以把它理解为两步：先用 $`\Psi^{-1}(t_0)`$ 将初始状态分解到所选解基，再用 $`\Psi(t)`$ 把同一组系数还原为当前状态。

**为什么基本解阵不唯一，而状态转移矩阵唯一？**

设另一个基本解阵为

```math
\Psi_2(t)=\Psi_1(t)P,
```

其中 $`P`$ 是可逆常阵，表示更换齐次解的基。于是

```math
\begin{aligned}
\Psi_2(t)\Psi_2^{-1}(t_0)
&=\Psi_1(t)P\,P^{-1}\Psi_1^{-1}(t_0)\\
&=\Psi_1(t)\Psi_1^{-1}(t_0).
\end{aligned}
```

更换基带来的因素完全抵消。也可直接检查 $`\Phi`$ 的微分方程和单位初值，由初值问题唯一性得到同一结论。

### 2.2.5 定常系统的表达式与性质

对应课件第 27–28 页，结论 2.14–2.15。

**结论 2.14【状态转移矩阵的形式】**

```math
\Phi(t)=e^{At}\quad(t_0=0),\qquad
\Phi(t-t_0)=e^{A(t-t_0)}\quad(t_0\ne0).
```

**结论 2.15【系统状态运动规律】**

```math
\boxed{
x(t)=\Phi(t-t_0)x_0+
\int_{t_0}^{t}\Phi(t-\tau)Bu(\tau)\,d\tau
}.
```

取 $`t_0=0`$ 即得到课件的另一种表达形式。

状态转移矩阵的性质：

```math
\begin{aligned}
\Phi(0)&=I,\\
\Phi^{-1}(t-t_0)&=\Phi(t_0-t),\\
\Phi(t_2-t_1)\Phi(t_1-t_0)&=\Phi(t_2-t_0),\\
\Phi(t_2+t_1)&=\Phi(t_2)\Phi(t_1),\\
\Phi(mt)&=[\Phi(t)]^m,\\
\frac{d}{dt}\Phi(t-t_0)&=A\Phi(t-t_0)=\Phi(t-t_0)A.
\end{aligned}
```

第三条是最重要的**转移组合性质**：先从 $`t_0`$ 走到 $`t_1`$，再从 $`t_1`$ 走到 $`t_2`$，等效于直接从 $`t_0`$ 走到 $`t_2`$。矩阵按实际作用顺序从右向左排列。

逆矩阵性质说明，在精确的齐次数学模型中，有限时间状态传播可逆。它并不意味着有噪声时可以无误差恢复过去，也不保证逆向计算数值稳定。

### 2.2.6 从已知轨迹反求状态转移矩阵和系统矩阵

对应课件第 29–30 页，例 2.4 及其推论。

把 $`n`$ 个线性无关初始状态和各自响应分别按列组合：

```math
X_0=[x^{(1)}(0),\ldots,x^{(n)}(0)],
\qquad
X(t)=[x^{(1)}(t),\ldots,x^{(n)}(t)].
```

因为 $`X(t)=e^{At}X_0`$，

```math
\boxed{\Phi(t)=X(t)X_0^{-1}},
\qquad
\boxed{A=\dot X(0)X_0^{-1}}.
```

也可以在其他时刻使用 $`A=\dot X(t)X^{-1}(t)`$；对于真正的定常系统，这个结果不应随时间变化。

若给出的初始状态线性相关，则 $`X_0`$ 不可逆，不能用这组轨迹唯一确定整个状态空间上的传播。这与第一章“足够多的独立方向才能组成坐标基”是同一件事。

### 2.2.7 基于特征结构的状态响应表达式

对应课件第 31 页。

在 $`A`$ 可对角化的条件下，

```math
P=[v_1,\ldots,v_n],\qquad
P^{-1}=
\begin{bmatrix}
\bar v_1^{\mathsf T}\\
\vdots\\
\bar v_n^{\mathsf T}
\end{bmatrix}.
```

这里 $`\bar v_i^{\mathsf T}`$ 沿用课件符号，表示与右特征向量配对的逆矩阵行向量，并非一般意义上的复共轭。它们满足

```math
\bar v_i^{\mathsf T}v_j=\delta_{ij},\qquad
\bar v_i^{\mathsf T}A=\lambda_i\bar v_i^{\mathsf T}.
```

所以

```math
\boxed{e^{At}=\sum_{i=1}^{n}v_i\bar v_i^{\mathsf T}e^{\lambda_it}}.
```

将 $`B=[b_1,\ldots,b_p]`$ 代入全响应，得到课件的一般初始时刻表达式：

```math
x(t)=\sum_{i=1}^{n}v_i\bar v_i^{\mathsf T}
\left[
x_0e^{\lambda_i(t-t_0)}
+\sum_{j=1}^{p}b_j\int_{t_0}^{t}
e^{\lambda_i(t-\tau)}u_j(\tau)\,d\tau
\right].
```

当 $`t_0=0`$，将式中 $`t_0`$ 全部置零即可。

更容易理解的模态形式是，令 $`z=P^{-1}x`$：

```math
\dot z_i=\lambda_i z_i+\bar v_i^{\mathsf T}Bu,\qquad
z_i(t_0)=\bar v_i^{\mathsf T}x_0.
```

于是

```math
z_i(t)=e^{\lambda_i(t-t_0)}\bar v_i^{\mathsf T}x_0+
\int_{t_0}^{t}e^{\lambda_i(t-\tau)}\bar v_i^{\mathsf T}Bu(\tau)\,d\tau,
\qquad
x(t)=\sum_i v_i z_i(t).
```

右特征向量 $`v_i`$ 给出模态在状态空间中的方向，配对左特征向量提取初始状态及输入对该模态的激发程度，$`e^{\lambda_it}`$ 则给出其时间变化。这把第一章的特征结构与本章的运动规律连接起来。

若 $`A`$ 不可对角化，必须改用约当块形式；不能只保留上式中的纯指数项而遗漏多项式因子。

## 2.3 连续时间线性时变系统运动规律分析

### 2.3.1 时变系统为什么需要两个时间变量

对应课件第 32–33 页。

考虑

```math
\dot x(t)=A(t)x(t)+B(t)u(t),\qquad x(t_0)=x_0.
```

定常系统的传播只取决于经过多久；时变系统还取决于这段时间内具体经历了哪些 $`A(t)`$。同样经过一秒，从不同起点出发可能经历不同的系统参数。

因此状态转移矩阵写成 $`\Phi(t,t_0)`$，分别保留观察时刻和起始时刻：

```math
\boxed{
\frac{\partial}{\partial t}\Phi(t,t_0)=A(t)\Phi(t,t_0),
\qquad \Phi(t_0,t_0)=I
}.
```

不要一般地将它写成 $`\Phi(t-t_0)`$。例 2.6 的转移矩阵包含 $`t^2-t_0^2`$，就清楚说明它并非只依赖时间差。

基本解阵方法仍然成立：

```math
\dot\Psi(t)=A(t)\Psi(t),\qquad
\boxed{\Phi(t,t_0)=\Psi(t)\Psi^{-1}(t_0)}.
```

单位传播、逆向传播和组合性质也成立：

```math
\begin{aligned}
\Phi(t_0,t_0)&=I,\\
\Phi^{-1}(t,t_0)&=\Phi(t_0,t),\\
\Phi(t_2,t_1)\Phi(t_1,t_0)&=\Phi(t_2,t_0).
\end{aligned}
```

但不能一般地交换乘积顺序，也不能一般地写成 $`\partial_t\Phi=\Phi A(t)`$。后一等式对定常系统成立，是因为矩阵指数与同一个 $`A`$ 可交换。

### 2.3.2 结论 2.16：时变系统状态响应

```math
\boxed{
x(t)=\Phi(t,t_0)x_0+
\int_{t_0}^{t}\Phi(t,\tau)B(\tau)u(\tau)\,d\tau
}.
```

与定常系统相比，响应的两项结构没有变化，但必须同时修改传播矩阵的时间变量和输入矩阵的时间参数。

**常数变易法推导。** 齐次解写成 $`\Psi(t)c`$；有输入时让系数变为 $`c(t)`$：

```math
x(t)=\Psi(t)c(t).
```

求导并代入：

```math
\dot\Psi c+\Psi\dot c=A(t)\Psi c+B(t)u.
```

利用 $`\dot\Psi=A(t)\Psi`$ 消去齐次部分：

```math
\dot c(t)=\Psi^{-1}(t)B(t)u(t).
```

积分并代回：

```math
\begin{aligned}
x(t)
&=\Psi(t)\Psi^{-1}(t_0)x_0\\
&\quad+\int_{t_0}^{t}\Psi(t)\Psi^{-1}(\tau)B(\tau)u(\tau)\,d\tau,
\end{aligned}
```

即得到结论 2.16。

**按课件建议做验证。** 对结果求导，注意积分上限也随 $`t`$ 变化：

```math
\begin{aligned}
\dot x(t)
&=A(t)\Phi(t,t_0)x_0
+\Phi(t,t)B(t)u(t)\\
&\quad+\int_{t_0}^{t}A(t)\Phi(t,\tau)B(\tau)u(\tau)\,d\tau\\
&=A(t)x(t)+B(t)u(t).
\end{aligned}
```

同时 $`t=t_0`$ 时积分区间长度为零，得到 $`x(t_0)=x_0`$。在系数或输入跳变点，按分段解及相应单侧导数理解。

### 2.3.3 状态转移矩阵的定义法计算

对应课件第 33–34 页。课件给出的是有序迭代积分级数，通常称为 Peano–Baker 级数。

先将矩阵微分方程积分：

```math
\Phi(t,t_0)=I+\int_{t_0}^{t}A(\tau_1)\Phi(\tau_1,t_0)\,d\tau_1.
```

不断将右边的 $`\Phi`$ 再用同一积分方程替换：

```math
\begin{aligned}
\Phi(t,t_0)
={}&I+\int_{t_0}^{t}A(\tau_1)\,d\tau_1\\
&+\int_{t_0}^{t}A(\tau_1)
\int_{t_0}^{\tau_1}A(\tau_2)\,d\tau_2\,d\tau_1\\
&+\int_{t_0}^{t}A(\tau_1)
\int_{t_0}^{\tau_1}A(\tau_2)
\int_{t_0}^{\tau_2}A(\tau_3)\,d\tau_3\,d\tau_2\,d\tau_1+\cdots.
\end{aligned}
```

第三重积分内层的上限必须是 $`\tau_2`$。它表示

```math
t_0\le\tau_3\le\tau_2\le\tau_1\le t,
```

较早作用的矩阵位于乘积右边，较晚作用的矩阵位于左边。这个时间顺序与状态传播的实际先后顺序一致。

也可递推定义：

```math
K_0(t)=I,\qquad
K_{k+1}(t)=\int_{t_0}^{t}A(\tau)K_k(\tau)\,d\tau,
\qquad
\Phi(t,t_0)=\sum_{k=0}^{\infty}K_k(t).
```

若 $`A`$ 为常阵，则

```math
K_k(t)=\frac{A^k(t-t_0)^k}{k!},
```

级数就退化为矩阵指数。因此时变公式包含定常公式，而不是另一套无关的方法。

### 2.3.4 为什么通常不能直接写指数的积分

对标量 $`\dot x=a(t)x`$，可写 $`x(t)=\exp(\int a)\,x_0`$；对矩阵，一般不能直接声称

```math
\Phi(t,t_0)=\exp\left(\int_{t_0}^{t}A(\tau)\,d\tau\right).
```

原因是不同时间的矩阵 $`A(\tau_1),A(\tau_2)`$ 可能不交换，普通指数会丢失有序传播的信息。

一个足够的条件是

```math
A(t_1)A(t_2)=A(t_2)A(t_1)
\quad\text{对区间内所有 }t_1,t_2\text{ 成立}.
```

满足该条件时，才可以使用上面的指数积分形式。例如 $`A(t)=a(t)A_0`$ 就满足条件。例 2.6 属于这种情形，例 2.5 则不属于；两题放在一起正好帮助理解差异。

特别注意：$`e^{A(t)(t-t_0)}`$ 把整段历史都当成当前参数，通常也不是精确状态转移矩阵。

### 2.3.5 时变系统的实际求解步骤

1. 写出齐次标量方程，检查是否具有三角结构，能否逐个积分。
2. 构造 $`n`$ 个独立齐次解，组成 $`\Psi(t)`$，求 $`\Phi(t,t_0)`$；也可直接按单位初值求转移矩阵各列，或用有序积分级数。
3. 求初始状态项 $`\Phi(t,t_0)x_0`$。
4. 在输入项中使用 $`\Phi(t,\tau)B(\tau)u(\tau)`$，按题给区间积分。
5. 相加后同时检查初始条件和原微分方程。

例 2.5 主要训练有序积分与三角方程求解；例 2.6 将基本解阵、任意初始时刻、阶跃输入和全响应组合起来。

## 本章小结与学习检查

对应课件第 38–40 页。

| 对象 | 定常系统 | 时变系统 |
| --- | --- | --- |
| 状态方程 | $`\dot x=Ax+Bu`$ | $`\dot x=A(t)x+B(t)u`$ |
| 状态转移矩阵 | $`\Phi(t-t_0)=e^{A(t-t_0)}`$ | $`\Phi(t,t_0)=\Psi(t)\Psi^{-1}(t_0)`$ |
| 零输入响应 | $`\Phi(t-t_0)x_0`$ | $`\Phi(t,t_0)x_0`$ |
| 零状态响应 | $`\int_{t_0}^{t}\Phi(t-\tau)Bu(\tau)\,d\tau`$ | $`\int_{t_0}^{t}\Phi(t,\tau)B(\tau)u(\tau)\,d\tau`$ |
| 全响应 | 零输入响应与零状态响应之和 | 零输入响应与零状态响应之和 |
| 时间依赖 | 传播只依赖时间差 | 传播一般依赖两个绝对时刻 |

课件结论的定位：

| 编号 | 主题 | 本文位置 |
| --- | --- | --- |
| 2.1 | 零输入响应 | 2.1.2 |
| 2.2 | 定义法 | 2.1.5 |
| 2.3 | 特征值法（1） | 2.1.5 |
| 2.4 | 特征值法（2） | 2.1.5 |
| 2.5 | 有限项展开法 | 2.1.5 |
| 2.6 | 预解矩阵法 | 2.1.5 |
| 2.7 | 零状态响应 | 2.1.6 |
| 2.8 | 系统状态运动规律 | 2.1.7 |
| 2.9 | 基本解阵不唯一性 | 2.2.3 |
| 2.10 | 基本解阵的构成 | 2.2.3 |
| 2.11 | 基本解阵的形式 | 2.2.3 |
| 2.12 | 状态转移矩阵与基本解阵的关系 | 2.2.4 |
| 2.13 | 状态转移矩阵的唯一性 | 2.2.4 |
| 2.14 | 定常系统状态转移矩阵的形式 | 2.2.5 |
| 2.15 | 基于状态转移矩阵的状态运动规律 | 2.2.5 |
| 2.16 | 时变系统状态响应 | 2.3.2 |

复习时应能独立解释：

1. $`Ax`$ 是变化率，为什么 $`e^{At}x_0`$ 才是状态？
2. 为什么全响应能拆开，零状态响应为什么仍有暂态？
3. 卷积中的 $`Bu(\tau)`$ 与 $`\Phi(t,\tau)`$ 分别做什么？
4. 约当块为什么产生 $`te^{\lambda t}`$ 和 $`t^2e^{\lambda t}/2`$？
5. 有限项展开为什么是精确表示，重根为什么需要导数条件？
6. 基本解阵可以换基，为什么实际状态传播结果不变？
7. 时变系统为何通常要保留 $`t,t_0`$ 两个自变量？
8. 怎样用“初值正确、导数满足方程”发现例题答案中的常数或符号错误？

计算时最常用的三组检查是：

```math
E(0)=I,\quad \dot E=AE;
\qquad
\Phi(t_0,t_0)=I,\quad \partial_t\Phi=A(t)\Phi;
\qquad
x(t_0)=x_0,\quad \dot x=A(t)x+B(t)u.
```

这些检查比只对照最后一行答案更可靠。例 2.6 的错误结果满足微分方程却不满足初值，说明两项必须同时检查。

## 课件校勘与适用条件

以下均根据页面图像核对，涉及数学判断的条目附有推导或检验依据。

| 页码 | 原课件问题或容易误解之处 | 本文处理 |
| --- | --- | --- |
| 3、6 | 矩阵指数展开首项有时写成标量 1；第 3 页“奇次方程” | 统一为单位矩阵 $`I`$；术语应为“齐次方程” |
| 5 | “自由运动趋于原点”的量词容易被忽略 | 明确为任意初始状态；单条轨迹不能判定系统渐近稳定 |
| 10、39 | 重根分支可能使读者误以为重根必然不可对角化 | 说明重根矩阵仍可能可对角化；约当型是一般处理 |
| 12 | 二阶导数行末项写成 $`\lambda_1^{n-2}`$；$`\lambda_2`$ 的函数值行第四项重复写成平方 | 分别校正为 $`\lambda_1^{n-3}`$ 和 $`\lambda_2^3`$，详见 S02 |
| 16 | 例 2.2 原题左下元素为 $`e^{-t}+e^{3t}`$，使零时刻矩阵不为 $`I`$；求导行另有指数符号错误 | 原题按字面无解；另按课件预期答案一致的负号修正版完整验算，详见 E02 |
| 19 | 第 5 项标题写“零输入响应更一般形式”，实际公式为零状态响应 | 校正标题；同时补充稳态依赖系统和输入、能达性是条件性性质 |
| 23 | “由于 $`H`$ 任意，所以解不唯一” | 区分基本解阵选择不唯一与固定初值后解唯一 |
| 24 | 指定任意 $`H`$ 后直接写 $`\Psi=e^{At}`$；旁注“齐次方程的通解+特解”不严谨 | 指定初值时应为 $`e^{A(t-t_0)}H`$；齐次通解是 $`\Psi c`$，基本解阵由独立齐次解构成 |
| 29 | 例 2.4 中将转移矩阵误记为 $`\Psi(t-t_0)`$ | 统一为 $`\Phi(t)`$，数值计算本身正确 |
| 31 | 未在模态展开旁重申可对角化条件；一般初始时刻一行仍标 $`t\ge0`$ | 补充条件；研究区间应为 $`t\ge t_0`$ |
| 33、34 | 有序积分级数第三重积分的最内层上限误写为 $`\tau_3`$ | 应为 $`\tau_2`$；依时间嵌套顺序修正 |
| 37 | 例 2.6 最后一行第二分量常数为 $`-5/6`$，漏加零输入响应中的 $`3/2`$ | 正确常数为 $`2/3`$；代入 $`t=1`$ 检查必须得到 2 |

## 附录：全部例题及详细计算

课件明确编号例题共 **6 题**，无跳号：例 2.1–2.6。另有 **2 组未编号重根演算示例**，作为 S01、S02 收录。正文中的一般证明、模态展开及第 30 页反求公式也均已保留，不能误计为遗漏的独立数值例题。

| 本文编号 | 原课件编号／内容 | 页码 |
| --- | --- | --- |
| E01 | 例 2.1：计算矩阵指数，特征值法与预解矩阵法 | 14–15 |
| E02 | 例 2.2：由矩阵指数反求 $`A`$，含原题一致性检查 | 16 |
| E03 | 例 2.3：单位阶跃零状态响应 | 18 |
| E04 | 例 2.4：由两组初值与轨迹反求 $`\Phi,A`$ | 29–30 |
| E05 | 例 2.5：时变三角系统状态转移矩阵 | 34 |
| E06 | 例 2.6：非零初始时刻的时变系统全响应 | 35–37 |
| S01 | 三阶和二阶约当块组合的指数 | 10 |
| S02 | 三重根、二重根的有限项展开系数方程 | 12 |

### E01　例 2.1：两种方法计算矩阵指数

**课件第 14–15 页。原题：**

```math
\dot x=Ax,\qquad A=\begin{bmatrix}0&1\\-2&-3\end{bmatrix},
\qquad \text{求 }e^{At}.
```

#### 解法一：特征值法

先求特征多项式：

```math
\det(\lambda I-A)
=\begin{vmatrix}\lambda&-1\\2&\lambda+3\end{vmatrix}
=\lambda^2+3\lambda+2
=(\lambda+1)(\lambda+2).
```

所以 $`\lambda_1=-1,\lambda_2=-2`$。

对 $`\lambda_1=-1`$：

```math
(-I-A)v_1=0
\quad\Longrightarrow\quad
-v_{11}-v_{12}=0.
```

取 $`v_{11}=1`$，得 $`v_1=[1,-1]^{\mathsf T}`$。对 $`\lambda_2=-2`$：

```math
(-2I-A)v_2=0
\quad\Longrightarrow\quad
-2v_{21}-v_{22}=0,
```

取 $`v_2=[1,-2]^{\mathsf T}`$。

```math
P=\begin{bmatrix}1&1\\-1&-2\end{bmatrix},
\qquad \det P=-1,
\qquad P^{-1}=\begin{bmatrix}2&1\\-1&-1\end{bmatrix}.
```

按相同特征值顺序组成对角矩阵：

```math
\begin{aligned}
e^{At}
&=P\begin{bmatrix}e^{-t}&0\\0&e^{-2t}\end{bmatrix}P^{-1}\\
&=\begin{bmatrix}e^{-t}&e^{-2t}\\-e^{-t}&-2e^{-2t}\end{bmatrix}
\begin{bmatrix}2&1\\-1&-1\end{bmatrix}\\
&=\boxed{
\begin{bmatrix}
2e^{-t}-e^{-2t}&e^{-t}-e^{-2t}\\
-2e^{-t}+2e^{-2t}&-e^{-t}+2e^{-2t}
\end{bmatrix}}.
\end{aligned}
```

#### 解法二：预解矩阵法

```math
sI-A=\begin{bmatrix}s&-1\\2&s+3\end{bmatrix},
```

```math
(sI-A)^{-1}
=\frac{1}{(s+1)(s+2)}
\begin{bmatrix}s+3&1\\-2&s\end{bmatrix}.
```

逐项部分分式展开：

```math
\begin{aligned}
\frac{s+3}{(s+1)(s+2)}&=\frac2{s+1}-\frac1{s+2},\\
\frac1{(s+1)(s+2)}&=\frac1{s+1}-\frac1{s+2},\\
\frac{-2}{(s+1)(s+2)}&=-\frac2{s+1}+\frac2{s+2},\\
\frac{s}{(s+1)(s+2)}&=-\frac1{s+1}+\frac2{s+2}.
\end{aligned}
```

例如第一式设分子为 $`a(s+2)+b(s+1)`$，比较系数得 $`a+b=1`$、$`2a+b=3`$，所以 $`a=2,b=-1`$。其余各项同理。

利用 $`\mathcal L^{-1}\{1/(s+a)\}=e^{-at}`$，得到与解法一完全相同的结果。

#### 补充解法：有限项展开法

这是为联系本章方法补充的解法，课件本题展示的是前两种。

设 $`e^{At}=\alpha_0I+\alpha_1A`$，代入两个特征值：

```math
\alpha_0-\alpha_1=e^{-t},\qquad
\alpha_0-2\alpha_1=e^{-2t}.
```

解得

```math
\alpha_1=e^{-t}-e^{-2t},\qquad
\alpha_0=2e^{-t}-e^{-2t}.
```

因此

```math
e^{At}
=
\begin{bmatrix}
\alpha_0&\alpha_1\\
-2\alpha_1&\alpha_0-3\alpha_1
\end{bmatrix},
```

再次得到相同答案。

**检查与解释。** $`E(0)=I`$；$`E'(0)=A`$；逐元素求导可验证 $`E'=AE`$。全部元素趋于零，两个模态分别以 $`e^{-t}`$ 和 $`e^{-2t}`$ 衰减，系统原点渐近稳定。

### E02　例 2.2：由矩阵指数反求系统矩阵

**课件第 16 页。原题按页面字面为：**

```math
E_{\mathrm{原}}(t)=
\begin{bmatrix}
\frac12(e^{-t}+e^{3t})&\frac14(-e^{-t}+e^{3t})\\
e^{-t}+e^{3t}&\frac12(e^{-t}+e^{3t})
\end{bmatrix},
\qquad E_{\mathrm{原}}(t)=e^{At},\quad \text{求 }A.
```

#### 第一步：先检查原题是否成立

在 $`t=0`$ 时，

```math
E_{\mathrm{原}}(0)=
\begin{bmatrix}1&0\\2&1\end{bmatrix}\ne I.
```

而任何矩阵 $`A`$ 都有 $`e^{A0}=I`$。因此：

```math
\boxed{\text{严格按原题所印矩阵，不存在满足要求的常矩阵 }A.}
```

不能跳过这个矛盾，直接对原题矩阵求导后把结果当作答案。

#### 第二步：按课件预期答案一致的修正版求解

课件给出的答案是

```math
A=\begin{bmatrix}1&1\\4&1\end{bmatrix}.
```

与该答案一致的最小题面修正，是将左下元素改为 $`-e^{-t}+e^{3t}`$：

```math
E(t)=
\begin{bmatrix}
\frac12(e^{-t}+e^{3t})&\frac14(-e^{-t}+e^{3t})\\
-e^{-t}+e^{3t}&\frac12(e^{-t}+e^{3t})
\end{bmatrix}.
```

这是基于课件答案及指数性质的校正，原题的加号已经在上一部分保留。

现在 $`E(0)=I`$，逐元素求导：

```math
E'(t)=
\begin{bmatrix}
-\frac12e^{-t}+\frac32e^{3t}&\frac14e^{-t}+\frac34e^{3t}\\
e^{-t}+3e^{3t}&-\frac12e^{-t}+\frac32e^{3t}
\end{bmatrix}.
```

因此

```math
\boxed{A=E'(0)=\begin{bmatrix}1&1\\4&1\end{bmatrix}}.
```

课件求导行左下位置还把 $`3e^{3t}`$ 印成了 $`3e^{-3t}`$，此处一并改正。虽然在零时刻两者数值相同，时间函数的导数表达式仍不能混写。

#### 第三步：反向验证修正版

写成

```math
A=I+K,\qquad
K=\begin{bmatrix}0&1\\4&0\end{bmatrix},
\qquad K^2=4I.
```

偶次幂与奇次幂分开求和：

```math
e^{At}
=e^t\left[\cosh(2t)I+\frac{\sinh(2t)}2K\right].
```

使用

```math
e^t\cosh(2t)=\frac{e^{3t}+e^{-t}}2,\qquad
e^t\sinh(2t)=\frac{e^{3t}-e^{-t}}2,
```

可恢复修正版 $`E(t)`$ 的全部四个元素。同时 $`\det(\lambda I-A)=(\lambda-3)(\lambda+1)`$，与公式中的 $`e^{3t},e^{-t}`$ 一致。

**本题要点：** 反求 $`A`$ 前要检查零时刻；“候选矩阵函数的导数在零时刻等于某矩阵”本身不能证明它是这个矩阵的指数。

### E03　例 2.3：单位阶跃输入下的零状态响应

**课件第 18 页。原题：**

```math
\begin{bmatrix}\dot x_1\\\dot x_2\end{bmatrix}
=
\begin{bmatrix}0&1\\-2&-3\end{bmatrix}
\begin{bmatrix}x_1\\x_2\end{bmatrix}
+\begin{bmatrix}0\\1\end{bmatrix}u(t),
\qquad x_1(0)=x_2(0)=0,
```

其中 $`u(t)`$ 为单位阶跃输入，求 $`t\ge0`$ 的系统响应。

初值为零，所以全响应等于零状态响应：

```math
x(t)=\int_0^t e^{A(t-\tau)}B\,d\tau.
```

采用 E01 的矩阵指数。因 $`B=[0,1]^{\mathsf T}`$，乘以 $`B`$ 就取指数矩阵的第二列：

```math
e^{A(t-\tau)}B
=
\begin{bmatrix}
e^{-(t-\tau)}-e^{-2(t-\tau)}\\
-e^{-(t-\tau)}+2e^{-2(t-\tau)}
\end{bmatrix}.
```

令 $`r=t-\tau`$，则

```math
\begin{aligned}
x_1(t)
&=\int_0^t(e^{-r}-e^{-2r})\,dr\\
&=(1-e^{-t})-\frac12(1-e^{-2t})\\
&=\frac12-e^{-t}+\frac12e^{-2t},
\end{aligned}
```

```math
\begin{aligned}
x_2(t)
&=\int_0^t(-e^{-r}+2e^{-2r})\,dr\\
&=-(1-e^{-t})+(1-e^{-2t})\\
&=e^{-t}-e^{-2t}.
\end{aligned}
```

最终

```math
\boxed{
x(t)=
\begin{bmatrix}
\frac12-e^{-t}+\frac12e^{-2t}\\
e^{-t}-e^{-2t}
\end{bmatrix},\qquad t\ge0
}.
```

**代回检查：**

```math
x(0)=0,\qquad
\dot x_1=e^{-t}-e^{-2t}=x_2,
```

```math
\dot x_2=-e^{-t}+2e^{-2t}
=-2x_1-3x_2+1.
```

**稳态检查：**

```math
x(\infty)=\begin{bmatrix}1/2\\0\end{bmatrix},
\qquad
A\begin{bmatrix}1/2\\0\end{bmatrix}+B=0.
```

本题的零状态响应包含 $`e^{-t}`$、$`e^{-2t}`$ 的暂态部分；恒定输入还引起非零平衡状态。第一式 $`\dot x_1=x_2`$ 表明两个状态是积分关联关系，是否能具体称为位移和速度，要看原系统的物理定义。

### E04　例 2.4：由两组轨迹反求传播规律

**课件第 29–30 页。原题：** 对定常自治系统 $`\dot x=Ax`$，给出

```math
x^{(1)}(0)=\begin{bmatrix}1\\-4\end{bmatrix},
\qquad
x^{(1)}(t)=\begin{bmatrix}e^{-3t}\\-4e^{-3t}\end{bmatrix},
```

```math
x^{(2)}(0)=\begin{bmatrix}2\\-1\end{bmatrix},
\qquad
x^{(2)}(t)=\begin{bmatrix}2e^{-2t}\\-e^{-2t}\end{bmatrix}.
```

求状态转移矩阵和系统矩阵 $`A`$。

#### 第一步：构造基本解阵并检查独立性

```math
\Psi(t)=
\begin{bmatrix}
e^{-3t}&2e^{-2t}\\
-4e^{-3t}&-e^{-2t}
\end{bmatrix},
\qquad
\Psi(0)=\begin{bmatrix}1&2\\-4&-1\end{bmatrix}.
```

```math
\det\Psi(0)=(-1)-(-8)=7\ne0,
\qquad
\Psi^{-1}(0)=\frac17\begin{bmatrix}-1&-2\\4&1\end{bmatrix}.
```

两组初始状态线性无关，可以构成基本解阵。

#### 第二步：求状态转移矩阵

```math
\begin{aligned}
\Phi(t)
&=\Psi(t)\Psi^{-1}(0)\\
&=\boxed{\frac17
\begin{bmatrix}
-e^{-3t}+8e^{-2t}&-2e^{-3t}+2e^{-2t}\\
4e^{-3t}-4e^{-2t}&8e^{-3t}-e^{-2t}
\end{bmatrix}}.
\end{aligned}
```

课件此行左侧误写为 $`\Psi(t-t_0)`$，实际计算的是 $`\Phi(t)`$。

一般初始时刻的转移矩阵只需将上式所有 $`t`$ 替换为 $`t-t_0`$。

#### 第三步：求系统矩阵

由轨迹直接求导：

```math
\dot\Psi(0)=
\begin{bmatrix}-3&-4\\12&2\end{bmatrix}.
```

于是

```math
\begin{aligned}
A
&=\dot\Psi(0)\Psi^{-1}(0)\\
&=\frac17
\begin{bmatrix}-3&-4\\12&2\end{bmatrix}
\begin{bmatrix}-1&-2\\4&1\end{bmatrix}\\
&=\boxed{\frac17\begin{bmatrix}-13&2\\-4&-22\end{bmatrix}}.
\end{aligned}
```

也可对 $`\Phi(t)`$ 求导再令 $`t=0`$，得到相同结果。

**检查：**

```math
A\begin{bmatrix}1\\-4\end{bmatrix}
=\begin{bmatrix}-3\\12\end{bmatrix}
=-3\begin{bmatrix}1\\-4\end{bmatrix},
```

```math
A\begin{bmatrix}2\\-1\end{bmatrix}
=\begin{bmatrix}-4\\2\end{bmatrix}
=-2\begin{bmatrix}2\\-1\end{bmatrix}.
```

原题两条轨迹恰好沿两个特征方向运动，衰减率分别为 $`-3,-2`$。此外 $`\Phi(0)=I`$、$`\dot\Phi=A\Phi`$，均满足要求。

第 30 页将本题推广到 $`n`$ 组独立轨迹：

```math
\Phi(t)=
[x^{(1)}(t),\ldots,x^{(n)}(t)]
[x^{(1)}(0),\ldots,x^{(n)}(0)]^{-1},
```

```math
A=
[\dot x^{(1)}(0),\ldots,\dot x^{(n)}(0)]
[x^{(1)}(0),\ldots,x^{(n)}(0)]^{-1}.
```

### E05　例 2.5：时变系统的状态转移矩阵

**课件第 34 页。原题：**

```math
\dot x=A(t)x,\qquad
A(t)=\begin{bmatrix}0&1\\0&t\end{bmatrix},
\qquad \text{求 }\Phi(t,0).
```

#### 解法一：按课件的有序积分级数计算

零阶项为 $`I`$，第一阶项：

```math
K_1(t)=\int_0^t
\begin{bmatrix}0&1\\0&\tau\end{bmatrix}d\tau
=\begin{bmatrix}0&t\\0&t^2/2\end{bmatrix}.
```

第二阶项，先计算矩阵乘积：

```math
A(\tau_1)K_1(\tau_1)
=
\begin{bmatrix}0&1\\0&\tau_1\end{bmatrix}
\begin{bmatrix}0&\tau_1\\0&\tau_1^2/2\end{bmatrix}
=
\begin{bmatrix}0&\tau_1^2/2\\0&\tau_1^3/2\end{bmatrix}.
```

所以

```math
K_2(t)=
\begin{bmatrix}0&t^3/6\\0&t^4/8\end{bmatrix}.
```

继续计算第三阶项：

```math
A(\tau_1)K_2(\tau_1)
=\begin{bmatrix}0&\tau_1^4/8\\0&\tau_1^5/8\end{bmatrix},
\qquad
K_3(t)=\begin{bmatrix}0&t^5/40\\0&t^6/48\end{bmatrix}.
```

因此

```math
\Phi(t,0)=
\begin{bmatrix}
1&t+\frac{t^3}{6}+\frac{t^5}{40}+\cdots\\
0&1+\frac{t^2}{2}+\frac{t^4}{8}+\frac{t^6}{48}+\cdots
\end{bmatrix}.
```

课件展示到第二阶项；这里补出第三阶项和精确求和，避免把省略号误认为有限项精确答案。

#### 解法二：直接解三角齐次方程并求出精确形式

标量方程为

```math
\dot x_2=tx_2,\qquad \dot x_1=x_2.
```

先解关于 $`x_2`$ 的独立方程：

```math
x_2(t)=e^{t^2/2}x_2(0).
```

再积分 $`\dot x_1=x_2`$：

```math
x_1(t)=x_1(0)+x_2(0)\int_0^t e^{\sigma^2/2}\,d\sigma.
```

故精确结果为

```math
\boxed{
\Phi(t,0)=
\begin{bmatrix}
1&\displaystyle\int_0^t e^{\sigma^2/2}\,d\sigma\\
0&e^{t^2/2}
\end{bmatrix}
}.
```

右上角积分保留为定积分即是精确表达，无须强行写成初等函数。将

```math
e^{\sigma^2/2}
=\sum_{k=0}^{\infty}\frac{\sigma^{2k}}{2^kk!}
```

逐项积分，得到

```math
\int_0^t e^{\sigma^2/2}\,d\sigma
=\sum_{k=0}^{\infty}\frac{t^{2k+1}}{(2k+1)2^kk!},
```

与前面的有序积分级数完全一致。

**直接检查：** 令 $`F(t)=\int_0^t e^{\sigma^2/2}d\sigma`$，则

```math
\dot\Phi(t,0)=
\begin{bmatrix}0&e^{t^2/2}\\0&te^{t^2/2}\end{bmatrix}
=A(t)\Phi(t,0),\qquad
\Phi(0,0)=I.
```

#### 补充：一般初始时刻与不能套用指数积分的原因

同样求解得到

```math
\Phi(t,t_0)=
\begin{bmatrix}
1&\displaystyle\int_{t_0}^{t}e^{(\sigma^2-t_0^2)/2}\,d\sigma\\
0&e^{(t^2-t_0^2)/2}
\end{bmatrix}.
```

比较两个时刻 $`a,b`$：

```math
A(a)A(b)=\begin{bmatrix}0&b\\0&ab\end{bmatrix},
\qquad
A(b)A(a)=\begin{bmatrix}0&a\\0&ab\end{bmatrix}.
```

当 $`a\ne b`$ 时一般不相等，因而不能直接用 $`\exp(\int_0^t A)`$ 代替正确传播矩阵。

事实上，错误套用指数积分所得右上角为

```math
\frac{2}{t}\bigl(e^{t^2/2}-1\bigr)
=t+\frac{t^3}{4}+\cdots\quad(t\ne0),
```

而正确右上角为 $`t+t^3/6+\cdots`$，从三次项起已经不同。在 $`t=0`$ 处按连续延拓理解该错误候选式。

### E06　例 2.6：时变系统在非零初始时刻的全响应

**课件第 35–37 页。原题：**

```math
\dot x=
\begin{bmatrix}0&0\\t&0\end{bmatrix}x
+\begin{bmatrix}1\\1\end{bmatrix}u(t),
\qquad
x(1)=\begin{bmatrix}1\\2\end{bmatrix},
\qquad t\in[1,10].
```

输入是从 $`t=1`$ 开始作用的单位阶跃 $`u(t)=1(t-1)`$。在本题积分区间内取 $`u(\tau)=1`$；阶跃点单个时刻的取值不影响积分。

#### 第一步：解自治／齐次状态方程

暂时令输入为零：

```math
\dot x_1=0,\qquad \dot x_2=tx_1.
```

所以

```math
x_1(t)=x_1(t_0),
```

```math
x_2(t)=x_2(t_0)+x_1(t_0)\int_{t_0}^{t}\sigma\,d\sigma
=x_2(t_0)+\frac{t^2-t_0^2}{2}x_1(t_0).
```

#### 第二步：按课件所选的两组初值构造基本解阵

取

```math
x^{(1)}(t_0)=\begin{bmatrix}0\\1\end{bmatrix},
\qquad
x^{(2)}(t_0)=\begin{bmatrix}2\\0\end{bmatrix}.
```

对应两条齐次解：

```math
\psi_1(t)=\begin{bmatrix}0\\1\end{bmatrix},
\qquad
\psi_2(t)=\begin{bmatrix}2\\t^2-t_0^2\end{bmatrix}.
```

因此

```math
\Psi(t)=\begin{bmatrix}0&2\\1&t^2-t_0^2\end{bmatrix},
\qquad
\Psi(t_0)=\begin{bmatrix}0&2\\1&0\end{bmatrix},
\qquad
\Psi^{-1}(t_0)=\begin{bmatrix}0&1\\1/2&0\end{bmatrix}.
```

这两组初值仅用于构造解基，并没有替换原题的实际初值 $`[1,2]^{\mathsf T}`$。

#### 第三步：计算状态转移矩阵

```math
\boxed{
\Phi(t,t_0)=\Psi(t)\Psi^{-1}(t_0)
=
\begin{bmatrix}
1&0\\
(t^2-t_0^2)/2&1
\end{bmatrix}
}.
```

输入积分中需要的是

```math
\Phi(t,\tau)=
\begin{bmatrix}1&0\\(t^2-\tau^2)/2&1\end{bmatrix}.
```

把 $`t_0=1`$ 代入初始状态项：

```math
x_{0u}(t)=
\begin{bmatrix}1&0\\(t^2-1)/2&1\end{bmatrix}
\begin{bmatrix}1\\2\end{bmatrix}
=
\begin{bmatrix}1\\t^2/2+3/2\end{bmatrix}.
```

#### 第四步：逐项计算零状态响应

```math
\begin{aligned}
x_{0x}(t)
&=\int_1^t
\begin{bmatrix}1&0\\(t^2-\tau^2)/2&1\end{bmatrix}
\begin{bmatrix}1\\1\end{bmatrix}\,d\tau\\
&=\int_1^t
\begin{bmatrix}1\\(t^2-\tau^2)/2+1\end{bmatrix}\,d\tau.
\end{aligned}
```

第一分量：

```math
x_{0x,1}(t)=t-1.
```

第二分量中 $`t`$ 是积分时固定的观察时刻：

```math
\begin{aligned}
x_{0x,2}(t)
&=\left(\frac{t^2}{2}+1\right)(t-1)
-\frac{t^3-1}{6}\\
&=\frac{t^3}{3}-\frac{t^2}{2}+t-\frac56.
\end{aligned}
```

故

```math
x_{0x}(t)=
\begin{bmatrix}
t-1\\
t^3/3-t^2/2+t-5/6
\end{bmatrix}.
```

#### 第五步：合成全响应并校正课件最后一行

```math
\begin{aligned}
x(t)
&=
\begin{bmatrix}1\\t^2/2+3/2\end{bmatrix}
+
\begin{bmatrix}t-1\\t^3/3-t^2/2+t-5/6\end{bmatrix}\\
&=\boxed{
\begin{bmatrix}
t\\
t^3/3+t+2/3
\end{bmatrix}},\qquad t\in[1,10].
\end{aligned}
```

常数项应相加为

```math
\frac32-\frac56=\frac96-\frac56=\frac23.
```

课件最后一行写成 $`t^3/3+t-5/6`$，漏加了零输入项中的 $`3/2`$。前一行的两项表达正确，错误发生在最后合并。

#### 第六步：双重验证

**初始条件：**

```math
x(1)=
\begin{bmatrix}
1\\1/3+1+2/3
\end{bmatrix}
=\begin{bmatrix}1\\2\end{bmatrix}.
```

课件错误结果在 $`t=1`$ 给出第二分量 $`1/2`$，与题设的 2 不符。

**微分方程：**

```math
\dot x(t)=\begin{bmatrix}1\\t^2+1\end{bmatrix},
\qquad
A(t)x(t)+B
=\begin{bmatrix}1\\t\cdot t+1\end{bmatrix}
=\begin{bmatrix}1\\t^2+1\end{bmatrix}.
```

还可直接求解带输入的标量方程作独立检查：

```math
\dot x_1=1,\quad x_1(1)=1
\quad\Longrightarrow\quad x_1=t,
```

```math
\dot x_2=tx_1+1=t^2+1,\quad x_2(1)=2,
```

```math
x_2(t)=2+\int_1^t(\tau^2+1)\,d\tau
=\frac{t^3}{3}+t+\frac23.
```

**与 E05 的区别。** 本题 $`A(t)=tN`$，其中

```math
N=\begin{bmatrix}0&0\\1&0\end{bmatrix},\qquad N^2=0.
```

不同时间的 $`A`$ 可交换，且任意两者的乘积为零，因此

```math
\Phi(t,t_0)
=\exp\left(N\frac{t^2-t_0^2}{2}\right)
=I+N\frac{t^2-t_0^2}{2}.
```

这个有限项结果是精确的。本题也展示了：$`A(t)`$ 每个时刻的特征值都为零，但传播矩阵仍有随时间增长的项，不能只看瞬时特征值就忽略实际运动。

### S01　第 10 页未编号演算：三阶、二阶约当块组合

**课件给定的结构：** $`n=5`$，$`\lambda_1`$ 代数重数为 3、几何重数为 1；$`\lambda_2`$ 代数重数为 2、几何重数为 1，且两者不同。

由于每个特征值的几何重数就是相应约当块个数，$`\lambda_1`$ 对应一个三阶块，$`\lambda_2`$ 对应一个二阶块：

```math
Q^{-1}AQ=\bar A=J=
\begin{bmatrix}
\lambda_1&1&0&0&0\\
0&\lambda_1&1&0&0\\
0&0&\lambda_1&0&0\\
0&0&0&\lambda_2&1\\
0&0&0&0&\lambda_2
\end{bmatrix}.
```

分别计算幂零部分：

```math
N_3=
\begin{bmatrix}0&1&0\\0&0&1\\0&0&0\end{bmatrix},
\qquad
N_3^2=
\begin{bmatrix}0&0&1\\0&0&0\\0&0&0\end{bmatrix},
\qquad N_3^3=0,
```

```math
N_2=\begin{bmatrix}0&1\\0&0\end{bmatrix},
\qquad N_2^2=0.
```

所以

```math
e^{J_3(\lambda_1)t}
=e^{\lambda_1t}(I+tN_3+t^2N_3^2/2),
\qquad
e^{J_2(\lambda_2)t}=e^{\lambda_2t}(I+tN_2).
```

组合并换回原坐标：

```math
\boxed{
e^{At}=Q
\begin{bmatrix}
e^{\lambda_1t}&te^{\lambda_1t}&\frac{t^2}{2}e^{\lambda_1t}&0&0\\
0&e^{\lambda_1t}&te^{\lambda_1t}&0&0\\
0&0&e^{\lambda_1t}&0&0\\
0&0&0&e^{\lambda_2t}&te^{\lambda_2t}\\
0&0&0&0&e^{\lambda_2t}
\end{bmatrix}
Q^{-1}
}.
```

课件没有提供具体 $`A,Q,\lambda_1,\lambda_2`$ 的数值，因此该结构表达式就是此示例能确定的完整结果，不另行虚构参数。

三阶块在约当坐标中的响应为

```math
\begin{aligned}
z_3(t)&=e^{\lambda_1t}z_3(0),\\
z_2(t)&=e^{\lambda_1t}[z_2(0)+tz_3(0)],\\
z_1(t)&=e^{\lambda_1t}[z_1(0)+tz_2(0)+t^2z_3(0)/2].
\end{aligned}
```

每向链前端传播一级，就多一次积分。这给出上三角多项式项的直观来源。

### S02　第 12 页未编号演算：重根有限项展开系数

**课件给定的结构：** $`\lambda_1`$ 为三重根、几何重数 1；$`\lambda_2`$ 为二重根、几何重数 1；其余 $`\lambda_3,\ldots,\lambda_{n-3}`$ 为不同的简单根。

设

```math
e^{At}=q_t(A),\qquad
q_t(z)=\alpha_0+\alpha_1z+\cdots+\alpha_{n-1}z^{n-1}.
```

#### 第一步：列出全部独立条件

三阶约当块需要函数值、一阶导数、二阶导数：

```math
q_t(\lambda_1)=e^{\lambda_1t},\qquad
q_t'(\lambda_1)=te^{\lambda_1t},\qquad
\frac12q_t''(\lambda_1)=\frac{t^2}{2}e^{\lambda_1t}.
```

二阶块需要函数值及一阶导数：

```math
q_t(\lambda_2)=e^{\lambda_2t},\qquad
q_t'(\lambda_2)=te^{\lambda_2t}.
```

其余简单根只需

```math
q_t(\lambda_i)=e^{\lambda_it},\qquad i=3,\ldots,n-3.
```

方程总数为 $`3+2+(n-5)=n`$，与 $`\alpha_0,\ldots,\alpha_{n-1}`$ 的数量相同。

#### 第二步：按课件的行顺序组成线性方程组

课件按“$`\lambda_1`$ 二阶导数除以 $`2!`$、一阶导数、函数值；$`\lambda_2`$ 一阶导数、函数值；其余函数值”的顺序排列。校正后的矩阵为

```math
V=
\begin{bmatrix}
0&0&1&3\lambda_1&\cdots&
\frac{(n-1)(n-2)}2\lambda_1^{n-3}\\
0&1&2\lambda_1&3\lambda_1^2&\cdots&
(n-1)\lambda_1^{n-2}\\
1&\lambda_1&\lambda_1^2&\lambda_1^3&\cdots&
\lambda_1^{n-1}\\
0&1&2\lambda_2&3\lambda_2^2&\cdots&
(n-1)\lambda_2^{n-2}\\
1&\lambda_2&\lambda_2^2&\lambda_2^3&\cdots&
\lambda_2^{n-1}\\
1&\lambda_3&\lambda_3^2&\lambda_3^3&\cdots&
\lambda_3^{n-1}\\
\vdots&\vdots&\vdots&\vdots&&\vdots\\
1&\lambda_{n-3}&\lambda_{n-3}^2&\lambda_{n-3}^3&\cdots&
\lambda_{n-3}^{n-1}
\end{bmatrix}.
```

于是

```math
\boxed{
\begin{bmatrix}\alpha_0\\\alpha_1\\\vdots\\\alpha_{n-1}\end{bmatrix}
=V^{-1}
\begin{bmatrix}
\frac{t^2}{2}e^{\lambda_1t}\\
te^{\lambda_1t}\\
e^{\lambda_1t}\\
te^{\lambda_2t}\\
e^{\lambda_2t}\\
e^{\lambda_3t}\\
\vdots\\
e^{\lambda_{n-3}t}
\end{bmatrix}
}.
```

当 $`n=5`$ 时不存在其他简单根，删除上述矩阵与向量中从 $`\lambda_3`$ 开始的各行，只保留前五行与五列。

#### 第三步：解释并校正两处幂次

第一行末项来自

```math
\frac1{2!}\frac{d^2}{dz^2}z^{n-1}
=\frac{(n-1)(n-2)}2z^{n-3}.
```

因此课件所写 $`\lambda_1^{n-2}`$ 应为 $`\lambda_1^{n-3}`$。

第五行是直接代入 $`z=\lambda_2`$，不做求导，所以第四个元素应为 $`\lambda_2^3`$，不是课件重复出现的 $`\lambda_2^2`$。

#### 第四步：为什么这组方程确实确定矩阵指数

若两个次数小于 $`n`$ 的多项式满足上述条件，它们之差在 $`\lambda_1`$ 有三重零点、在 $`\lambda_2`$ 有二重零点，在其他根各有一个零点，总重数为 $`n`$。次数小于 $`n`$ 的非零多项式不可能有这么多零点，所以差只能为零，系数唯一。

将所得多项式代入各约当块：

```math
q_t(J_r(\lambda))
=\sum_{\ell=0}^{r-1}
\frac{q_t^{(\ell)}(\lambda)}{\ell!}N^\ell
=e^{\lambda t}\sum_{\ell=0}^{r-1}\frac{t^\ell}{\ell!}N^\ell
=e^{J_r(\lambda)t}.
```

因此 $`q_t(A)=e^{At}`$。课件没有给出这些特征值的数值，本示例最终结果就是上述完整系数方程；它与 S01 的约当块指数表达的是同一个重根传播规律。
