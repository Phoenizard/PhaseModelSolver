# Wang–Du 双组分膜模型 — 数值格式推导

> **文档定位**：本文**推导数值更新格式**——从连续能量耗散律出发，建到一组干净的**时间离散更新公式**。走法：先在最难的子问题上把格式定下来，再逐项拼装到完整模型。方法方向来自 [doc/debate/](../debate/) 的辩论收敛。
> **设计次序**（本文只走到第 ②）：① 确定时间离散 → ② **数值格式（时间离散更新公式；算子 $\Delta,\Delta^2$ 保持连续形式，不做空间离散）** → ③ 空间离散（谱 / FFT）→ ④ 具体求解。**耗散在第 ② 步即可证**，故空间离散与求解一律留到格式定稿之后。
> **记号约定**：模型线张力常数原记 $\delta$，与变分记号 $\delta/\delta\phi$ 冲突，本文记 **$\delta_L$**。$\xi=\varepsilon$，$c_0\equiv0$，$d\in\{2,3\}$（先 2D）。

---

## 1. 连续问题

域 $\Omega=[-\pi,\pi]^d$，周期边界。两个相场 $\phi$（膜界面）、$\eta$（组分序参数）按 $L^2$ 梯度流演化，以极小化总能量 $E_M$：
$$
\phi_t=-\frac{\delta E_M}{\delta\phi},\qquad \eta_t=-\frac{\delta E_M}{\delta\eta}.
$$

**连续能量耗散律**（离散格式要复制的判据）：
$$
\frac{d}{dt}E_M=\Big(\frac{\delta E_M}{\delta\phi},\,\phi_t\Big)+\Big(\frac{\delta E_M}{\delta\eta},\,\eta_t\Big)=-\|\phi_t\|^2-\|\eta_t\|^2\le0.
$$

总能量 = 弯曲能 + 线张力 + 五个罚项：
$$
E_M=E+L+\tfrac12 M_1(V-v_d)^2+\tfrac12 M_2(A-a_0)^2+\tfrac12 M_3(D-a_d)^2+\tfrac12 M_4 N^2+\tfrac12 M_5 P^2.
$$
逐项展开（出现新符号即就地给出）。

**弯曲能 $E$** —— 唯一带变系数四阶主部，是全模型的难点：
$$
E=\int_\Omega\frac{k(\eta)}{2\varepsilon}\,G^2,\qquad
G=\varepsilon\Delta\phi+\tfrac1\varepsilon(\phi-\phi^3),\qquad
k(\eta)=k+c\tanh\tfrac\eta\xi\in[k_2,k_1].
$$
$G$ 里的 $\varepsilon\Delta\phi$ 使 $E$ 含 $\phi$ 的四阶项，$k(\eta)$ 是乘在其上的变系数。

**线张力能 $L$** —— 两个 Ginzburg–Landau 密度的乘积：
$$
L=\delta_L\int_\Omega W_\eta\,W_\phi,\qquad
W_\phi=\tfrac\varepsilon2|\nabla\phi|^2+\tfrac1{4\varepsilon}(\phi^2-1)^2,\quad
W_\eta=\tfrac\xi2|\nabla\eta|^2+\tfrac1{4\xi}(\eta^2-1)^2.
$$

**罚项里的约束泛函**：
$$
V=\int_\Omega\phi\ \ (\text{体积}),\qquad
A=\int_\Omega W_\phi\ \ (\text{总面积}),\qquad
D=\int_\Omega\tanh\tfrac\eta\xi\,W_\phi\ \ (\text{面积差}),
$$
$$
N=\int_\Omega\tfrac\varepsilon2|\nabla\phi\cdot\nabla\eta|^2\ \ (\text{正交}),\qquad
P=\int_\Omega p^2,\quad p:=\tfrac\xi2|\nabla\eta|^2-\tfrac1{4\xi}(\eta^2-1)^2\ \ (\eta\ \text{廓形}).
$$
（$A$ 复用 $W_\phi$；$D$ 是 $W_\phi$ 的 $\tanh$-加权；$P$ 是 $\eta$ 廓形泛函的平方。）

**阶数小结**：四阶主部只出现在 $\phi$ 方程、且只来自 $E$（变分得 $\varepsilon\Delta(k(\eta)\Delta\phi)$）；$\eta$ 方程最高 2 阶。**唯一结构性难点 = $\phi$ 的四阶变系数主部**——§3 先攻它。各项变分导数在 Step 4 拼装时逐项给出。

---

## 2. 时间离散主干

BDF2（二阶）。时间差商与二阶外推：
$$
\frac{3u^{n+1}-4u^n+u^{n-1}}{2\Delta t},\qquad u^\star:=2u^n-u^{n-1}.
$$
（$n{=}0{\to}1$ 用一步一阶自举。）算子 $\Delta,\Delta^2$ **保持连续形式**——本文只做到时间离散的"数值格式"，不部署空间离散。

---

## 3. 子问题：变系数四阶主部 → 时间更新格式

**子问题**：$\eta$ 取已知外推值 $\eta^\star$，记 $a:=k(\eta^\star)\in[k_2,k_1]$（已知正系数场）。弯曲主部的 $L^2$ 梯度流
$$
u_t=-\frac{\delta E_b}{\delta u},\qquad
E_b[u]=\frac{\varepsilon}{2}\int_\Omega a\,|\Delta u|^2,\qquad
\frac{\delta E_b}{\delta u}=\varepsilon\,\Delta\!\big(a\,\Delta u\big).
$$
$a$ 夹在两个 $\Delta$ 之间：整体隐式即变系数四阶算子，不可逐波求逆。

**IMEX 劈分**（推导与耗散证明见 [IMEX.md](IMEX.md)）：取中点 $\bar k=\tfrac12(k_1+k_2)$，拆
$$
\varepsilon\Delta(a\Delta u)=\underbrace{\varepsilon\bar k\,\Delta^2 u}_{\text{常系数,隐式}}+\underbrace{\varepsilon\Delta\big((a-\bar k)\Delta u\big)}_{\text{偏差,显式外推 }u^\star},\qquad
|a-\bar k|\le c:=\tfrac12(k_1-k_2),
$$
并加同阶稳定化 $\varepsilon S\,\Delta^2(u^{n+1}-u^\star)$（因 $u^{n+1}-u^\star=O(\Delta t^2)$，不损二阶精度）。

**时间更新格式**（BDF2，时间差商 = 右端算子）：
$$
\boxed{\;\frac{3u^{n+1}-4u^n+u^{n-1}}{2\Delta t}
=-\,\varepsilon\bar k\,\Delta^2 u^{n+1}
-\,\varepsilon S\,\Delta^2\big(u^{n+1}-u^\star\big)
-\,\varepsilon\,\Delta\!\big((a-\bar k)\,\Delta u^\star\big)\;}
$$
等价地（$u^{n+1}$ 在左）：$\big(\tfrac{3}{2\Delta t}+\varepsilon(\bar k+S)\Delta^2\big)u^{n+1}=\text{已知}$，左端常系数 $\Delta^2$ 算子，空间离散后逐波求逆（第 ④ 步）。

**参数取值**：
$$
\bar k=\tfrac12(k_1+k_2)\ \ (\text{中点,使 }|a-\bar k|\le c\text{ 最小}),\qquad
S\ge c=\tfrac12(k_1-k_2)\ \ (h,\varepsilon\text{ 无关}).
$$
下降量为修正能量 $\tilde E=\int\tfrac{\bar k}{2}\varepsilon|\Delta u^{n+1}|^2+\text{(BDF2 增量罚)}$；与真能量之差 $\int\tfrac{a-\bar k}{2}\varepsilon|\Delta u|^2$ 显式可算。完整 BDF2 耗散证明（$S$ 精确下界 + telescoping）在拼装定稿后单列。

---

## 4. 子问题：非线性势 → SAV（标量更新 + 线性场更新）

§3 解决了最高阶的刚性线性骨架，但弯曲能 $G^2$ 展开后还剩非线性势（$\phi^3,\phi^6$），$L,A,D,N,P$ 也都非线性。最小模型：在刚性线性部分上加一个非线性势
$$
u_t=-\Big(\varepsilon\bar k\,\Delta^2 u+f(u)\Big),\qquad f=F'(u),\qquad
\mathcal F[u]:=\int_\Omega F(u),
$$
$\mathcal F$ 下有界 ⟹ 存在常数 $C_0$ 使 $\mathcal F[u]+C_0\ge\delta>0$。难点：$f$ 非线性——全隐式要 Newton 迭代，全显式无能量控制。

**SAV 重写**：唯一新增的辅助变量是标量
$$
r:=\sqrt{\mathcal F[u]+C_0}\quad\Longrightarrow\quad f(u)=\frac{f(u)}{\sqrt{\mathcal F[u]+C_0}}\,r .
$$
梯度流等价改写为 $u$–$r$ 耦合系统（非线性因子 $\dfrac{f(u)}{\sqrt{\mathcal F[u]+C_0}}$ 始终显式写出，不另设符号）：
$$
u_t=-\,\varepsilon\bar k\,\Delta^2 u-\frac{f(u)}{\sqrt{\mathcal F[u]+C_0}}\,r,\qquad
r_t=\frac12\Big(\frac{f(u)}{\sqrt{\mathcal F[u]+C_0}},\ u_t\Big).
$$
第一式与 $\varepsilon\bar k\Delta^2u$ 项配对、第二式乘 $2r$，相加得**修正能量**耗散律
$$
\frac{d}{dt}\Big(\tfrac{\varepsilon\bar k}{2}\|\Delta u\|^2+r^2\Big)=-\|u_t\|^2\le0,
$$
其中 $r^2=\mathcal F[u]+C_0$，故修正能量 = 真能量 + 常数。

**时间更新格式**（BDF2，非线性因子显式外推到 $u^\star=2u^n-u^{n-1}$）：
$$
\boxed{\;
\begin{aligned}
\frac{3u^{n+1}-4u^n+u^{n-1}}{2\Delta t}&=-\,\varepsilon\bar k\,\Delta^2 u^{n+1}
-\frac{f(u^\star)}{\sqrt{\mathcal F[u^\star]+C_0}}\,r^{n+1},\\[4pt]
\frac{3r^{n+1}-4r^n+r^{n-1}}{2\Delta t}&=\frac12\Big(\frac{f(u^\star)}{\sqrt{\mathcal F[u^\star]+C_0}},\ \frac{3u^{n+1}-4u^n+u^{n-1}}{2\Delta t}\Big).
\end{aligned}\;}
$$
未知量是场 $u^{n+1}$ + **标量** $r^{n+1}$；非线性只进已知的 $\dfrac{f(u^\star)}{\sqrt{\mathcal F[u^\star]+C_0}}$，故场方程对 $u^{n+1}$ **线性、常系数**。

**解耦（标量更新 + 线性场更新）**：把第一式整理成"$u^{n+1}$ 在左"，左端是常系数算子（逐波可逆）：
$$
\Big(\tfrac{3}{2\Delta t}+\varepsilon\bar k\Delta^2\Big)u^{n+1}
=\frac{4u^n-u^{n-1}}{2\Delta t}-\frac{f(u^\star)}{\sqrt{\mathcal F[u^\star]+C_0}}\,r^{n+1}.
$$
右端对标量 $r^{n+1}$ 仿射，故 $u^{n+1}$ 也对 $r^{n+1}$ 仿射：
$$
u^{n+1}=\Big(\tfrac{3}{2\Delta t}+\varepsilon\bar k\Delta^2\Big)^{-1}\frac{4u^n-u^{n-1}}{2\Delta t}
-r^{n+1}\Big(\tfrac{3}{2\Delta t}+\varepsilon\bar k\Delta^2\Big)^{-1}\frac{f(u^\star)}{\sqrt{\mathcal F[u^\star]+C_0}}.
$$
代入第二式（乘 $2\Delta t$ 后 $3u^{n+1}=3u^{n+1}$），$r^{n+1}$ 由一个标量方程解出：
$$
r^{n+1}=\frac{\;4r^n-r^{n-1}
+\dfrac32\Big(\dfrac{f(u^\star)}{\sqrt{\mathcal F[u^\star]+C_0}},\ \big(\tfrac{3}{2\Delta t}+\varepsilon\bar k\Delta^2\big)^{-1}\dfrac{4u^n-u^{n-1}}{2\Delta t}\Big)
-\dfrac12\Big(\dfrac{f(u^\star)}{\sqrt{\mathcal F[u^\star]+C_0}},\ 4u^n-u^{n-1}\Big)\;}
{\;3+\dfrac32\Big(\dfrac{f(u^\star)}{\sqrt{\mathcal F[u^\star]+C_0}},\ \big(\tfrac{3}{2\Delta t}+\varepsilon\bar k\Delta^2\big)^{-1}\dfrac{f(u^\star)}{\sqrt{\mathcal F[u^\star]+C_0}}\Big)\;}.
$$
分母 $\ge3>0$（常系数算子 $\tfrac{3}{2\Delta t}+\varepsilon\bar k\Delta^2$ 正定 ⟹ 其逆的内积 $\ge0$）——**恒可解，无 Newton 迭代**。回代上式得 $u^{n+1}$。

下降量为修正能量（$G$-范数二次能量 + $r^2$），无条件单调；完整 BDF2 能量恒等式留到拼装后单列。

> **拼装时的并入**：完整模型左端实际算子是 §3 的 $\tfrac{3}{2\Delta t}I+\varepsilon(\bar k+S)\Delta^2$（含稳定化），SAV 机制原样骑在其上；多个非线性泛函 $L,A,D,N,P$ 各配一个标量（MSAV），$V$ 线性 → 精确隐式，无需辅助变量。

---

> **下一步（待本格式经审核后再写）**
> - **Step 4**：把弯曲（§3）、非线性势（§4）、线张力 $L$、罚项 $A,D,N,P$、线性 $V$ 逐项拼装成完整模型的 **φ-更新** 与 **η-更新**（MSAV：每个非线性泛函一个标量）。
> - 之后：离散耗散证明 → 空间离散（谱 / FFT）→ 高效求解 → 自适应步长 / 诊断 / Phase-B。
