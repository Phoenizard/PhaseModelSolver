# Wang–Du 双组分膜模型 — 数值方案（部署文档）

> **文档定位**：本文给出 2D 基线求解器的**完整、可实现**数值方案（"具体怎么部署"）。它是 [doc/debate/](../debate/) 五方三轮去偏见辩论的收敛产物（见 [synthesis 收敛点](../debate/round3.md)）。
> **状态**：算法设计定稿（开发链第 1 步）。耗散证明（第 2 步）与数值验收（第 4 步）的待办在 §11–§12 列出。
> **记号约定**：模型中线张力常数原记 $\delta$，与变分记号 $\delta/\delta\phi$ 冲突，本文线张力记 **$\delta_L$**。$d\in\{2,3\}$，先 2D。

---

## 0. 一页速览

| 能量项 | 阶数 / 变系数位置 | 处理机制 | 左端贡献 | 耗散类型 |
|---|---|---|---|---|
| 弯曲 $E$ — 四阶主部 $\varepsilon\,\Delta(k(\eta)\Delta\phi)$ | 4 阶，$k(\eta)$ 在主部 | **冻结 $\bar k$ 隐式 + 同阶 $\Delta^2$ 稳定化**；余项显式外推 | $\varepsilon(\bar k{+}S)\Delta^2$（常系数） | 代理能量，无条件、$h$ 无关 |
| 弯曲 $E$ — 低阶非线性 + $\delta E/\delta\eta$ | $\le$2 阶 | 显式外推 + $\beta(-\Delta){+}\gamma$ 稳定化 | $\beta(-\Delta){+}\gamma$ | 修正能量，条件（见 §11） |
| 线张力 $L$ | 2 阶，变系数 | SAV（$r_L=\sqrt{L{+}C_L}$） | 秩一更新 | 修正能量，无条件 |
| 体积 $V$（线性泛函） | 0 阶 | **精确隐式**（秩一，非 SAV） | 秩一更新 | 精确 |
| 面积 $A$ / 面积差 $D$ / 正交 $N$ / 廓形 $P$ 罚项 | 泛函的平方 | **线性辅助变量 MSAV**（$q_\bullet=$ 泛函值） | 秩一更新 | 修正能量，无条件 |

**主干**：BDF2（二阶）+ Fourier 伪谱 + 顺序解耦 $\phi\!\to\!\eta$。两个子步的左端算子都是**常系数、Fourier 逐模对角**，FFT 直接求逆 + Woodbury 处理秩一耦合。

**核心设计依据（辩论结论）**：硬约束3 强制 $k(\eta)$ 冻结外推 ⟹ 对真 $E_M$ 的无条件原能量谁都拿不到。最优分工 = **唯一的四阶变系数主部用冻结+同阶 $\Delta^2$ 稳定化（凸分裂式，$S$ 与 $h$ 无关）**，**其余低阶/罚项用 MSAV（吸收大 $M_i$、无需 Lipschitz 估计）**。原能量是近平衡 **Phase-B** 的可选忠实化（§10），不作基线主线。

---

## 1. 连续问题

域 $\Omega=[-\pi,\pi]^d$，周期边界。$L^2$ 梯度流
$$
\phi_t=-\frac{\delta E_M}{\delta\phi},\qquad \eta_t=-\frac{\delta E_M}{\delta\eta},\qquad
E_M=E+L+\tfrac12 M_1(V-v_d)^2+\tfrac12 M_2(A-a_0)^2+\tfrac12 M_3(D-a_d)^2+\tfrac12 M_4 N^2+\tfrac12 M_5 P^2.
$$

密度与泛函（$c_1=c_2=0\Rightarrow c_0\equiv0$，$\xi=\varepsilon$）：
$$
W_\phi=\tfrac\varepsilon2|\nabla\phi|^2+\tfrac1{4\varepsilon}(\phi^2-1)^2,\quad
W_\eta=\tfrac\xi2|\nabla\eta|^2+\tfrac1{4\xi}(\eta^2-1)^2,\quad
G=\varepsilon\Delta\phi+\tfrac1\varepsilon(\phi-\phi^3),
$$
$$
A=\!\int W_\phi,\ V=\!\int\phi,\ D=\!\int\tanh\tfrac\eta\xi\,W_\phi,\ N=\!\int\tfrac\varepsilon2|\nabla\phi\cdot\nabla\eta|^2,\ P=\!\int p^2\ (p:=\tfrac\xi2|\nabla\eta|^2-\tfrac1{4\xi}(\eta^2-1)^2),
$$
$$
E=\!\int\tfrac{k(\eta)}{2\varepsilon}G^2,\quad L=\delta_L\!\int W_\eta W_\phi,\quad k(\eta)=k+c\tanh\tfrac\eta\xi,\quad k'(\eta)=\tfrac c\xi\,\mathrm{sech}^2\tfrac\eta\xi.
$$

### 变分导数（实现 `energy.py` 的清单）

$$
\frac{\delta A}{\delta\phi}=-\varepsilon\Delta\phi+\tfrac1\varepsilon(\phi^3-\phi),\qquad \frac{\delta V}{\delta\phi}=1,
$$
$$
\boxed{\frac{\delta E}{\delta\phi}=\Delta\!\big(k(\eta)G\big)+\tfrac1{\varepsilon^2}(1-3\phi^2)\,k(\eta)G},\qquad
\frac{\delta E}{\delta\eta}=\frac{k'(\eta)}{2\varepsilon}G^2,
$$
四阶主部 $=\varepsilon\,\Delta\!\big(k(\eta)\Delta\phi\big)$；$\delta E/\delta\eta$ **零阶（乘性）**于 $\eta$（因 $c_0=0$，$G$ 不含 $\eta$）。
$$
\frac{\delta D}{\delta\phi}=-\varepsilon\nabla\!\cdot\!\big(\tanh\tfrac\eta\xi\nabla\phi\big)+\tfrac1\varepsilon\tanh\tfrac\eta\xi(\phi^3-\phi),\qquad
\frac{\delta D}{\delta\eta}=\tfrac1\xi\mathrm{sech}^2\tfrac\eta\xi\,W_\phi,
$$
$$
\frac{\delta N}{\delta\phi}=-\varepsilon\nabla\!\cdot\!\big((\nabla\phi\cdot\nabla\eta)\nabla\eta\big),\qquad
\frac{\delta N}{\delta\eta}=-\varepsilon\nabla\!\cdot\!\big((\nabla\phi\cdot\nabla\eta)\nabla\phi\big),
$$
$$
\frac{\delta P}{\delta\phi}=0,\qquad
\frac{\delta P}{\delta\eta}=-2\xi\nabla\!\cdot\!(p\,\nabla\eta)-\tfrac2\xi(\eta^2-1)\eta\,p,
$$
$$
\frac{\delta L}{\delta\phi}=-\varepsilon\nabla\!\cdot\!(\delta_L W_\eta\nabla\phi)+\tfrac1\varepsilon\delta_L W_\eta(\phi^3-\phi),\qquad
\frac{\delta L}{\delta\eta}=-\xi\nabla\!\cdot\!(\delta_L W_\phi\nabla\eta)+\tfrac1\xi\delta_L W_\phi(\eta^3-\eta).
$$

**阶数 / 变系数小结**：$\phi$ 方程最高 4 阶（仅来自 $E$，变系数 $k(\eta)$）；$\eta$ 方程最高 2 阶（来自 $L,D,N,P$，变系数依赖 $\phi$ 与 $\eta$）。**没有 $\eta$ 的四阶主部** —— $\eta$ 子步显著更便宜。

---

## 2. 空间离散

Fourier 伪谱。$N_g^d$ 网格，波数 $\mathbf K$。导数在谱空间（$\widehat{\partial_j u}=i k_j\hat u$，$\widehat{-\Delta u}=|\mathbf k|^2\hat u$，$\widehat{\Delta^2 u}=|\mathbf k|^4\hat u$），非线性密度在实空间 Fourier 节点求值（伪谱），乘积用 FFT 往返。**抗混叠**：弯曲 $G^2$、罚项 $(\cdot)^2$ 等可达**六次**多项式；3/2 规则只精确去**二次**乘积的混叠，六次需 **2 倍补零**（pad 到 $\ge 2N_g$）或取足够大 $N_g$。求解器对 $d$ 参数化（同一套代码 2D/3D）。

---

## 3. 时间离散主干

BDF2 + 二阶外推 + 顺序解耦。记
$$
D_t u^{n+1}=\frac{3u^{n+1}-4u^n+u^{n-1}}{2\Delta t},\qquad u^\star:=2u^n-u^{n-1}\ (\text{二阶外推}).
$$
**顺序**：先解 $\phi^{n+1}$（$\eta$ 取 $\eta^\star$），再解 $\eta^{n+1}$（$\phi$ 取**新值** $\phi^{n+1}$）。起步用一步一阶（BDF1/稳定化 Euler）自举 $n=0\to1$。

---

## 4. $\phi$ 子步

### 4.1 半隐分解

$$
D_t\phi^{n+1}=-\Big[\;\underbrace{\varepsilon\bar k\Delta^2\phi^{n+1}}_{\text{冻结主部(隐式)}}
+\underbrace{\varepsilon S\Delta^2(\phi^{n+1}-\phi^\star)}_{\text{同阶稳定(隐式-外推)}}
+\underbrace{\big(\beta(-\Delta)+\gamma\big)(\phi^{n+1}-\phi^\star)}_{\text{低阶稳定}}
+\underbrace{R_\phi^\star}_{\text{弯曲显式余项}}
+\underbrace{r_L^{n+1}u_{L,\phi}^\star}_{L\text{ 的SAV}}
+\underbrace{M_1(V^{n+1}-v_d)}_{V\text{ 精确隐式}}
+\underbrace{\textstyle\sum_{Q\in\{A,D,N\}}\!c_Q^{n+1}\,g_{Q,\phi}^\star}_{\text{罚项 MSAV}}\;\Big],
$$
其中
- $\bar k=\tfrac12(k_1+k_2)$（居中冻结，使余幅 $|k(\eta)-\bar k|\le|c|$ 最小），$c:=\tfrac12(k_1-k_2)$；稳定常数 **$S\ge|c|$**（取 $S=|c|$ 为紧；**与 $h,\varepsilon$ 无关**，依据见 §11.1）；
- $R_\phi^\star:=\big(\tfrac{\delta E}{\delta\phi}\big|_{(\phi^\star,\eta^\star)}-\varepsilon\bar k\Delta^2\phi^\star\big)$ —— 弯曲变分导数减去冻结主部，**全部用外推态显式求值**（含变系数四阶余项 $\varepsilon\Delta((k(\eta^\star)-\bar k)\Delta\phi^\star)$ 与低阶非线性）；
- $u_{L,\phi}^\star:=\dfrac{1}{\sqrt{L+C_L}}\dfrac{\delta L}{\delta\phi}\Big|_{(\phi^\star,\eta^\star)}$，$r_L$ 为线张力 SAV（§6）；
- $V^{n+1}=\int\phi^{n+1}$（线性，精确）；
- $c_A^{n+1}:=M_2(q_A^{n+1}-a_0)$，$c_D^{n+1}:=M_3(q_D^{n+1}-a_d)$，$c_N^{n+1}:=M_4 q_N^{n+1}$；$q_\bullet$ 为线性辅助变量（§6），$g_{Q,\phi}^\star:=\frac{\delta Q}{\delta\phi}\big|_{(\phi^\star,\eta^\star)}$。

### 4.2 左端算子（常系数，Fourier 对角）

$$
\boxed{\ \widehat{\mathcal A_\phi}(\mathbf k)=\frac{3}{2\Delta t}+\varepsilon(\bar k+S)\,|\mathbf k|^4+\beta\,|\mathbf k|^2+\gamma\ >0\quad\forall\mathbf k\ }
$$
$k(\eta)$ 完全不在左端。逐波数标量除法即得 $\mathcal A_\phi^{-1}$。

### 4.3 秩一耦合的消去（Woodbury）

$L,A,D,N$ 各自的 SAV/辅助变量使 $\phi^{n+1}$ 方程出现"$\mathcal A_\phi\phi^{n+1}=\text{RHS}+\sum_i \mathbf a_i\,\langle\mathbf b_i,\phi^{n+1}\rangle$"的**秩一更新和**（$V$ 亦然）。令 $m$=耦合数（此处 $\le5$）。解法：
1. 对 $\mathcal A_\phi$ 解 $m+1$ 个右端（FFT 各一次）：$\psi_0=\mathcal A_\phi^{-1}\text{RHS}$，$\psi_i=\mathcal A_\phi^{-1}\mathbf a_i$。
2. 解 $m\times m$ 标量线性系统得各 $\langle\mathbf b_i,\phi^{n+1}\rangle$。
3. $\phi^{n+1}=\psi_0+\sum_i\langle\cdots\rangle_i\psi_i$。
每步 $\phi$ 子步成本 $=O(m)$ 次 FFT 正逆 + $O(m^3)$ 标量代数（$m\le5$ 可忽略）。**无非线性迭代。**

---

## 5. $\eta$ 子步

拿到 $\phi^{n+1}$ 后（在所有 $\phi$-相关量中用之）。$\eta$ 最高 2 阶，无四阶主部。
$$
D_t\eta^{n+1}=-\Big[\big(\alpha(-\Delta)+\gamma_\eta\big)(\eta^{n+1}-\eta^\star)
+R_\eta^\star
+r_L^{n+1}u_{L,\eta}^\star
+\textstyle\sum_{Q\in\{D,N,P\}}c_Q^{n+1}\,g_{Q,\eta}^\star\Big],
$$
- $R_\eta^\star:=\frac{\delta E}{\delta\eta}\big|_{(\phi^{n+1},\eta^\star)}=\frac{k'(\eta^\star)}{2\varepsilon}G(\phi^{n+1})^2$（零阶，显式）；
- $g_{Q,\eta}^\star=\frac{\delta Q}{\delta\eta}\big|_{(\phi^{n+1},\eta^\star)}$；$c_P^{n+1}:=M_5 q_P^{n+1}$（此子步 $\delta P/\delta\eta\ne0$）；
- 左端 $\widehat{\mathcal A_\eta}(\mathbf k)=\frac{3}{2\Delta t}+\alpha|\mathbf k|^2+\gamma_\eta>0$，常系数对角。
- $L,D,N,P$ 的 SAV/辅助标量在本子步**补上 $\eta$ 部分的时间导数贡献**（见 §6，交叉项闭合的关键），同样 Woodbury 消秩一。

稳定常数 $\alpha$ 须盖住 $\eta$ 的变系数 2 阶项幅度（$\xi\delta_L\sup W_\phi$、$N$、$P$ 的 2 阶系数）；$\gamma_\eta$ 盖住零阶非线性 Lipschitz。具体下界在 §11 证明中定。

---

## 6. 辅助变量定义与更新

**线张力（经典开方 SAV）**：$r_L=\sqrt{L+C_L}$，$C_L$ 取使 $L+C_L>0$。
$$
D_t r_L^{n+1}=\tfrac12\!\int u_{L,\phi}^\star D_t\phi^{n+1}+\tfrac12\!\int u_{L,\eta}^\star D_t\eta^{n+1}\quad(\text{两子步各贡献一半}).
$$

**罚泛函（线性辅助变量，无开方、无下界）**：对 $Q\in\{A,D,N,P\}$ 令 $q_Q\approx Q$，
$$
D_t q_Q^{n+1}=\int g_{Q,\phi}^\star\,D_t\phi^{n+1}+\int g_{Q,\eta}^\star\,D_t\eta^{n+1}.
$$
（$A$ 仅 $\phi$ 部分；$P$ 仅 $\eta$ 部分；$D,N$ 两部分。）罚项能量写 $\tfrac12 M_Q(q_Q-\text{target})^2$。

**体积 $V$（线性，精确）**：$V^{n+1}=\int\phi^{n+1}$，不需辅助变量，直接作精确秩一隐式。

> **交叉项闭合（$N$ 是关键，对应辩论 SAV-Q4）**：$q_N$ 是**单个标量**，其 $D_t q_N^{n+1}$ 在 $\phi$ 子步收 $\int g_{N,\phi}^\star D_t\phi$、在 $\eta$ 子步收 $\int g_{N,\eta}^\star D_t\eta$；与两子步场方程中的 $-M_4 q_N^{n+1}g_{N,\bullet}^\star$ 逐项配对，能量恒等式里交叉项 $\phi$-$\eta$ 精确抵消。$L,D$ 同理。**这是顺序解耦下交叉项不漏的机制。**
>
> **实现/证明要点（§11.3 的真正风险）**：精确抵消要求两子步的力系数与 $q_\bullet$ 更新**用同一终值 $q_\bullet^{n+1}$**。但顺序解耦下解 $\phi^{n+1}$ 时 $q_\bullet^{n+1}$ 的 $\eta$ 贡献尚未知。可行构造：$\phi$ 子步对 $\eta$ 贡献用外推、对 $\phi$ 贡献 Woodbury 联立解 $(\phi^{n+1},\{q\})$；$\eta$ 子步以**同一** $q_\bullet^{n+1}$ 联立补 $\eta$ 贡献并回代。须证此构造仍精确抵消（否则退化半隐、留残差）。**不可**用半更新的 $q$ 充当终值。

---

## 7. 每步算法（伪代码）

```
输入: phi^n, phi^{n-1}, eta^n, eta^{n-1}, 辅助 {r_L,q_A,q_D,q_N,q_P}^{n,n-1}, dt
1. 外推: phi*=2phi^n-phi^{n-1}; eta*=2eta^n-eta^{n-1}
2. 实空间装配 (用 phi*,eta*): R_phi*, g_{A,phi}*,g_{D,phi}*,g_{N,phi}*, u_{L,phi}*
3. phi 子步: 组 RHS_phi(BDF2历史 + 稳定外推项 - 显式项); Woodbury 解 phi^{n+1}
            联立解 r_L,q_A,q_D,q_N 的 phi-贡献 (eta-贡献用外推; 终值规则见 §6 要点)
4. 实空间装配 (用 phi^{n+1},eta*): R_eta*, g_{D,eta}*,g_{N,eta}*,g_{P,eta}*, u_{L,eta}*
5. eta 子步: Woodbury 解 eta^{n+1}; 以同一终值 q^{n+1} 联立补 r_L,q_D,q_N,q_P 的 eta-贡献 (见 §6)
6. 诊断: 计算真 E_M^{n+1}, 修正能量 \tilde E^{n+1}, 间隙; 检查单调性
7. 自适应 dt (§9); 若违反则缩 dt 重做本步
输出: ^{n+1} 全量
```

---

## 8. 离散能量律与耗散诊断（**一等公民**）

**被证明单调的是修正/代理能量。** 在 BDF2 下，严格单调的对象是下式能量的 **G-范数提升**（精确恒等式见 §11.1）；其一致主部、亦即实际监控量为
$$
\widetilde E^{n+1}:=\underbrace{\int\tfrac{\bar k}{2\varepsilon}G(\phi^{n+1})^2}_{\text{弯曲代理：真 }k\text{ 换常数 }\bar k\text{，不含 }\eta}
+\tfrac{\varepsilon S}{2}\|\Delta(\phi^{n+1}-\phi^n)\|^2
+\tfrac12 (r_L^{n+1})^2-C_L
+\tfrac12 M_1(V^{n+1}-v_d)^2+\!\!\sum_{Q}\tfrac12 M_Q(q_Q^{n+1}-\text{target})^2.
$$
弯曲代理仅以常数 $\bar k$ 加权（与诊断 3 间隙 $\int\tfrac{\bar k-k}{2\varepsilon}G^2$ 自洽，**不重复计数**）。**BDF2 提升**把每个二次项 $\tfrac12\|\cdot^{n+1}\|^2$ 换成 G-范数 $\tfrac14(\|\cdot^{n+1}\|^2+\|2\cdot^{n+1}-\cdot^n\|^2)$；§4.1 的 $\phi^\star$（二阶差）稳定项经 §11.1 的恒等式恰好遗留上式的 $\tfrac{\varepsilon S}2\|\Delta(\phi^{n+1}-\phi^n)\|^2$（故 §4.1 与本式**相容**），单步耗散另含非负项 $\varepsilon S\|\Delta(\phi^{n+1}-2\phi^n+\phi^{n-1})\|^2\ge0$。目标：G-范数 $\widetilde E^{n+1}\le\widetilde E^n$（条件见 §11）。

**诊断量（`diagnostics.py`，每步记录）**：
1. 真能量 $E_M^{n+1}$（用**真实** $k(\eta^{n+1})$）—— 主验收对象。
2. 修正能量 $\widetilde E^{n+1}$ 与单调性 $\widetilde E^{n+1}-\widetilde E^n$。
3. **间隙** $\widetilde E-E_M$ 与各分量漂移：弯曲代理间隙 $\int\tfrac{\bar k-k(\eta)}{2\varepsilon}G^2$（显式可算）、SAV 漂移 $\tfrac12 r_L^2-(L+C_L)$、$q_Q-Q$。
4. 约束残差 $|V-v_d|,|A-a_0|,|D-a_d|,N,P$。

---

## 9. 自适应 $\Delta t$

以**能量下降**为准则（论文与项目惯例）：
- 若 $E_M^{n+1}>E_M^n+\text{tol}$（真能量上升超容差）⟹ $\Delta t\leftarrow\rho^-\Delta t$（$\rho^-=0.5$）重做。
- 若连续 $K$ 步平滑下降 ⟹ $\Delta t\leftarrow\min(\rho^+\Delta t,\Delta t_{\max})$（$\rho^+=1.2$）。
- 同时受稳定性/精度上限约束。BDF2 变步长用变系数公式或每变步重启二阶自举。

---

## 10. Phase-B 可选 add-on（基线走通后）

以下不进基线，作升级，对应 Notion 阶梯：

- **(A) 加权 SAV–Lagrange 的 $\omega$ 调度（近平衡原能量忠实化）**：仅对弯曲 $E$ 引权 $\omega_E$。瞬态 $\omega_E=1$（纯 SAV/本基线，大步）；近平衡（能量下降率 / 残差 $\|\mu\|$ 低于阈）切 $\omega_E\downarrow0$，在小步区开标量乘子 $\lambda$ 把弯曲冻结漂移折回真 $k$ 度量；某步 $g(\lambda)=0$ 无近 1 正根则 $\omega_E\uparrow1$ 回退，耗散不破防。**风险**：$g(\lambda)$ 变系数可解性无先例（§12）。
- **(B) ETD 刚部精确积分器**：仅当 $|c|\ll\bar k$ 且 $\varepsilon$ 小、$\bar k$ 大使线性刚部极硬时，用 $e^{\widehat{\mathcal L}\Delta t}$ 精确积分常系数刚部 $-\varepsilon\bar k\Delta^2$（$M_1,M_2$ 标量罚项可一并冻进对角 $\widehat{\mathcal L}=-\varepsilon\bar k|\mathbf k|^4-\tau|\mathbf k|^2-\sigma^n\varepsilon|\mathbf k|^2$），其余维持本基线。优势在大 $\Delta t$ 下刚模精确不失相。
- **(C) IEQ 逐点修补**：仅 $D$/$L$ 耦合的 $\mathrm{sech}^2(\eta/\xi)W_\phi$ 逐点权重项，若标量 MSAV 在该项交叉闭合数值上不足时启用。多数情况不需要。

---

## 11. 待证清单（开发链第 2 步 — 与设计构成循环）

1. **弯曲代理无条件性（BDF2 G-范数 + 变分型估计）** —— 三件事须写清：
   - **BDF2 恒等式**：$(3a{-}4b{+}c,a)=\tfrac12(\|a\|^2+\|2a{-}b\|^2-\|b\|^2-\|2b{-}c\|^2+\|a{-}2b{+}c\|^2)$（$a{=}u^{n+1},b{=}u^n,c{=}u^{n-1}$），给出隐式二次项以 $D_t u^{n+1}$ 测试时的 G-范数遗留 + 非负余项。
   - **$\phi^\star$ 稳定相容性**（回应"$\phi^\star$ vs $\phi^n$"）：$(\Delta^2(\phi^{n+1}{-}\phi^\star),D_t\phi^{n+1})=\tfrac1{\Delta t}\|\Delta\delta^2\phi^{n+1}\|^2+\tfrac1{2\Delta t}(\|\Delta(\phi^{n+1}{-}\phi^n)\|^2-\|\Delta(\phi^n{-}\phi^{n-1})\|^2)$，$\delta^2\phi^{n+1}{=}\phi^{n+1}{-}\phi^\star$。故二阶外推稳定（保 $O(\Delta t^2)$）与 §8 遗留项 $\tfrac{\varepsilon S}2\|\Delta(\phi^{n+1}{-}\phi^n)\|^2$ **同时成立**，外加非负耗散 $\varepsilon S\|\Delta\delta^2\phi^{n+1}\|^2$。
   - **变系数余项用变分型（非逐项展开）估计**：余项 $\varepsilon\Delta((k(\eta^\star){-}\bar k)\Delta\phi^\star)$ 在能量估计中经两次分部积分**保持为二次型** $\varepsilon\int(k(\eta^\star){-}\bar k)\,\Delta\phi^\star\,\Delta(\cdot)$，故 $S\ge|c|$ 即控、**与 $h,\varepsilon$ 均无关**。**切勿**展开成 $(k{-}\bar k)\Delta^2\phi+2\nabla k\!\cdot\!\nabla\Delta\phi+\Delta k\,\Delta\phi$ 逐项 Fourier 估计——那会引出虚假三阶交叉项（系数 $\sim c/\varepsilon$）与过悲观的 $\beta\sim c/\varepsilon^2$ 需求；变分型估计中 $\nabla k$ 永不单独露面。
2. **低阶稳定常数 $\beta,\gamma,\alpha,\gamma_\eta$ 的下界与 $\varepsilon$-标度**：盖住**真正低阶**的非线性项——弯曲的 $\Delta\!\big(k\cdot\tfrac1\varepsilon(\phi{-}\phi^3)\big)$ 与 $\tfrac1{\varepsilon^2}(1{-}3\phi^2)kG$（注意 $G\sim1/\varepsilon$，量级达 $1/\varepsilon^2\text{–}1/\varepsilon^3$，$\beta,\gamma$ 的 $\varepsilon$-标度据此定，**不可当 $O(1)$**）、以及 $L,D,N,P$ 变系数 2 阶/零阶项的 Lipschitz。**备选**：若这些 $\varepsilon$-病态低阶项的稳定常数过大损精度，可把弯曲低阶非线性也纳入 SAV（$r_E$ 仅作低阶余项，主部仍冻结+稳定化），以无条件性替代 Lipschitz 估计——留作 step-2 设计循环的分支。
3. **MSAV 子部无条件性 + 交叉项精确消去**：$r_L,q_A,q_D,q_N,q_P$ 的二次型与场方程配对，$\phi$-$\eta$ 交叉项（$L,D,N$）在顺序解耦下严格抵消（§6 机制的严格化）。
4. **整体修正能量单调** $\widetilde E^{n+1}\le\widetilde E^n$ 的合成条件。
5. **二阶时间精度**：稳定项 $\propto(\cdot^{n+1}-\cdot^\star)=O(\Delta t^2)$、外推 $O(\Delta t^2)$、BDF2 + SAV 同阶。
6. 复用文献框架：Shen–Xu–Yang(SAV)、Cheng–Shen(MSAV 囊泡 SISC2018)、Guan–Lowengrub–Wang–Wise(二阶凸分裂/同阶稳定化)。**注**：$k(\eta)$ 四阶主部 + 二阶 + 本方法的精确组合**无直接先例**，第 1 条是需自建的核心引理。

## 12. 待验清单（开发链第 4 步，2D）

- 离散耗散：长时单调性、漂移累积、$\Delta t$ 稳健性、退化构型。**真 $E_M$ 与修正 $\widetilde E$ 双轨记录。**
- 二阶收敛：时间细化下 $\phi,\eta$ 的 $O(\Delta t^2)$（Cauchy 或对参考解）。
- 稳定常数标定：$S=|c|$ 是否充分；$\beta,\gamma,\alpha,\gamma_\eta$ 实测下界与 $\varepsilon$-标度；刚度对比 $|c|/\bar k$ 大时的步长行为（是否触发 $h^4$ 类约束 → 决定是否上 Phase-B(B)）。
- 约束达成：$M_i$ 增大序列下 $V\to v_d,A\to a_0,D\to a_d$。
- （3D 阶段）形态复现 + 维度独立的耗散再验。

---

## 13. 代码模块映射（`src/`）

| 模块 | 职责 |
|---|---|
| `energy.py` | §1 全部变分导数、密度、$k(\eta),k'(\eta)$、真 $E_M$ 求值 |
| `operators.py` | 谱算子、$\widehat{\mathcal A_\phi},\widehat{\mathcal A_\eta}$ 的对角符号与逆、FFT 往返、去混叠 |
| `scheme.py` | §3–§7 BDF2 顺序解耦、Woodbury、辅助变量更新、自举 |
| `diagnostics.py` | §8 能量双轨、间隙、约束残差、收敛阶驱动 |
| `config.py` | $\varepsilon,\xi,\bar k,c,\delta_L,M_{1..5},v_d,a_0,a_d,S,\beta,\gamma,\alpha,\gamma_\eta,C_L$、网格、$\Delta t$ 策略 |
| `run.py` | 2D/3D 驱动（$d$ 参数化），输出能量轨迹与场快照 |

**默认参数起点**（论文参考，非要求）：$\varepsilon=\xi\approx0.17$，$N_g=64$（2D 可 128/256），$M_{1,2,3}\in[4\times10^3,3.2\times10^5]$，$M_{4,5}\approx10^4$，$\bar k=\tfrac12(k_1+k_2)$，$S=|c|,\ c=\tfrac12(k_1-k_2)$（按 §12 标定），$C_L$ 取使 $L+C_L>0$ 的最小裕度。
