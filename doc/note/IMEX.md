# IMEX时间格式对Bending Energy的应用

> 本文从零讲清：弯曲能的梯度流，为什么用 **IMEX（隐显结合）**，以及它一步步怎么得到时间格式、为什么它耗散。只做**时间离散**，算子 $\Delta,\Delta^2$ 保持连续形式。

**一句话理由**：变系数 $k(\eta)$ 是一个**场量（随空间变化），且出现在四阶项中**——四阶项必须隐式处理才能避开 $\Delta t\lesssim h^4$ 的步长限制，但场量 $k(\eta)$ 又让隐式求逆变贵；IMEX 把 $k(\eta)$ 冻成常数 $\bar k$，让四阶项既隐式又能廉价求逆，剩下的偏差显式。

记号：$(f,g):=\int_\Omega fg$ 是 $L^2$ 内积，$\|f\|^2:=(f,f)$。周期边界，分部积分不产生边界项。

---

## 问题与目标

把组分场 $\eta$ 冻结在已知值 $\eta^\star$（子问题的设定），记
$$
a(x):=k(\eta^\star(x))\in[k_2,k_1],\qquad 0<k_2\le k_1.
$$
即 $a(x)$ 是一个**已知正数，随空间变化**的Field。弯曲（主部）能量与其 $L^2$ 梯度流：
$$
E_b[u]=\frac{\varepsilon}{2}\int_\Omega a(x)\,|\Delta u|^2\,dx,\qquad
\frac{\delta E_b}{\delta u}=\varepsilon\,\Delta\!\big(a\,\Delta u\big),\qquad
u_t=-\frac{\delta E_b}{\delta u}.
$$

**连续能量耗散律**（时间格式要复制的目标）：
$$
\frac{d}{dt}E_b=\Big(\frac{\delta E_b}{\delta u},u_t\Big)=-\|u_t\|^2\le0.
$$

**时间格式的目标**：由 $u^n$, $u^{n-1}$ 算出 $u^{n+1}$，使

1. **离散能量不增** $E_b^{n+1}\le E_b^n$，最好对**任意** $\Delta t$ 成立；

2. 求 $u^{n+1}$ **便宜**（不要每步迭代解大系统）。

---

## 入门 CFL 分析

CFL 是**空间离散之后**才有的现象，主观倾向于使用FFT做空间离散化。是显式格式下 $\Delta t$ 超过某个由 $h$ 决定的上限，数值解就发散。在FFT谱分解下，把解拆成波 $e^{i\mathbf k x}$；算子作用在单个波(频率为 $\mathbf k$)上满足：
$$
\Delta\,e^{i\mathbf k x}=-|\mathbf k|^2 e^{i\mathbf k x},\qquad \Delta^2\,e^{i\mathbf k x}=+|\mathbf k|^4 e^{i\mathbf k x}.
$$

对于任意（两层、常系数）线性时间格式，在波 $\mathbf k$ 上都化成

$$
g(\mathbf k)\,\hat u^{n+1}=h(\mathbf k)\,\hat u^n\ \Longrightarrow\ \hat u^{n+1}=\frac{h(\mathbf k)}{g(\mathbf k)}\,\hat u^n .
$$

每步把振幅乘 $G:=h/g$。走 $N$ 步是 $G^N$ 我们希望：$|G|\le1$对任意$\mathbf k$，这样在若干步骤后还是稳定不发散的；相反如果存在某个 $\mathbf k$使得$|G|>1$，数值**发散**。

**重点: 一个项会不会被 CFL 限制，就看它放在 $u^{n+1}$ 侧还是 $u^n$ 侧**

**两个例子**（$L$ 是四阶项，在波 $\mathbf k$ 上乘 $\lambda_{\mathbf k}\sim h^{-4}$，**四阶 ⟹ 此数极大**）：
- **显式**（$L$ 在 $u^n$，进分子）：$g=1,\ h=1-\Delta t\lambda_{\mathbf k}$，$G=1-\Delta t\lambda_{\mathbf k}$。$\lambda$ 大 → $G<-1$ → 发散，要稳须 $\boxed{\Delta t\lesssim h^4}$,可是 $\lambda\sim h^{-4}$非常小，$\Delta t$ 不支持做数值计算。
- **隐式**（$L$ 在 $u^{n+1}$，进分母）：$g=1+\Delta t\lambda_{\mathbf k},\ h=1$，$G=\dfrac{1}{1+\Delta t\lambda_{\mathbf k}}\in(0,1]$ 对任意 $\Delta t$ 都数值安全。

### 两个朴素选择以及对应困难

- **全显式** $\dfrac{u^{n+1}-u^n}{\Delta t}=-\varepsilon\Delta(a\Delta u^n)$：四阶项作用在已知量，$\Delta t$ 不可用。
- **全隐式** $\dfrac{u^{n+1}-u^n}{\Delta t}=-\varepsilon\Delta(a\Delta u^{n+1})$：四阶项作用在未知量 → 无 CFL，但解 $u^{n+1}$ 要对 $\tfrac{I}{\Delta t}+\varepsilon\Delta(a\Delta\cdot)$ 求逆。$a(x)$ 随空间变 ⟹ $\Delta(a\Delta\cdot)$ 作用在一个波上**不再是同一个波**（$a$ 把它打散成许多波）⟹ 不能逐波相除 ⟹ 只能迭代解大系统，**贵**。

**当前问题"场量 $k(\eta)$ + 又在四阶项上"两件叠加**

---

## IMEX 核心想法

要同时拿到"隐式（无 CFL）"和"常系数（廉价求逆）"，把系数 $a(x)$ 拆成常数 + 偏差，让隐式只作用在常数部分。取常数 $\bar k$：
$$
a(x)=\bar k+\big(a(x)-\bar k\big)\ \Longrightarrow\
\varepsilon\Delta(a\Delta u)=\underbrace{\varepsilon\bar k\,\Delta^2 u}_{\text{常系数（隐式）}}+\underbrace{\varepsilon\Delta\big((a-\bar k)\Delta u\big)}_{\text{偏差（显式）}}.
$$
- **常系数主部** $\varepsilon\bar k\Delta^2$ 放**隐式**：四阶项作用在未知量,无 CFL, 且系数常数，没有a的参与,方便求逆。
- **偏差** $\varepsilon\Delta((a-\bar k)\Delta\cdot)$ 放**显式**（用已知 $u^n$）：进右端，不参与求逆。

只要 $\bar k$ 取得够大,如取**中点** $\bar k=\tfrac12(k_1+k_2)$，则
$$
|a(x)-\bar k|\le c:=\tfrac12(k_1-k_2)\quad(\forall x),
$$
偏差幅度 $\le c$，而隐式系数 $\bar k=\tfrac12(k_1+k_2)>c$（因 $k_2>0$）——隐式系数盖过显式偏差。

---

### 一阶 IMEX 格式

把劈分代进后向 Euler：
$$
\boxed{\;\frac{u^{n+1}-u^n}{\Delta t}=-\,\varepsilon\bar k\,\Delta^2 u^{n+1}\;-\;\varepsilon\Delta\big((a-\bar k)\,\Delta u^n\big)\;}
$$
- 左端时间差商；右端第一项隐式（含 $u^{n+1}$），第二项显式（全已知）。
- 整理成"$u^{n+1}$ 在左"：
$$
\big(\tfrac{I}{\Delta t}+\varepsilon\bar k\Delta^2\big)\,u^{n+1}=\tfrac{u^n}{\Delta t}-\varepsilon\Delta\big((a-\bar k)\Delta u^n\big).
$$

可证明一阶IMEX无条件耗散真能量 $E_b$

## 完整弯曲能的二阶 IMEX 格式

前面只取了主部 $(\varepsilon\Delta u)^2$。现把 $G$ 的非线性部分补回——完整弯曲能（$g(u):=u-u^3$）：
$$
E[u]=\int_\Omega\frac{a}{2\varepsilon}\,G^2,\qquad G=\varepsilon\Delta u+\tfrac1\varepsilon g(u).
$$

变分导数（两次分部积分得），按阶分组（$g'(u)=1-3u^2$）：
$$
\frac{\delta E}{\delta u}=\underbrace{\varepsilon\Delta(a\Delta u)}_{\text{四阶主部}}
+\underbrace{\tfrac1\varepsilon\Delta\!\big(a\,g(u)\big)+\tfrac a\varepsilon g'(u)\Delta u+\tfrac a{\varepsilon^3}g'(u)g(u)}_{\text{低阶非线性}} .
$$

四阶部分照冻结 $\bar k$ 隐式 + 同阶 $\Delta^2$ 稳定化；其余全部（变系数偏差 $\varepsilon\Delta((a-\bar k)\Delta u)$ 与三个低阶非线性项）**显式外推到 $u^\star$**。整条弯曲能由此走 IMEX

**二阶（BDF2）格式**（$a=k(\eta^\star)$、$u^\star=2u^n-u^{n-1}$ 均已知）：
$$
\boxed{\;\begin{aligned}
\frac{3u^{n+1}-4u^n+u^{n-1}}{2\Delta t}
=&-\varepsilon\bar k\Delta^2u^{n+1}-\varepsilon S\Delta^2(u^{n+1}-u^\star)\\
&-\varepsilon\Delta\!\big((a-\bar k)\Delta u^\star\big)
-\tfrac1\varepsilon\Delta\!\big(a\,g(u^\star)\big)
-\tfrac a\varepsilon g'(u^\star)\Delta u^\star
-\tfrac a{\varepsilon^3}g'(u^\star)g(u^\star)
\end{aligned}\;}
$$

## BDF2 修正能量耗散证明

### 修正能量何时被接受

原能量 = 弯曲能 $E$ 本身；其单调耗散往往证不出。退而求其次证一个**修正能量** $\mathcal E$，被接受需三条：

- **一致性**：$\mathcal E=E+$ 增量项，且增量项随 $\Delta t\to0$ 消失（故 $\mathcal E$ 是 $E$ 的离散影子，非另一泛函）。
- **下有界**：$\mathcal E\ge0$（否则"在减小"给不出界）。
- **出自格式**：增量项由离散能量恒等式导出，非手工指定。

三条成立则 $0\le E[u^n]\le\mathcal E^n\le\dots\le\mathcal E^0$，即原能量被与 $n$ 无关的常数一致压住、不发散——这是修正能量被接受的合理性所在。

### 构造与定理

Write the increments $\delta^m:=u^m-u^{m-1}$, $w^m:=\Delta\delta^m$ (so $\delta^{n+1}=u^{n+1}-u^n$, $w^n=\Delta(u^n-u^{n-1})$). The modified energy is
$$
\boxed{\;\mathcal E^n:=E[u^n]+\frac{1}{4\Delta t}\|\delta^n\|^2+\frac{\varepsilon S}{2}\|w^n\|^2\;}\qquad(\mathcal E^n\ge0).
$$

**Theorem.** Assume $\|u^m\|_{L^\infty}\le M$ for all $m$, and let $C=C(M,\varepsilon)$ be the constant obtained, via Young / interpolation, from the Lipschitz bound of the lower-order part of $\tfrac{\delta E}{\delta u}$ on $[-M,M]$. If
$$
S\ge\frac{2c^2}{k_2},\qquad \Delta t\le\frac{1}{4\,C(M,\varepsilon)},
$$
then $\mathcal E^{n+1}\le\mathcal E^n$.

**Remark.** Keeping only the fourth-order part (lower-order terms vanish) gives $C=0$, and the condition reduces to $S\ge c^2/k_2$, valid for any $\Delta t$.

### 证明

Start from the BDF2 scheme:
$$
\frac{3u^{n+1}-4u^n+u^{n-1}}{2\Delta t}
=-\varepsilon\bar k\Delta^2u^{n+1}-\varepsilon S\Delta^2(u^{n+1}-u^\star)-\varepsilon\Delta\big((a-\bar k)\Delta u^\star\big)
-\tfrac1\varepsilon\Delta\big(a\,g(u^\star)\big)-\tfrac a\varepsilon g'(u^\star)\Delta u^\star-\tfrac a{\varepsilon^3}g'(u^\star)g(u^\star).
$$
Take the $L^2$ inner product with $\delta^{n+1}$. From integrate by parts, we have $(\Delta^2X,\delta^{n+1})=(\Delta X,w^{n+1})$.

For LHS. With $3u^{n+1}-4u^n+u^{n-1}=2\delta^{n+1}+(\delta^{n+1}-\delta^n)$,
$$
\Big(\frac{3u^{n+1}-4u^n+u^{n-1}}{2\Delta t},\,\delta^{n+1}\Big)
=\frac1{2\Delta t}\Big(\tfrac52\|\delta^{n+1}\|^2-\tfrac12\|\delta^n\|^2+\tfrac12\|\delta^{n+1}-\delta^n\|^2\Big).
$$

For RHS, We have 3 terms: "Main + deviation", "Stabilization", and "Nonlinear" as follows, 

1. For Main + deviation. With $\bar k\Delta u^{n+1}+(a-\bar k)\Delta u^\star=a\Delta u^{n+1}-(a-\bar k)(w^{n+1}-w^n)$,

$$
-\varepsilon\bar k(\Delta u^{n+1},w^{n+1})-\varepsilon\big((a-\bar 	k)\Delta u^\star,w^{n+1}\big)=-\varepsilon(a\Delta u^{n+1},w^{n+1})+\varepsilon\big((a-\bar k)(w^{n+1}-w^n),w^{n+1}\big), 
$$

​	and with $\Delta u^{n+1}=\Delta u^n+w^{n+1}$,
$$
-\varepsilon(a\Delta u^{n+1},w^{n+1})=-\Big[\tfrac\varepsilon2(a\Delta u^{n+1},\Delta u^{n+1})-\tfrac\varepsilon2(a\Delta u^n,\Delta u^n)\Big]-\tfrac\varepsilon2(a\,w^{n+1},w^{n+1}).
$$

2. For Stabilization,

$$
-\varepsilon S(w^{n+1}-w^n,w^{n+1})=-\tfrac{\varepsilon S}2(\|w^{n+1}\|^2-\|w^n\|^2)-\tfrac{\varepsilon S}2\|w^{n+1}-w^n\|^2.
$$

3. Nonlinear. These three terms equal $\tfrac{\delta E}{\delta u}-\varepsilon\Delta(a\Delta u)$, i.e. the variation of $\int(\tfrac a\varepsilon\Delta u\,g+\tfrac a{2\varepsilon^3}g^2)$. Write $[\,\cdot\,]_v:=\tfrac1\varepsilon\Delta(a g(v))+\tfrac a\varepsilon g'(v)\Delta v+\tfrac a{\varepsilon^3}g'(v)g(v)$. Add and subtract at $u^{n+1}$:
	$$
	-([\,\cdot\,]_{u^\star},\delta^{n+1})
	=\underbrace{-([\,\cdot\,]_{u^{n+1}},\delta^{n+1})}_{\text{(i) Taylor}}
	+\underbrace{([\,\cdot\,]_{u^{n+1}}-[\,\cdot\,]_{u^\star},\delta^{n+1})}_{\text{(ii) gap}}.
	$$
	
	(i) $[\,\cdot\,]_v$ is the variation of that energy, so by Taylor
	
	$$
	-([\,\cdot\,]_{u^{n+1}},\delta^{n+1})=-\Big[\textstyle\int(\tfrac a\varepsilon\Delta u\,g+\tfrac a{2\varepsilon^3}g^2)\Big]_{u^n}^{u^{n+1}}+(\text{remainder quadratic in }\delta^{n+1}).
	$$
	(ii) Since $u^{n+1}-u^\star=\delta^{n+1}-\delta^n$, the difference $[\,\cdot\,]_{u^{n+1}}-[\,\cdot\,]_{u^\star}$ carries factors $\Delta(u^{n+1}-u^\star)=w^{n+1}-w^n$ and $u^{n+1}-u^\star=\delta^{n+1}-\delta^n$. Under $\|u\|_\infty\le M$, $g,g'$ are Lipschitz on $[-M,M]$ with bounds depending only on $M$; with these and the interpolation $\|\nabla\delta^{n+1}\|^2\le\eta\|w^{n+1}\|^2+\tfrac1{4\eta}\|\delta^{n+1}\|^2$ Young bounds the remainder of (i) and the gap (ii), giving

	$$
	-([\,\cdot\,]_{u^\star},\delta^{n+1})\le-\Big[\textstyle\int\big(\tfrac a\varepsilon\Delta u\,g+\tfrac a{2\varepsilon^3}g^2\big)\Big]_{u^n}^{u^{n+1}}+\tfrac{\varepsilon k_2}4\|w^{n+1}\|^2+C(\|\delta^{n+1}\|^2+\|\delta^{n+1}-\delta^n\|^2),
	$$
	
	where $C=C(M,\varepsilon)$ is the theorem's constant — the $g,g'$ Lipschitz bounds on $[-M,M]$ after Young/interpolation.

Combin LHS $=$ (1)$+$(2)$+$(3). The quartic difference in (1) and the energy difference in (3) add to $E^{n+1}-E^n$, because $E=\tfrac\varepsilon2(a\Delta u,\Delta u)+\int(\tfrac a\varepsilon\Delta u\,g+\tfrac a{2\varepsilon^3}g^2)$. So
$$
\frac1{2\Delta t}\Big(\tfrac52\|\delta^{n+1}\|^2-\tfrac12\|\delta^n\|^2+\tfrac12\|\delta^{n+1}-\delta^n\|^2\Big)
\le-(E^{n+1}-E^n)-\tfrac\varepsilon2(a\,w^{n+1},w^{n+1})
+\varepsilon\big((a-\bar k)(w^{n+1}-w^n),w^{n+1}\big)
-\tfrac{\varepsilon S}2(\|w^{n+1}\|^2-\|w^n\|^2)-\tfrac{\varepsilon S}2\|w^{n+1}-w^n\|^2
+\tfrac{\varepsilon k_2}4\|w^{n+1}\|^2+C(\|\delta^{n+1}\|^2+\|\delta^{n+1}-\delta^n\|^2).
$$
Two bounds on the right:
$$
\varepsilon\big((a-\bar k)(w^{n+1}-w^n),w^{n+1}\big)\le\tfrac{\varepsilon S}2\|w^{n+1}-w^n\|^2+\tfrac{\varepsilon c^2}{2S}\|w^{n+1}\|^2\quad(\text{Young, cancels the }\|w^{n+1}-w^n\|^2),\qquad
-\tfrac\varepsilon2(a\,w^{n+1},w^{n+1})\le-\tfrac{\varepsilon k_2}2\|w^{n+1}\|^2\ (a\ge k_2).
$$
Move $E^{n+1}-E^n$, $\tfrac1{4\Delta t}(\|\delta^{n+1}\|^2-\|\delta^n\|^2)$ and $\tfrac{\varepsilon S}2(\|w^{n+1}\|^2-\|w^n\|^2)$ to the left (forming $\mathcal E^{n+1}-\mathcal E^n$; the last cancels the stabilization telescoping). The $\delta$-terms collapse to $-\tfrac1{\Delta t}\|\delta^{n+1}\|^2-\tfrac1{4\Delta t}\|\delta^{n+1}-\delta^n\|^2$, giving
$$
\mathcal E^{n+1}-\mathcal E^n\le-\frac1{\Delta t}\|\delta^{n+1}\|^2-\frac1{4\Delta t}\|\delta^{n+1}-\delta^n\|^2-\frac{\varepsilon k_2}2\|w^{n+1}\|^2+\frac{\varepsilon c^2}{2S}\|w^{n+1}\|^2+\frac{\varepsilon k_2}4\|w^{n+1}\|^2+C(\|\delta^{n+1}\|^2+\|\delta^{n+1}-\delta^n\|^2).
$$

Conclude. The coefficient of $\|w^{n+1}\|^2$ is $-\tfrac{\varepsilon k_2}4+\tfrac{\varepsilon c^2}{2S}$; those of $\|\delta^{n+1}\|^2$ and $\|\delta^{n+1}-\delta^n\|^2$ are $-\tfrac1{\Delta t}+C$ and $-\tfrac1{4\Delta t}+C$. Take $S\ge\tfrac{2c^2}{k_2}$ and $\Delta t\le\tfrac1{4C}$; all three are $\le0$, hence $\mathcal E^{n+1}\le\mathcal E^n$. $\square$
