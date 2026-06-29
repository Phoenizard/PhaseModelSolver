# CLAUDE.md

Guidance for Claude Code working in this repository.

## What this project is

Build a **numerical solver** for the **Wang–Du two-component lipid-membrane phase-field model** (Wang & Du, *J. Math. Biol.* 56 (2008) 347–371). The deliverable is **simulation code**, and reproducing the paper's vesicle morphologies is the end goal.

Two scalar phase fields evolve under an $L^2$ gradient flow of a penalized energy $E_M(\phi,\eta)$:
- `phi` — membrane interface ($\{\phi=0\}$ is the membrane $\Gamma$)
- `eta` — composition order parameter ($\{\eta=0\}$ is the auxiliary surface $\Gamma_\perp$)

The single structural difficulty in this model: a **variable coefficient $k(\eta)$ multiplies the fourth-order principal part** of the bending energy.

**Full model spec (functionals, total energy $E_M$, variable coefficients, domain, BCs, sharp-interface limit, reference parameters) lives in Notion — do not duplicate it here, read it:**
→ Project Plan page: https://app.notion.com/p/38d5084901088169a000fb42310d2d3a
→ Project hub (status, roadmap): https://app.notion.com/p/38d50849010880b9934cd5b61a6c5ad5

## Hard requirements (non-negotiable)

1. **Energy dissipation.** The discrete energy must decrease monotonically. Original-energy dissipation is preferred; modified-energy (SAV-type, with drift) is acceptable. A scheme that does not demonstrably dissipate is not done.
2. **Second order in time.** First-order schemes are out of scope except as a scaffolding step.
3. **Constant-coefficient left-hand operator.** $k(\eta)$ must be handled by **explicit extrapolation to the right-hand side**, so the left operator stays constant-coefficient and FFT/DCT-invertible. Do not build per-step variable-coefficient elliptic solves on the main line.

## Method decisions already made (do not relitigate)

These were refined by a de-biased five-method debate (SAV/MSAV, IEQ/EIEQ, ETD/ETDRK, convex-splitting/stabilized, Lagrange/weighted-SAV) that did **not** assume SAV is the main line. The full transcript is in [doc/debate/](doc/debate/); the deployable scheme is in [doc/note/numerical_scheme.md](doc/note/numerical_scheme.md). Key finding: hard requirement 3 forces *every* method to freeze $k(\eta)\to\bar k$ and extrapolate the remainder, so the variable-coefficient fourth-order part is **not** a method differentiator, and **unconditional dissipation of the true original energy is not freely available to any method** under variable coefficient. The decision therefore is a **layered hybrid**, not a single method:

- **Main line: a layered hybrid.**
  - **Bending fourth-order variable-coefficient principal part** $\varepsilon\Delta(k(\eta)\Delta\phi)$: frozen $\bar k=\tfrac12(k_1+k_2)$ implicit biharmonic **+ same-order $\Delta^2$ stabilization** (convex-splitting style, stabilization constant $S=c=\tfrac12(k_1-k_2)$, **$h$-independent**); remainder explicitly extrapolated to the RHS.
  - **Low-order / penalty terms $A,D,N,P,L$**: **MSAV** (linear auxiliary variable per penalty functional, absorbing the large $M_i$; line tension via $\sqrt{\cdot}$ SAV). $V$ is linear → exact implicit, no auxiliary variable.
  - Sequential decoupled update ($\phi$ first, then $\eta$); both left-hand operators are constant-coefficient and FFT-diagonal; rank-one couplings via Woodbury.
- **Baseline energy is modified/proxy energy, and that is the accepted baseline** (consistent with hard requirement 1: original preferred, modified acceptable). Original-energy is **not** a separate main-line method but a **near-equilibrium Phase-B add-on**.
- **ETD, Lagrange/weighted-SAV, IEQ are conditional add-ons, not the main line.**
  - *ETD* — a stiff-part exact integrator, worth it only when $|c|\ll\bar k$ and the constant-coefficient stiff part is dominant; using it as the *primary* scheme to handle $k(\eta)$ directly remains structurally obstructed (it freezes $\bar k$ like everyone else, and the remainder re-enters explicitly with the same cost).
  - *Lagrange / weighted-SAV* — the Phase-B near-equilibrium original-energy fidelity upgrade (weight $\omega_E$ on the bending term only, with a scalar multiplier; falls back to SAV when the scalar equation has no admissible root). Open risk: solvability of that scalar equation under variable-coefficient fourth-order has no precedent.
  - *IEQ/EIEQ* — largely subsumed by SAV; the only residual niche is the $D$/$L$-coupling pointwise weight $\mathrm{sech}^2(\eta/\xi)W_\phi$, if scalar MSAV cross-closure proves numerically insufficient.
- **Multiple methods may be combined** — the hybrid above *is* such a combination; further mixing is allowed where it genuinely helps.

If a design choice on the main line turns out to block the dissipation proof, the correct move is to **revise the discretization** (auxiliary-variable construction, extrapolation order, stabilization) — design and proof form a loop. Switching the method itself (e.g. bringing in ETD or a Lagrange multiplier, or a hybrid) is allowed when there is a concrete, reasoned basis for it; flag the change and the reason rather than switching silently.

## Development order (current focus: 2D baseline)

The active milestone is the **2D baseline**, in four dependent steps. Each step's output feeds the next:

1. **Algorithm design** → complete algebraic form of the scheme (discretization, energy rewrite, decoupled update).
2. **Dissipation proof** → written proof + conditions (cross-term cancellation, single-field residual, modified-energy monotonicity). Loops with step 1.
3. **Code** → working 2D solver.
4. **Acceptance** → numerical evidence: discrete-energy monotonicity (long-time, drift accumulation, $\Delta t$ robustness, degenerate configs) + confirmed 2nd-order convergence.

**2D and 3D carry equal rigor on dissipation.** The dissipation *proof* is dimension-independent (holds in 2D ⟹ holds verbatim in 3D), but *numerical verification of discrete dissipation* must be done in each dimension separately — they are not interchangeable. 2D is the main arena for systematic dissipation/convergence testing because it is cheap; 3D is for morphology reproduction plus a final dissipation re-check. Do not treat 2D as a throwaway prototype.

## Code conventions

Stack: **Python**, `numpy` / `scipy`, FFT-based spectral spatial discretization on a periodic box.

- Domain: $\Omega=[-\pi,\pi]^d$, periodic BCs. Spatial derivatives in Fourier space via FFT; evaluate nonlinear terms in real space at Fourier nodes (pseudo-spectral).
- Keep the solver **dimension-parameterized** ($d\in\{2,3\}$) so the same scheme runs in 2D and 3D by changing one argument — do not fork separate 2D/3D codebases.
- Separate concerns into modules: energy/variational-derivative definitions, the time-stepping scheme, the spectral operators (FFT inversion of the constant-coefficient operator), and **diagnostics** (energy-trace logging, convergence-order driver). Dissipation diagnostics are first-class, not an afterthought.
- Set $\xi=\varepsilon$. Treat $\varepsilon,\xi$, grid size, $\Delta t$, and all $M_i$ as configurable parameters, not hard-coded constants. Paper reference values are in the Notion Project Plan page — they are starting points, not requirements.
- Adaptive $\Delta t$ chosen to enforce energy decay.

## Definition of done (per scheme)

A scheme is finished only when:
- the discrete energy trace is monotonically non-increasing across long-time runs and robust across $\Delta t$, **and**
- a convergence study confirms second-order temporal accuracy, **and**
- the dissipation proof and its conditions are written down.

"It runs" is not done. "It dissipates and is provably 2nd-order" is done.

## Working norms

- Figures/plots: **label everything in English** regardless of chat language.
- Keep responses concise; answer what is asked.
- Mathematics-first exposition: minimal prose, derive from first principles, don't introduce notation before it's needed, don't simplify for its own sake.
- Confirm mathematical content before writing it into files; read back after writing to confirm changes landed.
- Reuse mature, published analytical approaches before constructing novel proof chains — check the literature framework (Notion / `dissipative_schemes_references`) first.

## Notion as source of truth

The Notion hub and Project Plan page are the authoritative project record. When project direction, model spec, or roadmap is in question, **read Notion rather than relying on this file** — this file is a bootstrap pointer and a rules digest, not the full record. Keep this file in sync if hard requirements or method decisions change.
