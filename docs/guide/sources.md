# Sources

The formulas and the numbers the tests check, with where SETH takes them from. "via llm-wiki" means the statement was
read in the [llm-wiki](https://szaghi.github.io/llm-wiki/) bundles `weno-teno` and `high-order-capturing-schemes`,
which cite the paper; the provenance of every formula is in `notes/m2-sources.md`.

- Balsara, D.S., Shu, C.-W. (2000). Monotonicity preserving WENO schemes with increasingly high order of accuracy.
  *J. Comput. Phys.* 160:405. doi:10.1006/jcph.2000.6443 — the WENO-7/9 indicators, candidates, optimal weights and
  optimal combinations (§III.a–b, Tables I–II) (read in the paper; SETH's tables checked against it and against the
  indicators' definition by `scripts/check_papers.py`).
- Borges, R., Carmona, M., Costa, B., Don, W.S. (2008). An improved WENO scheme for hyperbolic conservation laws.
  *J. Comput. Phys.* 227:3191. doi:10.1016/j.jcp.2007.11.038 — WENO-Z, $\tau_5$, critical-point orders (via llm-wiki).
- Castro, M., Costa, B., Don, W.S. (2011). High order WENO-Z schemes. *J. Comput. Phys.* 230:1766.
  doi:10.1016/j.jcp.2010.11.028 — the rule for $\tau_{2r-1}$ (via Li et al. 2022).
- Don, W.S., Li, P., Wang, K.Y., Wang, Y. (2022). *J. Comput. Phys.* 448:110724 — scale-invariant weights (through
  ADAM).
- Fu, L., Hu, X.Y., Adams, N.A. (2016). A family of high-order targeted ENO schemes. *J. Comput. Phys.* 305:333.
  doi:10.1016/j.jcp.2015.10.037 — TENO, $C_T$, the TENO6 weights; the optimised backgrounds (Eqs. (31)–(33)) from
  which the generator derives the TENO5-opt and TENO6-opt weights of Tables 3–5 (read in the paper).
- Henrick, A.K., Aslam, T.D., Powers, J.M. (2005). Mapped WENO schemes. *J. Comput. Phys.* 207:542.
  doi:10.1016/j.jcp.2005.01.023 — accuracy at critical points (via llm-wiki).
- Jiang, G.-S., Shu, C.-W. (1996). Efficient implementation of WENO schemes. *J. Comput. Phys.* 126:202.
  doi:10.1006/jcph.1996.0130 — WENO-JS, indicators, composite profile.
- Li, Y., Fu, L., Adams, N.A. (2020). A low-dissipation shock-capturing framework with flexible nonlinear dissipation
  control. arXiv:2010.13289 — TENO restated (eqs. 8-13, 18-19, Table 1).
- Li, Q., Huang, X., Yan, P., Duan, Y. (2021). Third-order scale-independent WENO-Z scheme. arXiv:2107.13225 —
  WENO-Z weights, $\tau_3$.
- Li et al. (2022b). *Computers & Fluids*. doi:10.1016/j.compfluid.2022.105547 — Castro's $\tau$ rule, the order of
  $\tau_7$ (via llm-wiki).
- Martín, M.P., Taylor, E.M., Wu, M., Weirs, V.G. (2006). A bandwidth-optimized WENO scheme. *J. Comput. Phys.*
  220:270. doi:10.1016/j.jcp.2006.05.009 — modified wavenumber, resolution efficiency (via llm-wiki).
- NENO6: A six-point neuron-based ENO scheme (2022). arXiv:2207.08500 — TENO6 restated, the NENO6 weights.
- Pirozzoli, S. (2002). *J. Comput. Phys.* doi:10.1006/jcph.2002.7021 — the transfer function (via llm-wiki).
- Rathan, S., Naga Raju, G. (2016). An improved non-linear weights for seventh-order WENO scheme. arXiv:1611.06755 —
  WENO-Z weights, the two $\tau_7$.
- Tang, S. (2024). A high-efficiency adaptive TENO scheme with optimal accuracy order at critical points. *Comput.
  Phys. Commun.* 299:109158. doi:10.1016/j.cpc.2024.109158 — critical-point test functions and orders (Tables 1–3),
  TENO5-AZ: $\tau$ and the adaptive cutoff (Eqs. (19)–(20)) (read in the paper).
- Zhang et al. (2024b). *J. Comput. Phys.* doi:10.1016/j.jcp.2024.112793 — TENO6 orders at critical points (via
  llm-wiki).
- Liang, T., Fu, L. (2024b). A new high-order shock-capturing TENO scheme combined with skew-symmetric-splitting
  method for compressible gas dynamics and turbulence simulation. *Comput. Phys. Commun.* 302:109236 — TENO8's
  background: the face coefficients (Eq. (46)), the linear weights, $\tau_8$ (read in the paper).
- Hong, Z., Ye, Z., Ye, K. (2022). *J. Comput. Phys.* — $C_T = 10^{-7}$ for TENO6 and TENO8 (read in the paper).
- Fu, L., Hu, X.Y., Adams, N.A. (2019). Improved five- and six-point targeted essentially nonoscillatory schemes
  with adaptive dissipation. *AIAA J.* 57(3):1143. doi:10.2514/1.J057370 — TENO5-A, the original, Eqs. (20)–(22):
  the sensor, $C_r = 0.24$, $\xi = 10^{-3}$, $\alpha_1 = 10$, $\alpha_2 = 5$ (read in the paper, the Article in
  Advance; checked by `scripts/check_papers.py`). Its TENO5-A also has optimised linear weights and a sixth-order
  $\tau_5$ (§IV.A, Eq. (25)), and its TENO6-A ($C_r = 0.17$, $\alpha_1 = 10.5$, $\alpha_2 = 4.5$): not in SETH,
  whose TENO5-A is TENO5 with this cutoff (TENO-A as Tang 2024 defines it).
- Huang, H., Liang, T., Fu, L. (2023). A five-point TENO scheme with adaptive dissipation based on a new scale sensor.
  *Commun. Comput. Phys.* 33:1106. doi:10.4208/cicp.oa-2022-0286 — the standard TENO5-A, Eqs. (18)–(20): $C_r$, $\xi$,
  $\alpha_1$, $\alpha_2$ (read in the paper, arXiv:2303.10020; checked by `scripts/check_papers.py`). Its new
  five-point scale sensor is not in SETH.
- Fu, L. (2019). A hybrid method with TENO based discontinuity indicator for hyperbolic conservation laws. *Commun.
  Comput. Phys.* 26:973. doi:10.4208/cicp.oa-2018-0176 — the TENO-based detector, the original, §4.1–4.2, Eqs.
  (4.1)–(4.4): the four three-point stencils, $\beta_{3,3}$, $q = 6$, $C_T = 5\cdot10^{-4}$, $\varepsilon = 10^{-4}$
  in one dimension and $10^{-3}$ in several (read in the paper; SETH's indicators, evaluated from its source, and its
  two values of $\varepsilon$, `SETH_SENSOR_TENO_EPS_1D` and `_MD`, checked by `scripts/check_papers.py`).
- Takagi, S., Fu, L., Wakimura, H., Xiao, F. (2022). A novel high-order low-dissipation TENO-THINC scheme for
  hyperbolic conservation laws. *J. Comput. Phys.* 452:110899. doi:10.1016/j.jcp.2021.110899 — the THINC function
  and its face value (Eqs. (15)–(16), $\beta = 1.6$), the first-order guard (Eq. (17)), $\zeta$ and the detector
  (Eqs. (19)–(20)), Table 1 (read in the paper; `notes/m5-sources.md`; the face formula re-derived by
  `scripts/gen_nonpoly_tables.py`).
- Liang, T., Xiao, F., Shyy, W., Fu, L. (2022). A fifth-order low-dissipation discontinuity-resolving TENO scheme for
  compressible flow simulation. *J. Comput. Phys.* 467:111465. doi:10.1016/j.jcp.2022.111465 — TENO5-DR: the flags
  of Tables 1–3 (reproduced digit for digit by `seth_test_nonpoly`), $\zeta$ (Eq. (17)), the detectors (Eqs. (24),
  (32)), the two-steepness interpolation (Eq. (31), $A = 0.91, 1.23$), the guard (Eq. (33)), the static tests of §3.1
  and Table 4 (read in the paper; `notes/m5-sources.md`).
- Liang, T., Fu, L. (2024). A new type of non-polynomial based TENO scheme for hyperbolic conservation laws. *J.
  Comput. Phys.* 497:112618. doi:10.1016/j.jcp.2023.112618 — stencil targeting (§3.1, Fig. 2), the Gaussian-RBF flux
  rows (Eq. (25), every coefficient re-derived by `scripts/gen_nonpoly_tables.py`, which also finds the labels
  $R'_2$, $R'_3$ swapped against §3.1's stencils), the truncation errors (Eq. (29)), the sensor and $\lambda$ (Eqs.
  (30)–(32)) (read in the paper; `notes/m5-sources.md`).
- Fu, L. (2023). Review of the high-order TENO schemes for compressible gas dynamics and turbulence. *Arch. Comput.
  Methods Eng.* 30:2493. doi:10.1007/s11831-022-09877-7 — the TENO-based detector, §3.4, Eqs. (21)–(23): $q = 6$,
  $C_T = 5\cdot10^{-4}$ (read in the paper; checked by `scripts/check_papers.py`); its $\varepsilon$ is not given there
  (Fu 2019 gives it).
- Fu, L. (2021). Very-high-order TENO schemes with adaptive accuracy order and adaptive dissipation control. *Comput.
  Methods Appl. Mech. Engrg.* — the adaptive cutoff, Eqs. (17)–(19) (read in the paper).
- Li, Y., Chen, C., Ren, Y.-X. (2022). A class of high-order finite difference schemes with minimized dispersion and
  adaptive dissipation for solving compressible flows. *J. Comput. Phys.* 448:110770 — Ren's sensor at critical points
  (§3.2) (read in the paper).
- Wang, B.-S., Don, W.S. (2022). Affine-invariant WENO weights and operator. *Appl. Numer. Math.* 181:630 — the
  affine-invariant descaler (Eqs. (9)–(14), Remark 2) (read in the paper).
- Suresh, A., Huynh, H.T. (1997). Accurate monotonicity-preserving schemes with Runge–Kutta time stepping. *J.
  Comput. Phys.* 136:83 — the monotonicity-preserving bound (Eqs. (2.8)–(2.30)), $\alpha = 4$ (read in the paper).
- Wan, Y., Xia, Y. (2022). *J. Comput. Phys.* — the detector $\eta = \tau_5/(\varepsilon + \bar\beta)$ and its three
  regions (§3, Eqs. (3.14)–(3.15)) (read in the paper).
- De Michele, C., Coppola, G. (2023). Asymptotically entropy-conservative and kinetic-energy preserving numerical fluxes
  for compressible Euler equations. arXiv:2307.10918 — KEEP-PE's internal energy flux (Eq. (15)), in their notation of
  Shima et al. 2021 (read in the paper).
- Shima, N., Kuya, Y., Tamaki, Y., Kawai, S. (2021). Preventing spurious pressure oscillations in split convective form
  discretization for compressible flows. *J. Comput. Phys.* 427:110060. doi:10.1016/j.jcp.2020.110060 — KEEP-PE: the
  KEEP fluxes (Table 1) and the pressure-equilibrium internal energy flux (Eq. (17)) (read in the paper; SETH's kernels,
  evaluated from their source, checked by `scripts/check_papers.py`).
- Pirozzoli, S. (2010). Generalized conservative approximations of split convective derivative operators. *J.
  Comput. Phys.* 229:7180. doi:10.1016/j.jcp.2010.06.006 — the locally conservative form of any order, Eqs. (6),
  (11)–(13) (read in the paper; checked with Shima et al.'s flux).
- Pirozzoli, S. (2011). Stabilized non-dissipative approximations of Euler equations in generalized curvilinear
  coordinates. *J. Comput. Phys.* 230:2997. doi:10.1016/j.jcp.2011.01.001; and (2011) Numerical methods for high-speed
  flows. *Annu. Rev. Fluid Mech.* 43:163. doi:10.1146/annurev-fluid-122109-160718 — read; SETH takes no formula from
  them (curvilinear coordinates; a review).
- Kuya, Y., Totani, K., Kawai, S. (2018). Kinetic energy and entropy preserving schemes for compressible flows by split
  convective forms. *J. Comput. Phys.* 375:823 — the kinetic energy flux of KEP split forms (read in the paper).
- De Vanna, F. et al. (2023a). *Comput. Fluids* — the compression-only Ducros sensor (Eq. (9)) (read in the paper).
- Ducros, F., Ferrand, V., Nicoud, F., Weber, C., Darracq, D., Gacherieu, C., Poinsot, T. (1999). Large-eddy
  simulation of the shock/turbulence interaction. *J. Comput. Phys.* 152:517. doi:10.1006/jcph.1999.6238 — the Ducros
  sensor, Eq. (22), there a factor of Jameson's pressure sensor (Eqs. (23)–(24)) with $\varepsilon = 10^{-30}$ (read in
  the paper; `seth_ducros`, evaluated from its source, checked by `scripts/check_papers.py`; its $\varepsilon$ is the
  caller's).
