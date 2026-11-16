# Sources

The formulas and the numbers the tests check, with where SETH takes them from. "via llm-wiki" means the statement was
read in the [llm-wiki](https://szaghi.github.io/llm-wiki/) bundles `weno-teno` and `high-order-capturing-schemes`,
which cite the paper; the provenance of every formula is in `notes/m2-sources.md`.

- Balsara, D.S., Shu, C.-W. (2000). Monotonicity preserving WENO schemes with increasingly high order of accuracy.
  *J. Comput. Phys.* 160:405. doi:10.1006/jcph.2000.6443 — WENO-7/9 indicators (through ADAM's tables).
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
