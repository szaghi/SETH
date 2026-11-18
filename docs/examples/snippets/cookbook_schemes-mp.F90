call seth_linear9(u, u, rp, rm)                                   ! any reconstruction kernel
bp = rp ; bm = rm
call seth_mp_bound(u(-2:3), u(-2:3), SETH_MP_ALPHA, bp, bm)       ! the five cells around each upwind cell
