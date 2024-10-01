# computes demographic transition (and updates intervivo transfers accordingly)
compdemo = function() {

  # compute demography transition
  for (tt in 2:tend) {
    Nv[1, tt]  <<- NB[tt]
    for (i in 2:nag) {
      Nv[i, tt] <<- Nv[i - 1, tt - 1] * gamv[i - 1, tt - 1]
    }
  }

  Nz   <<- per2coh(Nv)
  N    <<- aggcoh2per(Nz)
  Nc[] <<- colSums(Nv[1:(fag - 1), ])

  # Compute neutral intervivo-transfers by rescaling received transfers
  for (tt in 1:tend) {
    ivgiven          = -sum(Nv[, tt] * ivv[, tt] * (ivv[, tt] < 0))
    ivreceived       = sum(Nv[, tt] * ivv[, tt] * (ivv[, tt] > 0))

    ivv[ivv[, tt] > 0, tt] <<- ivv[ivv[, tt] > 0, tt] * (ivgiven / ivreceived)

    if (abs(sum(ivv[, tt] * Nv[, tt])) > 1e-10) cat("ERROR IN RECOMPDEMO: Unbalanced intervivo transfers!")

  }

  ivz <<- per2coh(ivv)

}

# computes the further life expectancy
lifeexpect = function(gamv) {

  nag = length(gamv)

  lifeexpectageN = zeroscol(nag)
  for (a in (nag - 1):1) {
    lifeexpectageN[a] = (lifeexpectageN[a + 1] + 1) * gamv[a] + gamv[a + 1] * (1 - gamv[a])
  }

  return(lifeexpectageN)
}

# main routine that solves the transition path of the full model
solveOLG_GFT = function(starttime = 1, maxiter = 200L, tol = 1e-5, tune0 = 0.9, tune1 = 0.9) {

  # ===== demography ======#
  compdemo() # recomputes demographic transition

  # ===== final steady state =====#
  FSS()

  actual = gft(PATH, guess, starttime, maxtry = maxiter, tol = tol, tune0 = tune0, tune1 = tune1, iprint = 1L)

  # convert period-view variables to cohort-view variables
  Az       <<- per2coh(Av)
  Consz    <<- per2coh(Consv)
  Hz       <<- per2coh(Hv)
  Lambdaz  <<- per2coh(Lambdav)
  Qz       <<- per2coh(Qv)
  Savz     <<- per2coh(Savv)
  dis_totz <<- per2coh(dis_totv)
  ellz     <<- per2coh(ellv)
  rz       <<- per2coh(rv)
  wz       <<- per2coh(wv)
  pcz      <<- per2coh(pcv)
  yz       <<- per2coh(yv)

  cat(paste0("CHECK SOLUTION:\t\t", max(abs(edy) + abs(edl) + abs(edg) + abs(eda) + abs(ediv) + abs(edab))), "\n")

}

# computes a path of temporary equilibria for given guess of foresigh variables
PATH = function(guess, tstart = 1, tstop = (tend - 1)) {

  V                           <<- t(guess[, 1])
  Hv[fag:nag, ]      <<- t(guess[, (1:(nag - fag + 1)) + 1])
  Lambdav[fag:nag, ] <<- t(guess[, (nag - fag + 3):numforesight])

  tt                        <<- tstart

  xTEPATH0          = onescol(4)

  # TEMPORARY EQUILIBRIA
  while (tt <= tstop) {

    if (tt > 1) {

      xTEPATH0[1] = r[tt - 1] / r0
      xTEPATH0[2] = w[tt - 1] / w0

      # budget-instrument choice
      if (budget_bal == 1) xTEPATH0[3] = 1 + tauW[tt - 1] - tauW0
      if (budget_bal == 2) xTEPATH0[3] = 1 + tauF[tt - 1] - tauF0
      if (budget_bal == 3) xTEPATH0[3] = 1 + tauC[tt - 1] - tauC0
      if (budget_bal == 4) xTEPATH0[3] = 1 + taul[tt - 1] - taul0
      if (budget_bal == 5) xTEPATH0[3] = 1 + tauprof[tt - 1] - tauprof0
      if (budget_bal == 6) xTEPATH0[3] = 1 + cG[tt - 1] - cG0

      xTEPATH0[4] = ab[tt - 1] / ab0

    }

    xout        = multiroot(TE, xTEPATH0, rtol = 1e-12)
    if (abs(xout$estim.precis) > 1e-6) stop(paste("NEWTON METHOD DID NOT CONVERGE IN PERIOD ", tt, "\n", sep = ""))
    xTEPATH0 = xout$root
    tt     <<- tt + 1
    ncalls <<- ncalls + 1
  }

  # REVISED EXPECTATIONS
  actual[, 1]                        = t(V)
  actual[, (1:(nag - fag + 1)) + 1]        = t(Hv[fag:nag, ])
  actual[, (nag - fag + 3):numforesight] = t(Lambdav[fag:nag, ])

  return(actual)
}

# GFT -- GENERALIZED FAIR TAYLOR ALGORITHM (Wilcoxen, P. 1990: "A fast algorithm for solving rational expectations models")
gft = function(f, guess, tstart, maxtry = 200L, tol = 1e-5, tune0 = 0.9, tune1 = 0.9, iprint = 1L) {

  tstart_loop = proc.time()

  trys   <<- 0
  ncalls <<- 0
  nvar   = ncol(guess)
  nyrs   = nrow(guess)
  convergence_reached = F

  cat("\nRunning Generalized Fair Taylor (GFT) algorithm for Transition:\n\n")

  while (trys < maxtry) {

      tstart_iter = proc.time()

      trys <<- trys + 1

      actual = f(guess, tstart, nyrs - 1)

      mis     = max(abs(actual - guess))
      mistol1 = log(mis / tol)
      mistol2 = mis - tol

      tstop_iter = proc.time()

      if (iprint >= 1) cat(paste0("Iteration:  ", formatC(trys, width = 3), "   Time: ",  format_dec(tstop_iter[3] - tstart_iter[3], 3), " sec   log of err/tol: ", format_dec(mistol1, 8), "\n"))

      if (max(mistol2) < 0) {
        convergence_reached = T
        cat(paste0(rep(" ", 51), collapse = ""), "Convergence!\n")
        break
      }

      last_g = guess

      i = nyrs - 1
      while (i >= tstart) {
        tmp1 = actual[i, ] - tune1 * jacob %*% as.matrix(last_g[i + 1, ] - guess[i + 1, ], nvar)
        guess[i, ] = tune0 * tmp1 + (1 - tune0) * last_g[i, ]
        i = i - 1
      }
}

  if (convergence_reached == F) {
    actual = guess
    cat("\nIteration limit exceeded! no convergence!\n")
    cat(paste("Maximum number of iterations  = ", maxtry, "\n", sep = ""))
  }

  lastguess <<- guess

  cat(paste("\nNumber of intraperiod solutions = ", ncalls, "\n", sep = ""))
  cat(paste("Number of vector iterations     = ", trys, "\n", sep = ""))

  tstop_loop = proc.time()

  duration = tstop_loop[3] - tstart_loop[3]
  cat(paste("\nCompution time           = ", duration, " sec\n", sep = ""))

  return(actual)
}

# computation of ISS Jacobians for GFT
JACOB = function(f, guess) {

  tstart = proc.time()

  cat("\nComputing Jacobian...")

  nvar           = ncol(guess)
  jacob    = zerosmat(nvar, nvar)

  actual         = f(guess, 1, 2)
  save_act = actual[1, ]
  step     = 0.001 * guess[1, ]

  for (j in 1:nvar) {
    guess[2, j]   = guess[2, j] + step[j]
    actual           = f(guess, 1, 2)
    guess[2, j]   = guess[2, j] - step[j]
    jacob[, j]    = (actual[1, ] - save_act) / step[j]
  }

  cat("OK!\n")

  tstop = proc.time()

  duration = toString(tstop[3] - tstart[3])
  cat(paste("Compution time           = ", duration, " sec\n", sep = ""))

  tt <<- 1

  return(jacob)
}
