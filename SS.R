SS <- function(inputs) {
  
  resid = zeroscol(4);
  
  r[tt]     <<- r0*inputs[1];
  w[tt]     <<- w0*inputs[2];
  
  if (budget_bal == 1) tauWv[,tt]  <<- tauWv0+(inputs[3]-1);
  if (budget_bal == 2) tauF[tt]    <<- tauF0+(inputs[3]-1);
  if (budget_bal == 3) tauC[tt]    <<- tauC0+(inputs[3]-1);
  if (budget_bal == 4) taul[tt]    <<- taul0+(inputs[3]-1);
  if (budget_bal == 5) tauprof[tt] <<- tauprof0+(inputs[3]-1);
  if (budget_bal == 6) cGv[,tt]    <<- cGv0+(inputs[3]-1);
  
  # accidental bequest
  ab[tt]          <<- ab0*inputs[4];
  Nc[tt]          <<- sum(Nv[1:(fag-1),tt]);
  abv[fag:nag,tt] <<- ab[tt]/(N[tt]-Nc[tt])*onescol(nag-fag+1);
  pc[tt]          <<- 1+tauC[tt];
  
  CG[tt]          <<- sum(cGv[,tt]*Nv[,tt]);
  cG[tt]          <<- CG[tt]/N[tt];

  # RETIREMENT
  Nw[tt]          <<- sum(notretv[,tt]*Nv[,tt]); 		     # total nb workers
  Nr[tt] 	        <<- sum((1-notretv[,tt])*Nv[,tt]); 	   # nb retirees
  
  # HOURS SUPPLY
  ellv[fag:nag,tt]     <<- ((w[tt]*(1-tauWv[fag:nag,tt])*thetav[fag:nag,tt]/pc[tt])/parlv0[fag:nag])^sigL;
  dis_totv[fag:nag,tt] <<- (sigL/(1+sigL))*parlv0[fag:nag]*ellv[fag:nag,tt]^((1+sigL)/sigL)-parlv1[fag:nag];
  yv[fag:nag,tt]       <<- notretv[fag:nag,tt]*(w[tt]*(1-tauWv[fag:nag,tt])*ellv[fag:nag,tt]*thetav[fag:nag,tt])+(1-notretv[fag:nag,tt])*(1-tauWv[fag:nag,tt])*pv[fag:nag,tt]-taulv[fag:nag,tt];

  # TOTAL LABOR SUPPLY
  LS[tt]    <<- sum(notretv[,tt]*ellv[,tt]*thetav[,tt]*Nv[,tt]);
  
  # PRODUCTION
  uck[tt]   <<- (r[tt]+delta*(1-tauprof[tt]))/(1-tauprof[tt]);
  qTob[tt]  <<- (1-tauprof[tt])*alpha*Y[tt]/K[tt] + tauprof[tt]*delta + (1-delta);
  K[tt]     <<- LS[tt]*(alpha*TFP[tt]/uck[tt])^(1/(1-alpha));
  Y[tt]     <<- TFP[tt]*(K[tt]^alpha)*(LS[tt]^(1-alpha));
  LD[tt]    <<- Y[tt]*(1-alpha)/((1+tauF[tt])*w[tt]);
  
  # INVESTMENT
  Inv[tt]   <<- delta*K[tt];

  # FIRM VALUE
  TaxF[tt]  <<- tauprof[tt]*(Y[tt]-(1+tauF[tt])*w[tt]*LD[tt]-delta*K[tt]);
  Div[tt]   <<- Y[tt]-(1+tauF[tt])*w[tt]*LD[tt]-Inv[tt]-TaxF[tt];
  V[tt]     <<- (1+r[tt])/r[tt]*Div[tt];
  
  # CONSUMPTION
  Lambdav[nag,tt] <<- pc[tt]^(1-sigma);
  Hv[nag,tt]      <<- yv[nag,tt]-dis_totv[nag,tt]*pc[tt]+ivv[nag,tt]+abv[nag,tt];
  
  for (a in (nag-1):fag) {
    Lambdav[a,tt] <<- pc[tt]^(1-sigma)+((1/(1+rho)*gamv[a,tt])^sigma)*(1+r[tt])^(sigma-1)*Lambdav[a+1,tt];
    Hv[a,tt]      <<- yv[a,tt]-dis_totv[a,tt]*pc[tt]+ivv[a,tt]+abv[a,tt]+Hv[a+1,tt]/(1+r[tt]);
  }
  
  Omegav[,tt]     <<- Lambdav[,tt]/(pc[tt]^(1-sigma));
  
  # assets and consumption
  # solve forward in age
  Av[fag,tt]     <<- 0;
  Qv[fag,tt]     <<- (Av[fag,tt]+Hv[fag,tt])/(pc[tt]*Omegav[fag,tt]);
  Consv[fag,tt]  <<- Qv[fag,tt]+dis_totv[fag,tt];
  
  for (a in (fag+1):nag) {
    Av[a,tt]     <<- (1+r[tt])*(Av[a-1,tt]+yv[a-1,tt]+ivv[a-1,tt]+abv[a-1,tt]-pc[tt]*Consv[a-1,tt]);
    Qv[a,tt]     <<- (Av[a,tt]+Hv[a,tt])/(pc[tt]*Omegav[a,tt]);
    Consv[a,tt]  <<- Qv[a,tt]+dis_totv[a,tt];
  }
  Savv[,tt]      <<- Av[,tt]+yv[,tt]+ivv[,tt]+abv[,tt]-pc[tt]*Consv[,tt];
  
  # AGGREGATION
  A[tt]        <<- sum(Av[,tt]*Nv[,tt]);
  shareAv[,tt] <<- (Av[,tt]*Nv[,tt])/A[tt];
  Cons[tt]     <<- sum(Consv[,tt]*Nv[,tt]);
  P[tt]        <<- sum((1-notretv[,tt])*pv[,tt]*Nv[,tt]);
  tauW[tt]     <<- sum(tauWv[,tt]*notretv[,tt]*ellv[,tt]*thetav[,tt]*Nv[,tt])/LS[tt];
  TaxP         = sum((1-notretv[,tt])*tauWv[,tt]*pv[,tt]*Nv[,tt]);
  Taxl         = sum(taulv[,tt]*Nv[,tt]);

  # GOVERNMENT
  Exp[tt]   <<- CG[tt]+P[tt];
  Rev[tt]   <<- TaxF[tt]+(tauF[tt]*LD[tt]+tauW[tt]*LS[tt])*w[tt]+tauC[tt]*Cons[tt]+TaxP+Taxl;
  PB[tt]    <<- DG[tt]*r[tt]/(1+r[tt]);
  
  # EXCESS DEMANDS
  edy[tt]   <<- Inv[tt]+Cons[tt]+CG[tt]-Y[tt];
  edl[tt]   <<- LD[tt]-LS[tt];
  edg[tt]   <<- Rev[tt]-Exp[tt]-PB[tt];
  eda[tt]   <<- DG[tt]+V[tt]-A[tt];
  ediv[tt]  <<- -sum(ivv[,tt]*Nv[,tt]);
  edab[tt]  <<- sum((1-gamv[,tt])*Savv[,tt]*Nv[,tt])-ab[tt];
  
  # WALRAS LAW
  edw[tt]   <<- 1*edy[tt]+w[tt]*edl[tt]+r[tt]/(1+r[tt])*eda[tt]+edg[tt]+ediv[tt]+edab[tt];
  if (abs(edw[tt])> 1e-8) stop("Error in SS: Walras Law does not hold!\n");
  
  resid[1] = edy[tt];
  resid[2] = edl[tt];
  resid[3] = edg[tt];
  resid[4] = edab[tt];
  
  return(resid);

}

ISS <- function() {
  
  # COMPUTE ISS
  tt <<- 1;
  xout = multiroot(SS,onescol(4),rtol=1e-8);
  if (abs(xout$estim.precis) > 1e-6) {
    stop("NEWTON METHOD DID NOT CONVERGE!\n");
  } else {
    report("\nISS computed:",sum(abs(SS(xout$root))));
  }
}

FSS <- function() {
  
  # COMPUTE FSS
  tt <<- tend;
  xout = multiroot(SS,onescol(4),rtol=1e-8);
  if (abs(xout$estim.precis) > 1e-6) {
    stop("NEWTON METHOD DID NOT CONVERGE!\n");
  } else {
    report("\nISS computed:",sum(abs(SS(xout$root))));
  }
  
  # USE FSS RESULTS FOR GUESS
  guess     <<- kronecker(cbind(V[tend],t(Hv[fag:nag,tend]),t(Lambdav[fag:nag,tend])),onescol(tend));
}
