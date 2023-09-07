TE <- function(inputs) {
  
  resid = zeroscol(4);
  
  r[tt]     <<- r0*inputs[1];
  w[tt]     <<- w0*inputs[2];
  
  if (budget_bal == 1) tauWv[,tt]  <<- tauWv0+(inputs[3]-1);
  if (budget_bal == 2) tauF[tt]    <<- tauF0+(inputs[3]-1);
  if (budget_bal == 3) tauC[tt]    <<- tauC0+(inputs[3]-1);
  if (budget_bal == 4) taul[tt]    <<- taul0+(inputs[3]-1);
  if (budget_bal == 5) tauprof[tt] <<- tauprof0+(inputs[3]-1);
  if (budget_bal == 6) cGV[,tt]    <<- cGV0+(inputs[3]-1);
  
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
  Y[tt]        <<- TFP[tt]*(K[tt]^alpha)*(LS[tt]^(1-alpha));
  LD[tt]       <<- ((1-alpha)*TFP[tt]/((1+tauF[tt])*w[tt]))^(1/alpha)*K[tt];
  
  # INVESTMENT
  uck[tt+1] <<- (r[tt]+delta*(1-tauprof[tt+1]))/(1-tauprof[tt+1]);
  qTob[tt]  <<- (1-tauprof[tt])*alpha*TFP[tt]*(K[tt]/LD[tt])^(alpha-1)+tauprof[tt]*delta+(1-delta);
  Inv[tt]   <<- V[tt+1]/(1+r[tt])-(1-delta)*K[tt];
  K[tt+1]   <<- Inv[tt]+(1-delta)*K[tt];

  # FIRM VALUE
  TaxF[tt]  <<- tauprof[tt]*(Y[tt]-(1+tauF[tt])*w[tt]*LD[tt]-delta*K[tt]);
  Div[tt]   <<- Y[tt]-(1+tauF[tt])*w[tt]*LD[tt]-Inv[tt]-TaxF[tt];
  V[tt]     <<- Div[tt]+V[tt+1]/(1+r[tt]);
  
  # CONSUMPTION
  Lambdav[nag,tt]       <<- pc[tt]^(1-sigma);
  Lambdav[fag:(nag-1),tt] <<- pc[tt]^(1-sigma)+(1/(1+rho)*gamv[fag:(nag-1),tt])^sigma*(1+r[tt])^(sigma-1)*Lambdav[(fag+1):nag,tt+1];
  Omegav[,tt]           <<- Lambdav[,tt]/(pc[tt]^(1-sigma));
  
  Hv[nag,tt]           <<- yv[nag,tt]-dis_totv[nag,tt]*pc[tt]+ivv[nag,tt]+abv[nag,tt];
  Hv[fag:(nag-1),tt]   <<- yv[fag:(nag-1),tt]-dis_totv[fag:(nag-1),tt]*pc[tt]+ivv[fag:(nag-1),tt]+abv[fag:(nag-1),tt]+Hv[(fag+1):nag,tt+1]/(1+r[tt]);
  
  A[tt]          <<- V[tt]+DG[tt];
  Av[,tt]        <<- shareAv[,tt]/Nv[,tt]*A[tt];
  Qv[fag:nag,tt]        <<- (Av[fag:nag,tt]+Hv[fag:nag,tt])/(pc[tt]*Omegav[fag:nag,tt]);
  Consv[fag:nag,tt]     <<- Qv[fag:nag,tt]+dis_totv[fag:nag,tt];
  
  # assets and consumption
  # solve forward in age
  Av[fag,tt+1]         <<- 0;
  Av[(fag+1):nag,tt+1] <<- (1+r[tt])*(Av[fag:(nag-1),tt]+yv[fag:(nag-1),tt]+ivv[fag:(nag-1),tt]+abv[fag:(nag-1),tt]-pc[tt]*Consv[fag:(nag-1),tt]);
  
  Savv[,tt]      <<- Av[,tt]+yv[,tt]+ivv[,tt]+abv[,tt]-pc[tt]*Consv[,tt];

  # AGGREGATION
  A[tt+1]        <<- sum(Av[,tt+1]*Nv[,tt+1]);
  shareAv[,tt+1] <<- (Av[,tt+1]*Nv[,tt+1])/A[tt+1];
  Cons[tt]       <<- sum(Consv[,tt]*Nv[,tt]);
  P[tt]          <<- sum((1-notretv[,tt])*pv[,tt]*Nv[,tt]);                # expend pensions
  tauW[tt]       <<- sum(tauWv[,tt]*notretv[,tt]*ellv[,tt]*thetav[,tt]*Nv[,tt])/LS[tt];
  TaxP           = sum((1-notretv[,tt])*tauWv[,tt]*pv[,tt]*Nv[,tt]);
  Taxl           = sum(taulv[,tt]*Nv[,tt]);

  # GOVERNMENT
  PB[tt]    <<- DG[tt]-DG[tt+1]/(1+r[tt]);
  Exp[tt]   <<- CG[tt]+P[tt];
  Rev[tt]   <<- TaxF[tt]+(tauF[tt]*LD[tt]+tauW[tt]*LS[tt])*w[tt]+tauC[tt]*Cons[tt]+TaxP+Taxl;
  
  # EXCESS DEMANDS
  edy[tt]   <<- Inv[tt]+Cons[tt]+CG[tt]-Y[tt];
  edl[tt]   <<- LD[tt]-LS[tt];
  edg[tt]   <<- Rev[tt]-Exp[tt]-PB[tt];
  eda[tt]   <<- DG[tt]+V[tt]-A[tt];
  eda[tt+1] <<- DG[tt+1]+V[tt+1]-A[tt+1];
  ediv[tt]  <<- -sum(ivv[,tt]*Nv[,tt]);
  edab[tt]  <<- sum((1-gamv[,tt])*Savv[,tt]*Nv[,tt])-ab[tt];
  
  # WALRAS LAW
  edw[tt]  <<- edy[tt]+w[tt]*edl[tt]+edg[tt]+ediv[tt]+edab[tt]+eda[tt]-eda[tt+1]/(1+r[tt]);
  #if (abs(edw[tt])> 1e-8) stop("Error in TE: Walras Law does not hold!\n");
  
  #resid[1] = edy[tt];
  resid[1] = eda[tt+1]; # equivalent to eda[tt+1]
  resid[2] = edl[tt];
  resid[3] = edg[tt];
  resid[4] = edab[tt];
  
  return(resid);

}