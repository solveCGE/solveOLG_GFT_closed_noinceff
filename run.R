#------------------------------------------------------------------------------#
# "solveOLG_closed_noinceff_R"
# 
# Solves an AK-OLG-Model, closed economy, no income effects in R
# Philip Schuster, July, 2023
#
# run.R: main script, shock definition, then computes transition path
#------------------------------------------------------------------------------#

rm(list = ls());

library(rootSolve);
library(here);

# set working directory
setwd(here()); detach("package:here",unload = T);

cat("\nSIMPLE AUERBACH-KOTLIKOFF CLOSED ECONOMY MODEL (NO INCOME EFFECTS) IN R\n\n");

# load some handy auxiliary functions
source('loadfunctions.R');

# load main functions
source('ALGO.R');
source('SS.R');
source('TE.R');

# control center
tend            = 300;   # number of periods
nag             = 100;   # number of age groups (nag = x => max age = x-1)
budget_bal      = 3;     # budget closing instrument (1.. tauW, 2.. tauF, 3.. tauC, 4.. taul, 5.. tauprof, 6.. cG) to fulfill given debt path

# initializes all variables globally
source('initdata.R');

# run calibration routine
source('calib.R');   # <- change parameters in here

# REPLICATION CHECKS
tt = 1; report("\nREPLICATE SS in tt = 1:",sum(abs(SS(onescol(4)))));
tt = 2; SS(onescol(4));
tt = 1; report("REPLICATE TE in tt = 1:",sum(abs(TE(onescol(4)))));

# COMPUTE JACOBIAN
jacob = JACOB(PATH,guess);
    
# ISS
ISS();

##########################
# POLICY SHOCK SECTION
# ========================
# Note: it is typically easiest to introduce shocks to period-view variables 
#       and then convert them to cohort-view variables using per2coh()

## tax shocks
#tauprof = tauprof*0+0.15;                    # profit tax rate is set to 15%
#tauprof[10:tend] = tauprof[10:tend]*0+0.15;  # profit tax rate is set to 15%, starting in 10 years (announced today)
#tauprof[1:10] = tauprof[1:10]*0+0.15;        # profit tax rate is set to 15% temporarily for next 10 years
#tauWv = tauWv*1.02; tauWz = per2coh(tauWv);  # wage tax is increased by 2%

## delayed retirement (uncomment whole block)
#{rag[1:10] = seq(from=rag0,to=(rag0+2),length.out=10); rag[11:tend]=rag0+2; 
# notretv[] = 0; for (tt in 1:tend) {notretv[1:floor(rag[tt]),tt] = 1; notretv[floor(rag[tt])+1,tt] = rag[tt]-floor(rag[tt]);}
# notretz = per2coh(notretv);}                    # effective retirement age is increased linearly by 2 years over the next 10 years

## pension cut
#pv = pv*0.95; pz = per2coh(pv); # pensions are cut by 5%

## productivity shocks
#thetav = thetav*1.02; thetaz = per2coh(thetav);                       # individual productivity increases permanently by 2%
#thetav[,1:30] = thetav[,1:30]*1.02; thetaz = per2coh(thetav);         # individual productivity increases by 2% in the next 30 years
#TFP = TFP*1.02;                                                       # total factor productivity increases permanently by 2%

## fertility shocks
#NB = NB*1.02;             # 2% more newborns every year
NB[1:30] = NB[1:30]*1.02; # 2% more newborns every year over next 50 years

## mortality shocks
gamv[60:nag,] = 1-(1-gamv[60:nag,])*0.9; gamv[nag,]=0; gamz = per2coh(gamv);  # reduction of old-age mortality by 10%

## shock to the initial capital stock
#K[1] = K0*0.99;                                                       # 1% of capital stock is lost in first period

## change in targeted debt path
#DG[5:20] = seq(from=DG0,to=DG0*0.9,length.out=20-5+1); DG[21:tend] = DG0*0.9; # reduce public debt by 10% over 15 years starting in 5 years

##########################

# Solve transition path to new steady state
solveOLG_GFT(starttime = 1, maxiter = 100, tol = 1e-5);

# some transition plots
plot(0:tend,c(r0,r),type="l",xlab="time",ylab="real interest rate")
plot(0:tend,c(w0,w),type="l",xlab="time",ylab="wage rate")

y_min = min(c(N/N0,Y/Y0,Inv/Inv0,Cons/Cons0,CG/CG0,A/A0,P/P0));
y_max = max(c(N/N0,Y/Y0,Inv/Inv0,Cons/Cons0,CG/CG0,A/A0,P/P0));

plot(0:tend,c(1,N/N0)*100,type="l",ylim = c(y_min,y_max)*100,xlab="time",ylab="period 0 = 100")
lines(0:tend,c(1,Y/Y0)*100,col="blue")
lines(0:tend,c(1,Inv/Inv0)*100,col="red")
lines(0:tend,c(1,Cons/Cons0)*100,col="green")
lines(0:tend,c(1,CG/CG0)*100,col="orange")
lines(0:tend,c(1,A/A0)*100,col="yellow")
lines(0:tend,c(1,P/P0)*100,col="gray")

legend("topright", legend=c("population","GDP","investment","consumption","public consumption","aggregate assets","pension expend."), col=c("black","blue","red","green","orange","yellow","gray"),lty=1);

