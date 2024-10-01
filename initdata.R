# number of cohorts
ncoh = tend + (nag - 1)

# parameters and counters
vars_p  = c("alpha",       # capital-share parameter in production function
            "delta",       # depreciation rate
            "sigL",        # elasticity of hours-disutility
            "fag",         # first economically active age group
            "rho",         # discount rate households
            "sigma",       # elasticity of intertemporal substitution
            "tt")          # time index
initvars(timedep = F, agedep = F, vars_p)

# age-dependent variables
vars_a  = c("parlv0",      # multiplicative shift parameter hours-disutillity
            "parlv1")      # additive shift parameter hours-disutillity
initvars(timedep = F, agedep = T, vars_a)

# time-dependent variables
vars_t  =  c("A",          # aggregate assets
             "Cons",       # aggregate consumption
             "CG",         # aggregate public consumption
             "DG",         # public debt
             "Div",        # dividends
             "Exp",        # government expenditure
             "Inv",        # investment
             "K",          # capital stock
             "LD",         # aggregate labor demand
             "LS",         # aggregate labor supply
             "N",          # aggregate population
             "NB",         # number of newborns
             "Nc",         # number of children
             "Nr",         # number of retirees
             "Nw",         # number of workers
             "P",          # pension expenditure
             "PB",         # primary balance
             "Rev",        # government revenue
             "TaxF",       # payroll tax revenue
             "TFP",        # TFP stock
             "V",          # value of representative firm
             "Y",          # output (GDP)
             "ab",         # aggregate accidental bequests
             "cG",         # average per capita public consumption
             "eda",        # excess demand asset market
             "edab",       # (excess demand) accidental bequests
             "edg",        # (excess demand) government budget
             "ediv",       # (excess demand) intervivo transfers
             "edl",        # excess demand labor market
             "edw",        # (excess demand) Walras' Law
             "edy",        # excess demand goods market
             "iv",         # aggregate received intervivo transfers
             "pc",         # price of consumption
             "qTob",       # Tobin's q
             "r",          # real interest rate
             "rag",        # retirement age group
             "tauC",       # consumption tax rate
             "tauF",       # payroll tax rate
             "tauW",       # average wage tax rate
             "taul",       # lump-sum tax rate
             "tauprof",    # corporate tax rate
             "uck",        # user cost of capital
             "w")          # wage rate
initvars(timedep = T, agedep = F, vars_t)

# age-dependent and time-dependent variables (with v and z suffix)
vars_at =  c("A",          # assets
             "Cons",       # consumption
             "H",          # human wealth (discounted future non-capital income)
             "Lambda",     # ...
             "N",          # population mass
             "Omega",      # marginal propensity to consume (out of life-time wealth)
             "Q",          # compound consumption (incl. disutility of labor)
             "Sav",        # end-of-period savings
             "ab",         # accidental bequest
             "dis_tot",    # disutility of labor
             "cG",         # public consumption
             "ell",        # hours worked
             "gam",        # conditional survival probability
             "iv",         # intervivo transfers
             "notret",     # not retired indicator
             "p",          # pension income
             "pc",         # price of consumption
             "r",          # real interest rate
             "shareA",     # share in total assets per age group
             "tauC",       # consumption tax rate
             "tauW",       # wage tax rate
             "taul",       # lump-sum tax rate
             "theta",      # productivity (age-specific)
             "w",          # wage rate
             "y")          # per-period labor and pension income
initvars(timedep = T, agedep = T, vars_at)
