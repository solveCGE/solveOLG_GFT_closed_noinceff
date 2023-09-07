# formatted reporting
report <- function(reporttext,reportcalc) {
  cursorstart     = 45;
  
  countlinebreaks = 0;
  for (i in 1:nchar(reporttext)) {
    if (substring(reporttext,1,i) != "\n") break;
    countlinebreaks = countlinebreaks + 1;
  }
  
  cursorstart     = max(nchar(reporttext)-countlinebreaks,cursorstart)+3*(reportcalc>=0)+2*(reportcalc<0);
  
  charfill = paste(rep(" ",cursorstart-nchar(reporttext)+countlinebreaks),collapse="");
  
  cat(paste(reporttext,charfill,reportcalc,"\n",sep=""));
  
}

# initializes variables
initvars <- function(timedep,agedep,varnames) {
  
  for (varname in varnames) {
    
    if (timedep && agedep) {
      assign(paste0(varname,"v"), zerosmat(nag,tend), envir = .GlobalEnv);
      assign(paste0(varname,"z"), zerosmat(nag,ncoh), envir = .GlobalEnv);
      assign(paste0(varname,"v0"), zeroscol(nag), envir = .GlobalEnv);
    }
    
    if (!timedep && agedep) {
      assign(varname, zeroscol(nag), envir = .GlobalEnv);
    }
    
    if (timedep && !agedep) {
      assign(varname, zerosrow(tend), envir = .GlobalEnv);
      assign(paste0(varname,"0"), 0, envir = .GlobalEnv);
    }
    
    if (!timedep && !agedep) {
      assign(varname, 0, envir = .GlobalEnv);
    }
  }
}

# fill time-dependent variables with calibration values
fillvars <- function(agedep,varnames) {
  
  for (varname in varnames) {
    
    if (agedep) {
      if (!exists(paste0(varname,"z")) || !exists(paste0(varname,"v")) || !exists(paste0(varname,"v0"))) stop("Variable ", paste0(varname), " is not properly initialized!");
      assign(paste0(varname,"v"), kronecker(eval(parse(text=paste0(varname,"v0"))),onesrow(tend)), envir = .GlobalEnv);
      assign(paste0(varname,"z"), kronecker(eval(parse(text=paste0(varname,"v0"))),onesrow(ncoh)), envir = .GlobalEnv);
    }
    
    if (!agedep) {
      if (!exists(varname) || !exists(paste0(varname,"0"))) stop("Variable ", paste0(varname), " is not properly initialized!");
      assign(varname, eval(parse(text=paste0(varname,"0")))*onesrow(tend), envir = .GlobalEnv);
    }
  }
}

# convert variable from cohort-view to period-view
coh2per <- function(inmat) {
  maxage = nrow(inmat);
  numcoh = ncol(inmat);
  numper = numcoh-(maxage-1);
  
  if (numper<=0) stop("coh2per: insufficient number of columns in input matrix");
  
  outmat = matrix(nrow=maxage,ncol=numper);
  
  for (a in 1:maxage) {
    outmat[a,] = inmat[a,(maxage:numcoh)-(a-1)];
  }
  
  return(outmat);
}

# convert variable from period-view to cohort-view
per2coh <- function(inmat, numrow = NULL, calibvec = NULL) {
  numper = ncol(inmat);
  
  if (!is.null(numrow)) {
    inmat = matrix(inmat,byrow = T, nrow=numrow,ncol=numper);
  }
  if (is.null(calibvec)) {
    calibvec = inmat[,1];
  }
  
  maxage = nrow(inmat);
  numcoh = numper+(maxage-1);
  
  outmat = matrix(nrow=maxage,ncol=numcoh);
  
  for (a in 1:maxage) {
    if (a < maxage) outmat[a,1:(maxage-a)] = rep(calibvec[a],maxage-a);
    outmat[a,(maxage:numcoh)-(a-1)] = inmat[a,];
    if (a > 1) outmat[a,(numcoh-(a-2)):numcoh] = rep(inmat[a,numper],a-1);
  }
  
  return(outmat);
}

# convert variable from cohort-view to period-view and aggregate over age
aggcoh2per <- function(inmat) {
  return(matrix(colSums(coh2per(inmat)),nrow=1));
}

zeroscol <- function(dim1) {
  return(matrix(rep(0,dim1),nrow=dim1,ncol=1));
}

onescol <- function(dim1) {
  return(matrix(rep(1,dim1),nrow=dim1,ncol=1));
}

zerosrow <- function(dim1) {
  return(matrix(rep(0,dim1),ncol=dim1,nrow=1));
}

onesrow <- function(dim1) {
  return(matrix(rep(1,dim1),ncol=dim1,nrow=1));
}

zerosmat <- function(dim1,dim2) {
  return(matrix(rep(0,dim1*dim2),nrow=dim1,ncol=dim2));
}

onesmat <- function(dim1,dim2) {
  return(matrix(rep(1,dim1*dim2),nrow=dim1,ncol=dim2));
}

zeroscube <- function(dim1,dim2,dim3) {
  return(array(rep(0,dim1*dim2*dim3),dim=c(dim1,dim2,dim3)));
}

onescube <- function(dim1,dim2,dim3) {
  return(array(rep(1,dim1*dim2*dim3),dim=c(dim1,dim2,dim3)));
}

format_dec <- function(x, k) trimws(format(round(x, k), nsmall=k));
format_sci <- function(numb) formatC(numb, format = "e", digits = 2);