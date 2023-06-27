
# LMS function (adapted from gamlss: "Generalised Additive Models for Location Scale and Shape" package)
# (https://rdrr.io/cran/gamlss/src/R/lms.R) 

#-------------------------------------------------------------------------------
#  the LMS family of distributions
LMS <- c("BCCGo","NO",  "BCPEo", "BCTo") #o stands for original
# the SHASH
theSHASH <-  c("NO", "SHASHo")
#-------------------------------------------------------------------------------
lms_w_GAIC <- function(y, x, n.cyc = 300, 
                       families = LMS,
                       data = NULL, 
                       k=4, #higher penalty for the smoothness
                       cent = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6), 
                       calibration = TRUE,
                       trans.x = FALSE,
                       fix.power = NULL,
                       max.df = NULL,
                       lim.trans = c(0, 1.5),   
                       prof = FALSE,
                       step = 0.1, 
                       legend = FALSE,
                       mu.df = NULL,
                       sigma.df = NULL,
                       nu.df = NULL,
                       tau.df = NULL,
                       c.crit = 0.01,
                       method.pb = "GAIC",
                       ... 
)  
{
  #-------------------------------------------------------------------------------
  # local function 
  findPower <- function(y, x, data = NULL,  lim.trans = c(0, 1.5), prof=FALSE, k=2,  c.crit = 0.01, step=0.1)  
  {
    cat("*** Checking for transformation for x ***", "\n") 
    ptrans<- function(x, p) if (abs(p)<=0.0001) log(x) else I(x^p)
    fn <- function(p) GAIC(gamlss(y~pb(ptrans(x,p)), c.crit = c.crit, trace=FALSE), k=2)
    if (prof) # profile dev
    {
      pp <- seq(lim.trans[1],lim.trans[2], step) 
      pdev <- rep(0, length(pp)) 
      for (i in 1:length(pp)) 
      {
        pdev[i] <- fn(pp[i])  
        #   cat(pp[i], pdev[i], "\n")
      }
      plot(pdev~pp, type="l")
      points(pdev~pp,col="blue")
      par <- pp[which.min(pdev)]
      cat('*** power parameters ', par,"***"," \n") 
    } else
    {
      fn <- function(p) GAIC(gamlss(y~pb(ptrans(x,p)), c.crit = c.crit, trace=FALSE), k=2)
      par <- optimise(fn, lower=lim.trans[1], upper=lim.trans[2])$minimum
      cat('*** power parameters ', par,"***"," \n") 
    }  
    par
  }
  #-------------------------------------------------------------------------------
  ptrans<- function(x, p) if (p==0) log(x) else I(x^p)
  #-------------------------------------------------------------------------------
  # end of local function
  #-------------------------------------------------------------------------------
  ## the families to fit
  FAM <- families      
  ## which method
  method.pb <- match.arg(method.pb)
  ## get the variables  
  ylab <- deparse(substitute(y))
  xlab <- deparse(substitute(x))
  y <- if (!is.null(data)) get(deparse(substitute(y)), envir=as.environment(data)) else y
  x <- if (!is.null(data)) get(deparse(substitute(x)), envir=as.environment(data)) else x
  ## ----------------------------------------------------------------------------- 
  ## if need to check for transformation in x    
  if (is.null(fix.power))
  {
    if (trans.x) # if x^p
    {
      par <- findPower(y, x,   lim.trans = lim.trans, prof=prof, k=2,  c.crit = c.crit, step=0.1)    
      ox <- x
      x <-  ptrans(x,par)
    } 
  } else
  { par <- fix.power
  cat('*** power parameters fixed at ', par,"***"," \n")  
  ox <- x
  x <-  ptrans(x,par)
  }
  ##  starting  model for fitted values for mu 
  ##  Note no sigma is fitted here
  ##  fit the model -------------------------------------------------------------- 
  cat('*** Initial  fit***'," \n")   
  switch(method.pb, 
         "GAIC"= {m0 <- gamlss(y~pb(x, method="GAIC", k=k), sigma.formula=~1, data=data, c.crit = 0.01)}) ## initial fit  finish
  
  ## creating lists etc ----------------------------------------------------------
  failed <- list() 
  fits <- list()
  aic <- AIC(m0, k=2)
  fits <- c(fits, aic) 
  whichdist <- 0
  ## fitting the different models in FAM ------------------------------------------       
  for (i in 1:length(FAM)) 
  {
    cat('*** Fitting', FAM[i], "***","\n")  
    switch(method.pb, 
           "GAIC"= { m1 <- try(gamlss(y ~ pb(x,  method="GAIC", k=k, df=mu.df, max.df = max.df),
                                      sigma.formula = ~pb(x,  method="GAIC", k=k, df=sigma.df, max.df = max.df),
                                      nu.formula = ~pb(x,  method="GAIC", k=k, df=nu.df, max.df = max.df), 
                                      tau.formula = ~pb(x,  method="GAIC", k=k, df=nu.df, max.df = max.df), 
                                      family = FAM[i], data = data, n.cyc = n.cyc,
                                      ...), silent=TRUE)
           })      
    if (any(class(m1)%in%"try-error")) # if fitting failed
    {
      cat(FAM[i], " failed", "\n")
      failed <- c(failed, FAM[i]) 
    }
    else
    {
      aic <- AIC(m1, k=2)
      names(aic) <- FAM[i]
      fits <- c(fits, aic)
      if (AIC(m1, k=2) < AIC(m0, k=2)) 
      {
        m0 <-m1 
        whichdist <- i
      }
    }
  }
  if(whichdist==0) 
  { # refitting the Normal with sigma  if not of the models is any good
    cat('*** Refitting', "NO", "***","\n")  
    m0 <-  switch(method.pb, 
                  "GAIC"= {m0 <- gamlss(y~pb(x, method="GAIC", k=k, max.df = max.df), 
                                        sigma.formula=~pb(x, method="GAIC", k=k, max.df = max.df), 
                                        data  =data, c.crit = 0.01, n.cyc = n.cyc)}) ## initial fit  finish
  }
  ## changing the call t look better in the output -------------------------------
  m0$call$mu.start <- NULL # this works OK
  m0$call$data <- substitute(data) # this is OK
  m0$call$family <- if(whichdist==0) "NO" else FAM[whichdist] # this is OK
  if (is.null(mu.df))    m0$call$formula <- formula(y~pb(x))
  if (is.null(sigma.df)) m0$call$sigma.formula <- formula(~pb(x))
  
  if (is.null(nu.df))    m0$call$nu.formula <- formula(~pb(x))
  if (is.null(tau.df))   m0$call$tau.formula <- formula(~pb(x))
  # convert to string
  stringCall <- toString(deparse(m0$call))
  stringCall <- gsub("x", xlab, stringCall)
  stringCall <- gsub("y", ylab, stringCall)
  # convert bact to call  
  m0$call <- str2lang(stringCall)   
  ## transformation needed -------------------------------------------------------        
  if (trans.x||!is.null(fix.power))   
  { 
    x <- ox
    m0$power <- par 
  }
  ## save the rest information ---------------------------------------------------       
  m0$failed <- failed
  fits <- unlist(fits)
  m0$fits <- fits[order(fits)] 
  m0$xvar <- x#with(DaTa,x)
  m0$y <- y#with(DaTa,y)
  m0$ylab <- ylab
  m0$xlab <- xlab
  if (!is.null(data)) m0$call$data  <- substitute(data)
  ## calibration -----------------------------------------------------------------
  if (calibration)
  {
    calibration(m0, xvar=x, cent=cent, pch = 15, cex = 0.5, col = gray(0.7), ylab=ylab, xlab=xlab, legend=legend)	
  } 
  else 
  {
    centiles(m0, xvar=x, cent=cent, pch = 15, cex = 0.5, 
             col = gray(0.7), ylab=ylab, xlab=xlab, legend=legend)		
  }
  ## saving the fitted functions for mu sigma nu and tau  for prediction --------
  if ("mu"%in%m0$par)
  {
    muFun <- splinefun(x, fitted(m0,"mu"), method="natural")
    m0$muFun <- muFun
  }
  if ("sigma"%in%m0$par)
  {
    sigmaFun <- splinefun(x, fitted(m0,"sigma"), method="natural")
    m0$sigmaFun <- sigmaFun
  }
  if ("nu"%in%m0$par)
  {
    nuFun <- splinefun(x, fitted(m0,"nu"), method="natural")
    m0$nuFun <- nuFun
  }
  if ("tau"%in%m0$par)
  {
    tauFun <- splinefun(x, fitted(m0,"tau"), method="natural")
    m0$tauFun <- tauFun
  }
  #------------------------------------------------------------------------------
  class(m0) <- c("lms", class(m0))
  m0  # save the last model
}