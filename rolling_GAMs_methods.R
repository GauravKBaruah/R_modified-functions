###########################################################################################
##functino to calculate the slopes, derivatives and ci's of the data, based on the Gam method from Burthe et al.:
require(mgcv)
is.there = function(x=0, L.CI, U.CI){
  pos<- ifelse(x<U.CI, 1, -1) 
  negs<-ifelse(x>L.CI, -1, 1)
  return(pos+negs)}
##round(length(timeseries)/4)	
gam_smoothing<-function(years, timeseries,knots){
  if(length(which(timeseries<=0))==0){
    gam1<-gam(timeseries~s(as.vector(years), bs="cs", k=knots), family=gaussian(link="log"))}else{
      gam1<-gam(timeseries~s(as.vector(years), bs="cs", k=knots), family=gaussian)}
  time.series.fit<-predict(gam1, newdata=data.frame(years=years), type="response")
  
  X0<-predict(gam1, newdata=data.frame(years=years), type= "lpmatrix")
  X1<-predict(gam1, newdata=data.frame(years=years), type= "lpmatrix")
  Xi<-(X1-X0)
  df <- Xi%*%coef(gam1)              ## ith smooth derivative 
  df.sd <- rowSums(Xi%*%gam1$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  #plot(years,df,type="l",ylim=range(c(df+2*df.sd,df-2*df.sd)))##plot 'em
  #lines(years,df+2*df.sd,lty=2);lines(years,df-2*df.sd,lty=2)
  splines<-data.frame(years=years,deriv=df,U.CI=df+2*df.sd, L.CI=df-2*df.sd)	
  splines$sign<-is.there(0, splines$L.CI, splines$U.CI)/2
  splines$fit<-time.series.fit
  return(splines)}

###################slidingn window
slindingsd<-function(Pop, size,roll) {
  rolwn<-round(length(Pop)*roll/100);  #rolling window should be roughly 50 percent of the data
  lenrol<-length(Pop)-rolwn+1;
  mat<-matrix(NA, nrow=rolwn,ncol=lenrol);
  bsize<-matrix(NA,nrow=rolwn,ncol=lenrol);
  msz<-matrix(NA,nrow=rolwn,ncol=lenrol)
  
  for (i in 1:lenrol){
    nys<-Pop[i:(i+rolwn-1)];
    sys<-size[i:(i+rolwn-1)];
    bsize[,i]<-sys
    mat[,i]<-nys ;
  }
  
  #standard deviation
  sdsize<-numeric()
  Ndd<- numeric()
  
  #autoregressive coefficient at lag 1
  for (i in 1:ncol(mat)){
    sdsize[i]<-sd(bsize[,i])/3
    Ndd[i]<-sd(mat[,i])/max(Pop)
  }
  return(cbind(sdsize,Ndd))}


##############################
rolling.window<-function(phenotype,N)
  
{
  roll.Nd<-numeric()
  roll.ssd<-numeric()
  time<-seq(495,550)
  
 # ss<-loess(phenotype[450:650]~time,span=1)
#  nn<-loess(N[450:650]~time,span =1)
 # sloes<-predict(ss)
  #nloes<-predict(nn)
  sdata<-gam_smoothing(time,phenotype[495:550],-1)
  ndata<-gam_smoothing(time,N[495:550],-1)
  
  for (i in 1:length(ndata[,6])){
    
    roll.Nd[i]<-sd(ndata[1:(1+i),6])/mean(ndata[1:(1+i),6])
    roll.ssd[i]<-sd(sdata[1:(1+i),6])/mean(sdata[1:(1+i),6])
    
    
  }
  #return(cbind(roll.ssd[1:150],roll.Nd[1:150]))
  return(cbind(roll.ssd,roll.Nd,time))
  
}
####

CI_sd<-function(z,N,time){
  
  shift.no<-numeric()
  shift.trait<-numeric()
  z<-na.spline(z)
  N<-na.spline(N)
  end.z<-length(z)
  end.N<-length(N)
  day<-time #[12:length(time)]
  
  for (i in 1:length(z)){
  
  shift.no[i]<-sd(N[1:(i)])/N[1]
  shift.trait[i]<-sd(z[1:(i)])/z[1]

  }
  sdata<-gam_smoothing(day,shift.trait,-1)
  ndata<-gam_smoothing(day,shift.no,-1)
  
  shift.data<-data.frame(Phenotype=sdata[,6],Abundance=ndata[,6],Time=day,Shift.t=shift.trait,shift.n=shift.no)
  return(shift.data)
  
  #return(cbind(shift.no,shift.trait,sdata[,6],ndata[,6]))
  
}



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}

shift <- function(z,N){

shift.no<-numeric()
shift.trait<-numeric()

time<-seq(495,550)

shift.trait<-(z[495]-z[495:550])/z[495]
shift.no<-(N[495]-N[495:550])/N[495]

phenotype<-gam_smoothing(time,shift.trait,-1)
popsize<-gam_smoothing(time,shift.no,-1)

return(cbind(phenotype[,6],popsize[,6],shift.trait,shift.no))
  
}

###########

shift.exp <- function(z,N,time){
  
  shift.no<-numeric()
  shift.trait<-numeric()
  z<-na.spline(z)
  N<-na.spline(N)
  end.z<-length(z)
  end.N<-length(N)
  day<-time[1:length(time)]

  shift.trait<-(z[1]-z[1:end.z])/z[1]
  shift.no<-(N[1]-N[1:end.N])/N[1]
  
  phenotype<-gam_smoothing(day,shift.trait,-1)
  popsize<-gam_smoothing(day,shift.no,-1)
  
  shift.data<-data.frame(Phenotype=phenotype[,6],Abundance=popsize[,6],Time=day,Shift.t=shift.trait,shift.n=shift.no)
  return(shift.data)
  #return(cbind(Phenotype=phenotype[,6],Popsize=popsize[,6],shift.trait,shift.no))
  
}

#####################

genericEWS<-function (timeseries, winsize = 50, detrending = c("no", "gaussian", 
                                                   "loess", "linear", "first-diff"), bandwidth = NULL, span = NULL, 
          degree = NULL, logtransform = FALSE, interpolate = FALSE, 
          AR_n = FALSE, powerspectrum = FALSE) 
{
  timeseries <- data.matrix(timeseries)
  if (dim(timeseries)[2] == 1) {
    Y = timeseries
    timeindex = 1:dim(timeseries)[1]
  }
  else if (dim(timeseries)[2] == 2) {
    Y <- timeseries[, 2]
    timeindex <- timeseries[, 1]
  }
  else {
    warning("not right format of timeseries input")
  }
  if (interpolate) {
    YY <- approx(timeindex, Y, n = length(Y), method = "linear")
    Y <- YY$y
  }
  else {
    Y <- Y
  }
  if (logtransform) {
    Y <- log(Y + 1)
  }
  detrending <- match.arg(detrending)
  if (detrending == "gaussian") {
    if (is.null(bandwidth)) {
      bw <- round(bw.nrd0(timeindex))
    }
    else {
      bw <- round(length(Y) * bandwidth/100)
    }
    smYY <- ksmooth(timeindex, Y, kernel = "normal", bandwidth = bw, 
                    range.x = range(timeindex), x.points = timeindex)
    nsmY <- Y - smYY$y
    smY <- smYY$y
  }
  else if (detrending == "linear") {
    nsmY <- resid(lm(Y ~ timeindex))
    smY <- fitted(lm(Y ~ timeindex))
  }
  else if (detrending == "loess") {
    if (is.null(span)) {
      span <- 25/100
    }
    else {
      span <- span/100
    }
    if (is.null(degree)) {
      degree <- 2
    }
    else {
      degree <- degree
    }
    smYY <- loess(Y ~ timeindex, span = span, degree = degree, 
                  normalize = FALSE, family = "gaussian")
    smY <- predict(smYY, data.frame(x = timeindex), se = FALSE)
    nsmY <- Y - smY
  }
  else if (detrending == "first-diff") {
    nsmY <- diff(Y)
    timeindexdiff <- timeindex[1:(length(timeindex) - 1)]
  }
  else if (detrending == "no") {
    smY <- Y
    nsmY <- Y
  }
  mw <- round(length(Y) * winsize/100)
  omw <- length(nsmY) - mw + 1
  low <- 6
  high <- omw
  nMR <- matrix(data = NA, nrow = mw, ncol = omw)
  x1 <- 1:mw
  for (i in 1:omw) {
    Ytw <- nsmY[i:(i + mw - 1)]
    nMR[, i] <- Ytw
  }
  nARR <- numeric()
  nSD <- numeric()
  nSK <- numeric()
  nKURT <- numeric()
  nACF <- numeric()
  nDENSITYRATIO <- numeric()
  nSPECT <- matrix(0, nrow = omw, ncol = ncol(nMR))
  nCV <- numeric()
  smARall <- numeric()
  smARmaxeig <- numeric()
  detB <- numeric()
  ARn <- numeric()
  nSD <- apply(nMR, 2, sd, na.rm = TRUE)
  for (i in 1:ncol(nMR)) {
    nYR <- ar.ols(nMR[, i], aic = FALSE, order.max = 1, dmean = FALSE, 
                  intercept = FALSE)
    nARR[i] <- nYR$ar
    nSK[i] <- abs(moments::skewness(nMR[, i], na.rm = TRUE))
    nKURT[i] <- moments::kurtosis(nMR[, i], na.rm = TRUE)
    nCV[i] <- nSD[i]/mean(nMR[, i])
    ACF <- acf(nMR[, i], lag.max = 1, type = c("correlation"), 
               plot = FALSE)
    nACF[i] <- ACF$acf[2]
    spectfft <- spec.ar(nMR[, i], n.freq = omw, plot = FALSE, 
                        order = 1)
    nSPECT[, i] <- spectfft$spec
    nDENSITYRATIO[i] <- spectfft$spec[low]/spectfft$spec[high]
    if (AR_n) {
      ARall <- ar.ols(nMR[, i], aic = TRUE, order.max = 6, 
                      demean = F, intercept = F)
      smARall[i] <- ARall$ar[1]
      ARn[i] <- ARall$order
      roots <- Mod(polyroot(c(rev(-ARall$ar), 1)))
      smARmaxeig[i] <- max(roots)
      detB[i] <- (prod(roots))^(2/ARn[i])
    }
  }
  nRETURNRATE = 1/nARR
  timevec <- seq(1, length(nARR))
  KtAR <- cor.test(timevec, nARR, alternative = c("two.sided"), 
                   method = c("kendall"), conf.level = 0.95)
  KtACF <- cor.test(timevec, nACF, alternative = c("two.sided"), 
                    method = c("kendall"), conf.level = 0.95)
  KtSD <- cor.test(timevec, nSD, alternative = c("two.sided"), 
                   method = c("kendall"), conf.level = 0.95)
  KtSK <- cor.test(timevec, nSK, alternative = c("two.sided"), 
                   method = c("kendall"), conf.level = 0.95)
  KtKU <- cor.test(timevec, nKURT, alternative = c("two.sided"), 
                   method = c("kendall"), conf.level = 0.95)
  KtDENSITYRATIO <- cor.test(timevec, nDENSITYRATIO, alternative = c("two.sided"), 
                             method = c("kendall"), conf.level = 0.95)
  KtRETURNRATE <- cor.test(timevec, nRETURNRATE, alternative = c("two.sided"), 
                           method = c("kendall"), conf.level = 0.95)
  KtCV <- cor.test(timevec, nCV, alternative = c("two.sided"), 
                   method = c("kendall"), conf.level = 0.95)
  out <- data.frame(timeindex[mw:length(nsmY)], nARR, nSD, 
                    nSK, nKURT, nCV, nRETURNRATE, nDENSITYRATIO, nACF, round(KtSD$estimate, digits = 3), round(KtAR$estimate,digits = 3),round(KtDENSITYRATIO$estimate,digits=3) )
  colnames(out) <- c("timeindex", "ar1", "sd", "sk", "kurt", 
                     "cv", "returnrate", "densratio", "acf1","KTauSD","KTauAr1","KTauDR")
  return(out)
}


