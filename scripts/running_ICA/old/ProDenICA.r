
G0 <-
function(s,a=1)
### kurtosos negentropy
list(Gs = (s^4)/4, gs = s^3, gps = 3 * s^2)

G1 <-
function(s, a = 1)
###Cosh negentropy function + derivatives
list(Gs = logb(cosh(a * s))/a, gs = tanh(a * s), gps = a * (1 - tanh(a * s)^
	2))

GPois <-
function(x, df = 6, B = 500, order = 1, widen = 1.2, density.return=FALSE, ...)
{
	### x must be a unit vector with mean zero
	### Note; if x has outliers, this will make the gaussian more robust
	### order and widen are robustness parameters in computing the grid
	### This function also computes derivitives
	x <- drop(scale(x))
	n <- length(x)
	rangex <- range(x)
	if(order == 1)
		rx <- rangex
	else {
		rx <- sort(x)[c(order, n - order + 1)]
	}
	rx <- ylim.scale(rx, diff(rx) * widen)
	xg <- seq(from = rx[1], to = rx[2], length = B)
	gaps <- diff(rx)/(B - 1)
	xcuts <- c(min(rangex[1], rx[1]) - gaps/2, xg[ - B] + gaps/2, max(
		rangex[2], rx[2]) + gaps/2)
	ys <- as.vector(table(cut(x, xcuts)))
	gxg <- dnorm(xg)
	bigdata <- list(xg=xg, gxg=gxg, ys=ys)
	##  assign("bigdata",bigdata,frame=0)
	## assign("df", df, frame = 0)
	pois.fit <- gam(ys ~ s(xg, df) + offset(logb(gxg)), family = poisson,
		data = bigdata, ...)
	## Now to get the derivitives
	Gs <- predict(pois.fit)-logb(gxg)
        if(density.return){
        ## package up the function to return for plotting
        density=list(x=xg,y=exp(Gs+logb(gxg)))
      }
	Gs <- Gs + logb(sum(gxg)/sum(fitted(pois.fit)))
	resp <- Gs + residuals(pois.fit, type = "working")
	weights <- pois.fit$weights
	df <- B - pois.fit$df.residual
	pois.refit <- smooth.spline(xg, resp, weights, df)
	Gs <- predict(pois.refit, x, deriv = 0)$y
	gs <- predict(pois.refit, x, deriv = 1)$y
	gps <- predict(pois.refit, x, deriv = 2)$y
	rl=list(Gs = Gs, gs = gs, gps = gps)
        if(density.return)rl$density=density
        rl
}

ICAorthW <-
function(W)
{
	sW <- svd(W)
	sW$u %*% t(sW$v)
}

amari <-
function(V, W, orth = FALSE)
{
### a metric between two orthonormal matrices
	# V and W should be square
	if(orth) A <- abs(t(V) %*% W) else A <- abs(solve(V, W))
	rmax <- apply(A, 1, max)
	cmax <- apply(A, 2, max)
	rsum <- apply(A, 1, sum)
	csum <- apply(A, 2, sum)
	(sum(rsum/rmax - 1) + sum(csum/cmax - 1))/(2 * nrow(A))
}

ProDenICA <-
function(x, k=p, W0=NULL, whiten=FALSE, maxit = 500, thresh = 1e-7, restarts = 0,
	trace = FALSE, Gfunc=G1, eps.rank=1e-7, ...)
{
  this.call=match.call()
  p <- ncol(x)
  n <- nrow(x)
  x <- scale(x, T, F)## x should have mean zero
  if(whiten){## First sphere the data
    sx <- svd(x)
    ##Get the effective rank
    condnum=sx$d;condnum=condnum/condnum[1]
    good=condnum >eps.rank
    rank=sum(good)
    if(k>rank){
      warning(paste("Rank of x is ",rank,"; k reduced from",k," to ",rank,sep=""))
      k=rank
    }
    x <- sqrt(n) * sx$u[,good]# no need to rotate via  %*% t(sx$v[,good])
    whitener=sqrt(n)*scale(sx$v[,good],FALSE,sx$d[good])
  }
  else whitener=NULL
### Get a random start if needed
        if(is.null(W0))	W0 <- matrix(rnorm(p * k), p, k) else   k=ncol(W0)
	W0 <- ICAorthW(W0)
###Initialization
	GS <- matrix(0., n, k)
	gS <- GS
	gpS <- GS
	s <- x %*% W0
	flist <- as.list(1.:k)
	for(j in 1.:k)
		flist[[j]] <- Gfunc(s[, j], ...)
	flist0 <- flist
	crit0 <- mean(sapply(flist0, "[[", "Gs"))
### can try some better starts; only evaluated at first iteration
	while(restarts) {
		W1 <- matrix(rnorm(p * k), p, k)
		W1 <- ICAorthW(W1)
		s <- x %*% W1
		for(j in 1.:k)
			flist[[j]] <- Gfunc(s[, j], ...)
		crit <- mean(sapply(flist, "[[", "Gs"))
		if(trace)
			cat("old crit", crit0, "new crit", crit, "\n")
		if(crit > crit0) {
			crit0 <- crit
			W0 <- W1
			flist0 <- flist
		}
		restarts <- restarts - 1.
	}
###Here is the loop
	nit <- 0
	nw <- 10
	repeat {
		nit <- nit + 1
		gS <- sapply(flist0, "[[", "gs")
		gpS <- sapply(flist0, "[[", "gps")
		t1 <- t(x) %*% gS/n
		t2 <- apply(gpS, 2, mean)
		W1 <- t1 - scale(W0, F, 1/t2)
		W1 <- ICAorthW(W1)
		nw <- amari(W0, W1)
		if(trace)
			cat("Iter", nit, "G", crit0, "crit", nw, "\n")
		W0 <- W1
		if((nit > maxit) | (nw < thresh))
			break
		s <- x %*% W0
		for(j in 1:k)
			flist0[[j]] <- Gfunc(s[, j], ...)
		crit0 <- mean(sapply(flist0, "[[", "Gs"))
	}
  rl=list(W = W0, negentropy = crit0,s= x %*% W0,whitener=whitener,call=this.call)
  rl$density=lapply(flist0,"[[","density")
  class(rl)="ProDenICA"
  rl
}
