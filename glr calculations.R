# GLR calculations for toxicity and activity


logit = function(p) log(p/(1-p))

expit = function(u) 1/(1+exp(-u))

# single-dose likelihood
lik.1 = function(p, x, n) (p^x)*((1-p)^(n-x))

# generalized likelihood ratio for comparing
# H1: p <= p0 versus H2: p > p0
# based on single-dose likelihood
glr.1 = function(x, n, p0) {
	p = x/n
	lik.max = lik.1(p,x,n)
	lik.0 = lik.1(p0,x,n)
	lr.obs = lik.max/lik.0
	ifelse(p<p0, lr.obs, 1/lr.obs)
}

# GLR for comparing H1: p <= p0 versus H2: p > p0
# based on single-dose likelihood
# in a different format
glr.i = function(i, xx, nn, p0) glr.1(xx[i], nn[i], p0)

# GLRs for all doses based on single-dose likelihood
glr.all = function(xx, nn, p0) {
	d = length(nn)
	glr = numeric(d)
	for (i in 1:d) glr[i] = glr.i(i, xx, nn, p0)
	p.hat = xx/nn
	list(glr=glr, p.hat=p.hat)
}

# log-likelihood for one binomial observation
# with some logical checks to avoid technicalities
log.lik.1 = function(p, x, n) {
	ll = 0
	if (n>0) {
		if (x==0) ll = n*log(1-p) else {
			if (x==n) ll = n*log(p) else ll = x*log(p)+(n-x)*log(1-p)
		}
	}
	ll
}

# derivative of log.lik.1
der.log.lik.1 = function(p, x, n) {
	der = 0
	if (n>0) {
		if (x==0) der = -n/(1-p) else {
			if (x==n) der = n/p else der = (x/p)-(n-x)/(1-p)
		}
	}
	der
}

# log-likelihood for a vector of binomial observations
# with unrelated success probabilities
# some elements of n may be zero
log.lik = function(pp, xx, nn) {
	d = length(pp)
	ll = 0
	for (i in 1:d) ll = ll + log.lik.1(pp[i],xx[i],nn[i])
	ll
}

# gradient of log.lik
gr.log.lik = function(pp, xx, nn) {
	d = length(pp)
	gr = rep(0, d)
	for (i in 1:d) gr[i] = der.log.lik.1(pp[i],xx[i],nn[i])
	gr
}

# identity matrix
ident = function(d) diag(rep(1,d))

# maximum log-likelihood under an isotonic regression model
# all elements of nn must be > 0
max.log.lik.iso = function(xx, nn, epsilon=0.25) {
	d = length(nn)
	if (d==1) {
		p.hat = xx/nn
		mll = log.lik.1(p.hat, xx, nn)
	} else {
		p.bar = sum(xx)/sum(nn)
		if (p.bar>0&p.bar<1) {
			f = function(q) -log.lik(expit(q),xx,nn)
			g = function(q) {
				p = expit(q)
				-p*(1-p)*gr.log.lik(p,xx,nn)
			}
			e = (1:d)*epsilon
			e = e-mean(e)
			q.0 = rep(logit(p.bar), d) + e # result$getValue(p)
			U = cbind(-ident(d-1),rep(0,d-1)) + cbind(rep(0,d-1),ident(d-1))
			c = rep(0, d-1)
			opt.rst = constrOptim(q.0, f, g, U, c)
			mll = -opt.rst$value
			p.hat = expit(opt.rst$par)
		} else {
			mll = 0
			p.hat = rep(p.bar, d)
		}
	}
	list(mll=mll, p.hat=p.hat)
}

# maximum log-likelihood under an isotonic regression model
# and an additional constraint: p[i] = p0
# all elements of nn must be > 0
max.log.lik.iso.i = function(i, xx, nn, p0, epsilon=0.25) {
	d = length(nn)
	q0 = logit(p0)
	q.ext = function(q) {
		q1 = numeric(d)
		q1[i] = q0
		q1[-i] = q
		q1
	}
	f = function(q) -log.lik(expit(q.ext(q)), xx, nn)
	g = function(q) {
		p = expit(q.ext(q))
		-(p*(1-p)*gr.log.lik(p, xx, nn))[-i]
	}
	q.0 = (q0 + ((1:d)-i)*epsilon)[-i]
	U.0 = cbind(-ident(d-1),rep(0,d-1)) + cbind(rep(0,d-1),ident(d-1))
	U = cbind(U.0[,-i])
	c = rep(0,d-1) - q0*U.0[,i]
	opt.rst = constrOptim(q.0, f, g, U, c)
	list(mll=-opt.rst$value, p.hat=expit(q.ext(opt.rst$par)))
}

# GLR value for H1: p[i]<=p0 versus H2: p[i]>p0
glr.iso.i = function(i, xx, nn, p0) {
	mla = max.log.lik.iso(xx,nn)
	mll = mla$mll
	p.hat = mla$p.hat
	mll.0 = max.log.lik.iso.i(i,xx,nn,p0)$mll
	rat = exp(mll-mll.0)
	glr = ifelse(p.hat[i]<=p0, rat, 1/rat)
	list(glr=glr, mll=mll, p.hat=p.hat)
}

# GLR values for all doses based on joint isotonic likelihood
glr.iso.all = function(xx, nn, p0) {
	mla = max.log.lik.iso(xx,nn)
	mll = mla$mll
	p.hat = mla$p.hat
	d = length(p.hat)
	glr = numeric(d)
	for (i in 1:d) {
		mll.0 = max.log.lik.iso.i(i,xx,nn,p0)$mll
		rat = exp(mll-mll.0)
		glr[i] = ifelse(p.hat[i]<=p0, rat, 1/rat)
	}
	list(glr=glr, mll=mll, p.hat=p.hat)
}

# maximum log-likelihood under a two-sided isotonic constraint with mode set to i in 1:d
# all elements of nn must be > 0
max.log.lik.iso2.i = function(i, xx, nn, epsilon=0.25) {
	p.hat = xx/nn
	mll = log.lik(p.hat, xx, nn)
	d = length(nn)
	inc = p.hat[-1] - p.hat[-d]
	sgn = c(rep(1,i-1), rep(-1,d-i))
	min.sgn.inc = min(sgn*inc)
	if (d>1&!is.na(min.sgn.inc)&min.sgn.inc<0) {
		p.bar = sum(xx)/sum(nn)
		if (p.bar==0|p.bar==1) {
			mll = 0
			p.hat = rep(p.bar, d)
		} else {
			f = function(q) -log.lik(expit(q),xx,nn)
			g = function(q) {
				p = expit(q)
				-p*(1-p)*gr.log.lik(p,xx,nn)
			}
			e = -abs((1:d)-i)*epsilon
			e = e-mean(e)
			q.0 = rep(logit(p.bar), d) + e # result$getValue(p)
			B = cbind(-ident(d-1),rep(0,d-1)) + cbind(rep(0,d-1),ident(d-1))
			U = rbind(c(rep(1,i-1), rep(-1,d-i))*B)
			c = rep(0,d-1)
			opt.rst = constrOptim(q.0, f, g, U, c)
			mll = -opt.rst$value
			p.hat = expit(opt.rst$par)
		}
	}
	list(mll=mll, p.hat=p.hat)
}

# maximum log-likelihood under a two-sided isotonic constranit with unspecified mode
# all elements of nn must be > 0
max.log.lik.iso2 = function(xx, nn) {
	d = length(nn)
	mlls = rep(NA, d); p.hats = matrix(NA, d, d)
	for (i in 1:d) {
		mll.rst = max.log.lik.iso2.i(i,xx,nn)
		mlls[i] = mll.rst$mll; p.hats[i,] = mll.rst$p.hat
	}
	mll = max(mlls)
	j.hat = min(which(mlls==mll))
	glr.d = rep(1, d) # GLRs for de-escalation, comparing H1: OBD < i with H2: OBD >= i
	if (d>1) { for (i in 2:d) glr.d[i] = exp(max(mlls[1:(i-1)])-max(mlls[i:d])) }
	glr.s = rep(1, d) # GLRs for staying, comparing H1: OBD = i with H2: OBD != i
	if (d>1) { for (i in 1:(d-1)) glr.s[i] = exp(mlls[i]-max(mlls[-i])) }
	list(mll=mll, mlls=mlls, j.hat=j.hat, glr.d=glr.d, glr.s=glr.s, p.hats=p.hats)
}

nn = rep(3, 6)
xx = c(0, 1, 2, 3, 3, 2)
max.log.lik.iso2(xx, nn)

# end of program