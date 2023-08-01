# simulation study comparing different designs

library(Iso)
source("glr calculations.R")

# new dose based on toxicity and activity data
# toxicity assessment may be based on single-dose or isotonic GLRs
# activity assessment is based on joint likelihood under a two-sided isotonic regression model
# glr.fct = glr.all (single dose) or glr.iso.all (isotonic)
new.dose.glr = function(cur.dose, xx, yy, nn, p0, glr.fct, k1, k2, k3=3.87, k4=NULL, k5=NULL, L=9) {
	max.dose = length(nn)
	max.tried.dose = max(which(nn>0))
	glr.function = ifelse(max.tried.dose==1, glr.all, glr.fct)
	tox.glr = glr.function(xx[1:max.tried.dose], nn[1:max.tried.dose], p0)$glr
	act.mla.rst = max.log.lik.iso2(yy[1:max.tried.dose], nn[1:max.tried.dose])
	act.glr.d = act.mla.rst$glr.d
	act.glr.s = act.mla.rst$glr.s
	if (act.glr.d[cur.dose]>=k4) new.dose = cur.dose - 1 else {
		if (act.glr.s[cur.dose]>=k5) {
			if (tox.glr[cur.dose]<=1/k2) new.dose = cur.dose - 1 else new.dose = cur.dose
		} else {
			new.dose = cur.dose
			if (tox.glr[cur.dose]>=k1) new.dose = cur.dose + 1
			if (tox.glr[cur.dose]<=1/k2) new.dose = cur.dose - 1
		}
	}
	new.dose = min(max.dose, max(1,new.dose))
	elim.doses = which(tox.glr<=(1/k3)) # |act.glr.d>=k6)
	max.avail.dose = ifelse(length(elim.doses)>0, min(elim.doses)-1, max.tried.dose+1)
	new.dose = min(new.dose, max.avail.dose)
	if (new.dose>0&nn[new.dose]>=L) {
	  avail.doses = 1:min(max.avail.dose, max.tried.dose)
	  nad = length(avail.doses)
	  if (nad==0) new.dose = 0 else {
	    avail.p.hat = act.mla.rst$p.hats[act.mla.rst$j.hat,avail.doses]
	    avail.p.total = sum(avail.p.hat)
	    avail.prob = ifelse(avail.p.hat>0, avail.p.hat/avail.p.total, rep(1/nad,nad))
	    new.dose = ifelse(nad==1, avail.doses, sample(avail.doses, 1, prob=avail.prob))
	  }
	}
	new.dose
}

# new dose in the isotonic design of Zang, Lee and Yuan (2014)
new.dose.glr.zly = function(cur.dose, xx, yy, nn, p0, tox.prob.tol=0.95) {
	max.dose = length(nn)
	max.tried.dose = max(which(nn>0))
	j.hat = max.log.lik.iso2(yy[1:max.tried.dose], nn[1:max.tried.dose])$j.hat
	if (j.hat>cur.dose) new.dose = cur.dose + 1
	if (j.hat<cur.dose) new.dose = cur.dose - 1
	if (j.hat==cur.dose) new.dose = ifelse(cur.dose==max.tried.dose, cur.dose+1, cur.dose)
	new.dose = min(max.dose, max(1,new.dose))
	tox.prob = pava(1-pbeta(p0, xx+1, nn-xx+1))
	elim.doses = which(tox.prob>tox.prob.tol)
	max.avail.dose = ifelse(length(elim.doses)>0, min(elim.doses)-1, max.tried.dose+1)
	new.dose = min(new.dose, max.avail.dose)
	new.dose
}

# boundary of BOIN design for given phi and phi.tol
boin.bdry = function(phi, phi.tol=0.4) {
	phi.1 = (1-phi.tol)*phi
	phi.2 = (1+phi.tol)*phi
	lam.1 = log((1-phi.1)/(1-phi))/(log(phi*(1-phi.1)/(phi.1*(1-phi))))
	lam.2 = log((1-phi)/(1-phi.2))/(log(phi.2*(1-phi)/(phi*(1-phi.2))))
	c(lam.1, lam.2)
}

# boundary of TEQR design with given phi and epsilon
teqr.bdry = function(phi, epsilon=0.05) {
	lam.1 = phi-epsilon
	lam.2 = phi+epsilon
	c(lam.1, lam.2)
}

# unit probability mass values based on (p, n)
# under an mTPI design with parameters (phi, epsilon, alpha, beta)
upm = function(p, n, phi, epsilon=0.05, alpha=1, beta=1) {
	a = alpha + n*p
	b = beta + n*(1-p)
	cp.1 = pbeta(phi-epsilon,a,b)
	cp.2 = pbeta(phi+epsilon,a,b)
	upm.lo = cp.1/(phi-epsilon)
	upm.mid = (cp.2-cp.1)/(2*epsilon)
	upm.hi = (1-cp.2)/(1-phi-epsilon)
	c(upm.lo, upm.mid, upm.hi)
}

# bisection method to find x in [lo, hi] that solves f(x) = 0
# for an increasing function f
zero.bisc = function(f, lo, hi, y.tol=1e-6) {
	l = lo; h = hi; x = (l+h)/2; y = f(x)
	while (abs(y)>y.tol) {
		if (y>0) h = x else l = x
		x = (l+h)/2
		y = f(x)
	}
	x
}

# boundary of mTPI design with parameters (n, phi, e, a, b)
mtpi.bdry = function(n, phi, e=0.05, a=1, b=1) {
	f.1 = function(p) {
		upm.val = upm(p,n,phi,epsilon=e,alpha=a,beta=b)
		(upm.val[2]/upm.val[1])-1
	}
	p.1 = zero.bisc(f.1,0,phi+e)
	f.2 = function(p) {
		upm.val = upm(p,n,phi,epsilon=e,alpha=a,beta=b)
		(upm.val[3]/upm.val[2])-1
	}
	p.2 = zero.bisc(f.2,phi-e,1)
	c(p.1, p.2)
}

# boundary of i3+3 design with given parameters (n, phi, epsilon)
i3plus3.bdry = function(n, phi, epsilon=0.05) {
	lam.1 = phi-epsilon
	lam.2 = max(phi+epsilon, lam.1+(1/n))
	c(lam.1, lam.2)
}

# new dose in a generic interval design
new.dose.int = function(cur.dose, xx, nn, p0, tox.prob.tol=0.95, method="BOIN", ...) {
	max.dose = length(nn)
	max.tried.dose = max(which(nn>0))
	if (method=="BOIN") lam = boin.bdry(p0, ...)
	if (method=="TEQR") lam = teqr.bdry(p0, ...)
	if (method=="mTPI") lam = mtpi.bdry(nn[cur.dose], p0, ...)
	if (method=="i3+3") lam = i3plus3.bdry(nn[cur.dose], p0, ...)
	pi.hat = xx[cur.dose]/nn[cur.dose]
	new.dose = cur.dose
	if (pi.hat<lam[1]) new.dose = min(max.dose, cur.dose+1)
	if (pi.hat>lam[2]) new.dose = max(1, cur.dose-1)
	tox.prob = 1-pbeta(p0, xx+1, nn-xx+1)
	elim.doses = which(tox.prob>tox.prob.tol)
	max.avail.dose = ifelse(length(elim.doses)>0, min(elim.doses)-1, max.tried.dose+1)
	new.dose = min(new.dose, max.avail.dose)
	new.dose
}

# dose estimation based on isotonic regression for toxicity
# and two-sided isotonic regression for activity
dose.est.iso2 = function(xx, yy, nn, p0) {
	max.tried.dose = max(which(nn>0))
	tox.p.hat = max.log.lik.iso(xx[1:max.tried.dose], nn[1:max.tried.dose])$p.hat
	safe.doses = which(tox.p.hat<=p0)
	est.mtd = ifelse(length(safe.doses)>0, max(safe.doses), 0)
	act.rst = max.log.lik.iso2(yy[1:max.tried.dose], nn[1:max.tried.dose])
	est.mfad = act.rst$j.hat
	act.q.hat = act.rst$p.hats[est.mfad,]
	if (length(safe.doses)==0) est.obd = 0 else {
		q.hat.safe = act.q.hat[safe.doses]
		max.q.hat.safe = max(q.hat.safe)
		est.obd = min(which(q.hat.safe==max.q.hat.safe))
	}
	list(est.obd=est.obd, est.mtd=est.mtd, est.mfad=est.mfad)
}

# simulate one trial based on GLRs
# with zero correlation between toxicity and activity
# the option "ZLY" refers to the isotonic design of Zang, Lee and Yuan (2014)
# four-letter options are for dose escalation in the AGLL design
sim.one.glr = function(act.pp, tox.pp, p0, n, M, glr.fct=glr.all, k1=1.5, k2=1.05, k3=3.87, k4=NULL, k5=k4, method="GLR") {
	safe.doses = which(tox.pp<=p0)
	tru.mtd = ifelse(length(safe.doses)>0, max(safe.doses), 0)
	max.act = max(act.pp)
	tru.mfad = min(which(act.pp==max.act))
	if (length(safe.doses)==0) tru.obd = 0 else {
		act.pp.safe = act.pp[safe.doses]
		max.act.safe = max(act.pp.safe)
		tru.obd = min(which(act.pp.safe==max.act.safe))
	}
	d = length(act.pp)
	nn = rep(0, d); xx = nn; yy = nn
	m = 1
	new.dose = 1
	while (m<=M&new.dose>0) {
		nn[new.dose] = nn[new.dose] + n
		xx[new.dose] = xx[new.dose] + rbinom(1,n,tox.pp[new.dose])
		yy[new.dose] = yy[new.dose] + rbinom(1,n,act.pp[new.dose])
		cur.dose = new.dose
		if (method=="GLR") new.dose = new.dose.glr(cur.dose, xx, yy, nn, p0, glr.fct, k1, k2, k3=k3, k4=k4, k5=k5)
		if (method=="ZLY") new.dose = new.dose.glr.zly(cur.dose, xx, yy, nn, p0)
		if (nchar(method)>3) new.dose = new.dose.int(cur.dose, xx, nn, p0, method=method)
		m = m + 1
	}
	dose.est = dose.est.iso2(xx, yy, nn, p0)
	est.obd = dose.est$est.obd
	est.mtd = dose.est$est.mtd
	obd.right = (est.obd==tru.obd)
	N = sum(nn); X = sum(xx); Y = sum(yy)
	N.mtd.over = ifelse(tru.mtd<d, sum(nn[(tru.mtd+1):d]), 0)
	obd.rst = rep(0, d)
	obd.rst[est.obd] = 1
	c(obd.right, obd.rst, N.mtd.over, X, Y, N)
}

# simulate many trials based on GLRs
sim.many.glr = function(nr, act.pp, tox.pp, p0, n, M, glr.fct, k1, k2, k3=3.87, k4=k1, k5=k4) {
	sim.rst = matrix(NA, nr, length(act.pp)+5)
	for (r in 1:nr) sim.rst[r,] = sim.one.glr.safe(pp, p0, n, M, glr.fct, k1, k2, k3, k4, k5)
	sim.rst
}

# safe version of a generic function
safe.run = function(fct, ...) {
	rst = try(fct(...))
	if (is.character(rst)) rst = rep(NA, 11)
	rst
}

# summary of simulation results
sim.smry = function(sim.rst) {
	good = !(rowSums(is.na(sim.rst))>0)
	prop.good = mean(good)
	cm = colMeans(sim.rst,na.rm=TRUE)
	100*c(cm[1],cm[8:10]/cm[11],cm[2:7],cm[11]/100,prop.good)
}

p0 = 0.25 # 0.28 # target toxicity rate
tox.pp.mx = rbind(c(0.01, 0.02, 0.04, 0.10, 0.20, 0.40),
			c(0.05, 0.15, 0.30, 0.45, 0.60, 0.75),
			c(0.15, 0.30, 0.45, 0.60, 0.75, 0.90))
tox.pp =  tox.pp.mx[3,] # toxicity profile 
d = length(tox.pp)
act.pp.mx = rbind(c(0.10, 0.20, 0.30, 0.33, 0.35, 0.36),
			c(0.10, 0.20, 0.30, 0.30, 0.30, 0.35),
			c(0.10, 0.20, 0.25, 0.30, 0.30, 0.35),
			c(0.10, 0.20, 0.25, 0.30, 0.30, 0.30),
			c(0.10, 0.20, 0.30, 0.30, 0.30, 0.30),
			c(0.10, 0.20, 0.30, 0.30, 0.28, 0.25),
			c(0.10, 0.20, 0.30, 0.28, 0.25, 0.20),
			c(0.10, 0.20, 0.30, 0.28, 0.30, 0.32),
			c(0.10, 0.20, 0.30, 0.25, 0.35, 0.35))
snrs = 1:nrow(act.pp.mx) # activity scenarios
ns = length(snrs)
M = d # 2*d # total number of cohorts
n = 3 # cohort size
nr = 10^3 # number of simulation replicates
nm = 1 # number of methods to compare
all.smry = NULL # all simulation summaries
file.out = "sim rst 230724b hier-tox M=6.csv"

for (s in 1:ns) {
	snr = snrs[s]
	act.pp = act.pp.mx[snr,]
	sim.rst = array(NA, dim=c(nm,nr,d+5))
	smry = matrix(NA, nm, d+6)
	for (r in 1:nr) {
		if (r%%100==0) print(paste(s, r, date()))
		sim.rst[1,r,] = safe.run(sim.one.glr, act.pp, tox.pp, p0, n, M, method="ZLY")
		sim.rst[2,r,] = safe.run(sim.one.glr, act.pp, tox.pp, p0, n, M, method="TEQR")
		sim.rst[3,r,] = safe.run(sim.one.glr, act.pp, tox.pp, p0, n, M, method="mTPI")
		sim.rst[4,r,] = safe.run(sim.one.glr, act.pp, tox.pp, p0, n, M, k4=Inf)
		sim.rst[5,r,] = safe.run(sim.one.glr, act.pp, tox.pp, p0, n, M, k4=3)
		sim.rst[6,r,] = safe.run(sim.one.glr, act.pp, tox.pp, p0, n, M, k4=2)
		sim.rst[7,r,] = safe.run(sim.one.glr, act.pp, tox.pp, p0, n, M, k4=1.5)
		sim.rst[8,r,] = safe.run(sim.one.glr, act.pp, tox.pp, p0, n, M, k4=1.05)
	}
	for (m in 1:nm) smry[m,] = sim.smry(sim.rst[m,,])
	print(smry)
	all.smry = rbind(all.smry, smry)
}

write.csv(all.smry, file=file.out)

# end of program