# simulation study for BOIN-ET

library(boinet)

sim.boinet = function(act.pp, tox.pp, p0, q0, n, M, nr) {
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
	sim.rst = boinet(n.dose=d, start.dose=1, size.cohort=n, n.cohort=M, toxprob=tox.pp, effprob=act.pp, phi=p0,
		delta=q0, tau.T=1, tau.E=1, accrual=1e8, estpt.method="multi.iso", obd.method="max.effprob", n.sim=nr)
	pct.obd = ifelse(tru.obd==0, 0, sim.rst$prop.select[tru.obd])
	N.ave = sum(sim.rst$n.patient)
	prop.trt = sim.rst$n.patient/N.ave
	pct.mtd.over = 100*ifelse(tru.mtd==d, 0, sum(prop.trt[(tru.mtd+1):d]))
	pct.dlt = 100*sum(prop.trt*tox.pp)
	pct.act = 100*sum(prop.trt*act.pp)
	c(pct.obd, pct.mtd.over, pct.dlt, pct.act, sim.rst$prop.select, N.ave, nr)
}

p0 = 0.25 # 0.28 # 
tox.pp.mx = rbind(c(0.01, 0.02, 0.04, 0.10, 0.20, 0.40),
			c(0.05, 0.15, 0.30, 0.45, 0.60, 0.75),
			c(0.15, 0.30, 0.45, 0.60, 0.75, 0.90))
tox.pp =  tox.pp.mx[1,]
d = length(tox.pp)
q0s = c(0.6, 0.3, 0.1)
nm = length(q0s)
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
M = d # 2*d # 
n = 3
nr = 10^3
# nm = 3
sim.smry = matrix(NA, ns*nm, d+6)
file.out = "sim rst 230709b BOIN-ET M=6 hier-tox.csv"

for (s in 1:ns) {
	print(paste(s, date()))
	snr = snrs[s]
	act.pp = act.pp.mx[snr,]
	for (m in 1:nm) sim.smry[(s-1)*nm+m,] = sim.boinet(act.pp, tox.pp, p0, q0s[m], n, M, nr)
}

print(sim.smry)
write.csv(sim.smry, file=file.out)

# end of program