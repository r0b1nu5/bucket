using PyPlot

include("volumes.jl")

n = 83
qm = floor(Int64,n/4)

mi = .77
mf = .29
si = .25
sf = .16
li = .20
lf = .13

qis = 0:qm
qfs = -qm:qm

pf_nor = zeros(length(qis),length(qfs))
pf_exp = zeros(length(qis),length(qfs))

figure()

for i in 1:length(qis)
	qi = qis[i]
	for j in 1:length(qfs)
		qf = qfs[j]
		pf_nor[i,j] = vol_qf_normal(qf,qi,n,mf,sf,mi,si)
		pf_exp[i,j] = vol_qf_exp(qf,qi,n,mf,lf,mi,li)
	end
	subplot(2,1,1)
	PyPlot.plot(qfs,pf_nor[i,:],"-o",label="qi = $qi")
	subplot(2,1,2)
	PyPlot.plot(qfs,pf_exp[i,:],"-o",label="qi = $qi")
end

subplot(2,1,1)
xlabel("qf")
ylabel("pf")
legend()
subplot(2,1,2)
xlabel("qf")
ylabel("pf")
legend()



