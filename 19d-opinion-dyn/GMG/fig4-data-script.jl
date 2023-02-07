using DelimitedFiles, NPZ

p = 3

xxx = npzread("MULTI-PARTY(RM)/EPS_CRITICAL/Ord_all_with_eq_area_0_$p.npz")

nsim,neps,junk = size(xxx["NUM_VOTES"])
nagents = 2001

eff = zeros(p,neps)
n_1st = zeros(p,neps)
n_2nd = zeros(p,neps)

for sim in 1:nsim
	for eps in 1:neps
		votes = xxx["NUM_VOTES"][sim,eps,:]
		m,i1 = findmax(votes)
		m,i2 = findmax([votes[1:i1-1];-1000.;votes[i1+1:p]])
		n_1st[i1,eps] += 1
		n_2nd[i2,eps] += 1
		e = xxx["NUM_AGENTS_INFLUENCED"][sim,eps,1]
		eff[i2,eps] += e
	end
end

for i in 1:p
	writedlm("add_data/fig4a-$i-$p-rm.csv",eff[i,:]'/nsim/nagents,',')
end

writedlm("add_data/n-1st-$p-rm.csv",n_1st',',')
writedlm("add_data/n-2nd-$p-rm.csv",n_2nd',',')
writedlm("add_data/tot-n-agents-influenced-over-$nsim-$p-rm.csv",eff',',')

writedlm("add_data/fig4b-1st-$p-rm.csv",n_1st'./nsim,',')
writedlm("add_data/fig4b-2nd-$p-rm.csv",n_2nd'./nsim,',')




