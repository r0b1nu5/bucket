using PyPlot, DelimitedFiles

h = .001

for ntw in ["euroroad","usairports"]
	local t = readdlm("line_dist_data/"*ntw*"_t.csv",',')
	local ts = readdlm("line_dist_data/"*ntw*"_ts.csv",',')
	local te = vec(readdlm("line_dist_data/"*ntw*"_te.csv",','))
	local d = readdlm("line_dist_data/"*ntw*"_d.csv",',')
	local ds = readdlm("line_dist_data/"*ntw*"_ds.csv",',')
	local de = vec(readdlm("line_dist_data/"*ntw*"_de.csv",','))
	local p = readdlm("line_dist_data/"*ntw*"_p.csv",',')
	local ps = readdlm("line_dist_data/"*ntw*"_ps.csv",',')
	local pe = vec(readdlm("line_dist_data/"*ntw*"_pe.csv",','))
	local ijs = vec(Int.(readdlm("line_dist_data/"*ntw*"_ijs.csv",',')))
	ids = setdiff(Array(1:length(te)),ijs)

	local mit = ts[1,:]
	local d05t = ts[2,:]
	local d25t = ts[3,:]
	local met = ts[4,:]
	local d75t = ts[5,:]
	local d95t = ts[6,:]
	local mat = ts[7,:]
	local mid = ds[1,:]
	local d05d = ds[2,:]
	local d25d = ds[3,:]
	local med = ds[4,:]
	local d75d = ds[5,:]
	local d95d = ds[6,:]
	local mad = ds[7,:]
	local mip = ps[1,:]
	local d05p = ps[2,:]
	local d25p = ps[3,:]
	local mep = ps[4,:]
	local d75p = ps[5,:]
	local d95p = ps[6,:]
	local map = ps[7,:]
	
	local T1 = 500
	local T2 = length(mit)
	
	figure(ntw*": time series ")

	subplot(3,2,1)
	PyPlot.fill(h*[T1:T2;T2:-1:T1],[mit[T1:T2];mat[T2:-1:T1]],"k",alpha=.2)
	PyPlot.fill(h*[T1:T2;T2:-1:T1],[d05t[T1:T2];d95t[T2:-1:T1]],"k",alpha=.3)
	PyPlot.fill(h*[T1:T2;T2:-1:T1],[d25t[T1:T2];d75t[T2:-1:T1]],"k",alpha=.3)
	PyPlot.plot(h*(T1:T2),met[T1:T2],"k")
	for i in 1:length(t[:,1])
		PyPlot.plot(h*(T1:T2),t[i,T1:T2])
	end
	ylabel("θ")
	subplot(3,2,3)
	PyPlot.fill(h*[T1:T2;T2:-1:T1],[mid[T1:T2];mad[T2:-1:T1]],"k",alpha=.2)
	PyPlot.fill(h*[T1:T2;T2:-1:T1],[d05d[T1:T2];d95d[T2:-1:T1]],"k",alpha=.3)
	PyPlot.fill(h*[T1:T2;T2:-1:T1],[d25d[T1:T2];d75d[T2:-1:T1]],"k",alpha=.3)
	PyPlot.plot(h*(T1:T2),med[T1:T2],"k")
	for i in 1:length(d[:,1])
		PyPlot.plot(h*(T1:T2),d[i,T1:T2])
	end
	ylabel("θ'")
	subplot(3,2,5)
	PyPlot.fill(h*[T1:T2;T2:-1:T1],[mip[T1:T2];map[T2:-1:T1]],"k",alpha=.2)
	PyPlot.fill(h*[T1:T2;T2:-1:T1],[d05p[T1:T2];d95p[T2:-1:T1]],"k",alpha=.3)
	PyPlot.fill(h*[T1:T2;T2:-1:T1],[d25p[T1:T2];d75p[T2:-1:T1]],"k",alpha=.3)
	PyPlot.plot(h*(T1:T2),mep[T1:T2],"k")
	for i in 1:length(p[:,1])
		PyPlot.plot(h*(T1:T2),p[i,T1:T2])
	end
	xlabel("t")
	ylabel("ψ")
	
	#figure(ntw*": snapshot")
		
	subplot(3,2,2)
	PyPlot.plot(ids,te[ids],".",color="C7")
	for i in ijs
		PyPlot.plot(i,te[i],"o",markeredgecolor="k")
	end
	#ylabel("θ")
	subplot(3,2,4)
	PyPlot.plot(ids,de[ids],".",color="C7")
	for i in ijs
		PyPlot.plot(i,de[i],"o",markeredgecolor="k")
	end
	#ylabel("θ'")
	subplot(3,2,6)
	PyPlot.plot(ids,pe[ids],".",color="C7")
	for i in ijs
		PyPlot.plot(i,pe[i],"o",markeredgecolor="k")
	end
	xlabel("id")
end


	


