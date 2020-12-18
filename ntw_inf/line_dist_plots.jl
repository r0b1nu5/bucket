using PyPlot, DelimitedFiles

h = .001

for ntw in ["euroroad","polblogs"]
	local t = readdlm("line_dist_data/"*ntw*"_t.csv",',')
	local ts = readdlm("line_dist_data/"*ntw*"_ts.csv",',')
	local d = readdlm("line_dist_data/"*ntw*"_d.csv",',')
	local ds = readdlm("line_dist_data/"*ntw*"_ds.csv",',')
	local p = readdlm("line_dist_data/"*ntw*"_p.csv",',')
	local ps = readdlm("line_dist_data/"*ntw*"_ps.csv",',')
	
	local mit = ts[1,:]
	local met = ts[2,:]
	local mat = ts[3,:]
	local mid = ds[1,:]
	local med = ds[2,:]
	local mad = ds[3,:]
	local mip = ps[1,:]
	local mep = ps[2,:]
	local map = ps[3,:]
	
	local T1 = 1
	local T2 = length(mit)
	
	figure()
	
	subplot(3,1,1)
	PyPlot.fill(h*[T1:T2;T2:-1:T1],[mit[T1:T2];mat[T2:-1:T1]],"k",alpha=.3)
	PyPlot.plot(h*(T1:T2),met[T1:T2],"k")
	for i in 1:length(t[:,1])
		PyPlot.plot(h*(T1:T2),t[i,T1:T2])
	end
	ylabel("x")
	
	subplot(3,1,2)
	PyPlot.fill(h*[T1:T2;T2:-1:T1],[mid[T1:T2];mad[T2:-1:T1]],"k",alpha=.3)
	PyPlot.plot(h*(T1:T2),med[T1:T2],"k")
	for i in 1:length(d[:,1])
		PyPlot.plot(h*(T1:T2),d[i,T1:T2])
	end
	ylabel("x'")

	subplot(3,1,3)
	PyPlot.fill(h*[T1:T2;T2:-1:T1],[mip[T1:T2];map[T2:-1:T1]],"k",alpha=.3)
	PyPlot.plot(h*(T1:T2),mep[T1:T2],"k")
	for i in 1:length(p[:,1])
		PyPlot.plot(h*(T1:T2),p[i,T1:T2])
	end
	ylabel("ψ")
	xlabel("t")
end


	


