using PyPlot

function check_new(A1::Array{Float64,2}, A1h::Array{Float64,2}, a2::Array{Float64,1}, a2h::Array{Float64,1}, g::Float64, gh::Array{Float64,1}, l::Int64, nu::Float64, T::Float64, Ls::Array{Float64,2}, ls::Tuple{Int64,Int64,Int64}, ks::Tuple{Int64,Int64,Int64})
	n = size(A1)[1]

	lmin,lmax,dl = ls
	kmin,kmax,dk = ks

	ma1 = max(maximum(A1),maximum(A1h))
	ma2 = max(maximum(a2),maximum(a2h))
	mi1 = min(minimum(A1),minimum(A1h))
	mi2 = min(minimum(a2),minimum(a2h))
	
	figure()
	subplot(1,2,1)
	PyPlot([mi1,ma1],[mi1,ma1],"--k")
	xlabel("L_{ij}")
	ylabel("hat{L}_{ij}")
	subplot(1,2,2)
	PyPlot.plot([mi2,ma2],[mi2,ma2],"--k")
	xlabel("d_i")
	ylabel("hat{d}_i")

	for i in 1:n
		subplot(1,2,2)
		PyPlot.plot(a2[i],a2h[i],".",color="C1")
		for j in i:n
			subplot(1,2,1)
			PyPlot.plot(A1[i,j],A1h[i,j],".",color="C0")
		end
	end

	x,y = size(Ls)
	mi = minimum(Ls)
	ma = maximum(Ls)
	
	figure()
	subplot(1,2,1)
	PyPlot.plot([nu,nu],[mi,ma],"--k")
	xlabel("nu")
	ylabel("L")
	subplot(1,2,2)
	PyPlot.plot([l,l],[mi,ma],"--k")
	xlabel("l")
	ylabel("L")
	for i in 1:x
		subplot(1,2,1)
		PyPlot.plot((kmin:dk:kmax)/T,Ls[i,:])
	end
	for i in 1:y
		subplot(1,2,2)
		PyPlot.plot(lmin:dl:lmax,Ls[:,i])
	end

end









end




