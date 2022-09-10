using PyPlot

include("kuramoto.jl")

its = 50

#L = [1. 0. -1.;0. 1. -1.;-1. -1. 2.]
L = [2. -1. -1.;-1. 2. -1.;-1. -1. 2.]
ep = 0.99
w = ep*[1., 1., -2.]

#tfs = Array{Array{Float64,1},1}()

for i in 1:its
	@info "$i"

	th0 = 2pi*rand(3)

	ths,dth,it = kuramoto(L,w,th0,true,false,10000,1e-6,1e-2)
	thf = ths - repeat(ths[end,:]',3,1)

	col = "C0"
	if it < 10000
		test = true
		c = 1
		while test && c <= length(tfs)
			if norm(tfs[c] - thf[:,end]) < 1e-3
				col = "C$c"
				test = false
			end
			c += 1
		end

		if test
			push!(tfs,thf[:,end])
			col = "C$c"
		end
	end

	te = vec(prod(-2pi .< thf[1:end-1,:] .< 4pi,dims=1) .> 0)
	T = maximum(te.*(1:size(thf)[2]))
	
	figure("$(w)'")
	PyPlot.plot(thf[1,1:T],thf[2,1:T],color=col)
end







