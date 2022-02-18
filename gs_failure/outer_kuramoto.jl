using LinearAlgebra, DelimitedFiles

include("../L2B.jl")


function outer_kuramoto(B::Matrix{Float64}, η::Vector{Float64}, y0::Vector{Float64}, save_history::Bool=true, ϵ::Float64=1e-6, max_iter::Int64=10000, h::Float64=.01)
	n,m = size(B)
	Lo = B'*B

	y1 = y0
	y2 = y0
	ys = y0

	err = 1000.
	iter = 0

	while err > ϵ && iter < max_iter
		iter += 1
		if iter%1000 == 0
			c = round(Int64,iter/1000)
			writedlm("temp/ys_$c.csv",ys[:,1:end-1],',')
			ys = ys[:,end]

			@info "iter: $iter"
		end

		k1 = η - Lo*sin.(ys[:,end])
		k2 = η - Lo*sin.(ys[:,end] + h/2*k1)
		k3 = η - Lo*sin.(ys[:,end] + h/2*k2)
		k4 = η - Lo*sin.(ys[:,end] + h*k3)

		dy = (k1 + 2*k2 + 2*k3 + k4)/6

		ys = [ys (ys[:,end]+h*dy)]

		err = maximum(abs.(dy))
	end

	Ys = zeros(m,0)
	for i in 1000:1000:iter
		c = round(Int64,i/1000)
		Ys = [Ys readdlm("temp/ys_$c.csv",',')]
		rm("temp/ys_$c.csv")
	end
	Ys = [Ys ys]
	
	return Ys
end




