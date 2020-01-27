using PyPlot, SparseArrays

include("kuramoto.jl")
include("ntw_inf.jl")

figure()


for (ntw,co) in [("uk_w","C0"),("er_w","C1"),("sw_w","C2")]
	@info ntw

	Lsp = readdlm("ntws_data/"*ntw*"_lap_mat_sp.csv",',')
	L = sparse(Lsp[:,1],Lsp[:,2],Lsp[:,3])
	
	n = size(L)[1]

	ls = eigvals(Array(L))
	l2 = ls[2]
	
	a0 = .2
	ws = LinRange(.001,l2,50)
	p0 = 0.
	
	T = 1000
	Ttot = 10000
	h = .1
	
	nJ = Array{Float64,1}()
	
	c = 0 

	for w0 in ws
		c += 1
		@info "c = $c, a0 = $a0, w0 = $w0, p0 = $p0"
		
		Ldh = zeros(n,n)
	
		for i in 1:n
			@info "i = $i"
	
			a = zeros(n)
			a[i] = a0
			w = zeros(n)
			w[i] = w0
			p = zeros(n)
			p[i] = p0
	
			t = min(T,round(Int,pi/(2*w0*h)))
	
			thij = kuramoto_sine(L,zeros(n),zeros(n),a,w,p,T,Ttot)
			for j in i:n
				Ldh[i,j] = ntw_inf_sine(thij[j,t],t,n,a0,w0,p0)
				Ldh[j,i] = Ldh[i,j]
			end
		end
		
		Ldhh = Ldh - Ldh*ones(n)*ones(1,n)
		Lh = pinv(Ldhh)
	
		push!(nJ, norm(Array(L) - Lh)/norm(Array(L)))

		writedlm("data1/"*ntw*"_$(c)_Ldh.csv",Ldhh,',')
	end
	
	PyPlot.plot(ws./l2,nJ,color=co)
	
end



