using PyPlot, LinearAlgebra

function compare_barbier(A::Matrix{Float64}, κ::Vector{Float64}, μsS::Float64, σsS::Float64, fig::String="barbier")
	S = size(A)[1]
	Id = diagm(0 => ones(S))

	N0 = inv(Id + μsS*ones(S,S) + σsS*A)*κ
	Ns = (1 - μsS*sum(N0))/(1 - μsS) # Barbier21, Eq. (3)
	N2 = (ones(S,S) - Id)*(N0.^2)

	Δ = -((N0 .- Ns)./N2)*N0'
	E = (1 - μsS)*Δ .+ μsS # Barbier21, Eq. (4)

	C = zeros(S,S,S) # Barbier21, Eq. (5)
	for i in 1:S
		C[i,:,:] = -(N0*N0')./N2[i]
	end

	figure(fig)
	subplot(1,2,1)
	PyPlot.plot(vec(Δ),vec(A),".",color="C0")
	xlabel("Δ")
	ylabel("β")

	subplot(1,2,2)
	for i in 1:S
		for j in 1:S
			if j != i
				PyPlot.plot(C[i,j,[1:j-1;j+1:S]],exp.(A[i,j]).*exp.(A[i,[1:j-1;j+1:S]]),".",color="C1")
			end
		end
	end
	xlabel("corr")
	ylabel("exp(Aij)*exp(Aik)")
end

function compare_barbier(A::Matrix{Float64}, κ::Float64=1., μsS::Float64=5., σsS::Float64=2.7, fig::String="barbier")
	return compare_barbier(A,κ*ones(size(A)[1]),μsS,σsS,fig)
end


function gen_barbier(η::Vector{Float64}, n_iter::Int64=1000)
	S = length(η)
	Id = diagm(0 => ones(S))
	βs = zeros(S,S,n_iter)
	for i in 1:n_iter
		n = zeros(S)
		nmin = -1000.
		while nmin < 1e-15
			A = randn(S,S)
			n = inv(Id + .1*A)*ones(S)
			nmin = minimum(n)
		end
		βs[:,:,i] = (Id + A)*diagm(0 => n./η)
	end

	return βs,η
end

function gen_barbier(S::Int64=21, n_iter::Int64=1000)
	η = rand(S)
	return gen_barbier(η)
end



