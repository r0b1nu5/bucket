function basis_expansion_3rd(X, K, TYPE, NODE)
    # basis_expansion(X, K, TYPE, NODE) generates a multidimensional array of
    # basis expansions evaluated on all points of a multivariate time series.
    #
    # Parameters
    # ------------------
    # X:    Matrix containing N time series of M time points.
    # K:    Maximum order of the basis expansion.
    # TYPE: Type of basis function employed. In this file, we only
    #       implement expansions up to pairwise interactions. Expansions
    #       availables are: polynomial (x_{j}^{k}), polynomial_diff
    #       ((x_{j}-x_{NODE})^{k}), power_series (x_{j}^{k1}*x_{i}^{k2}),
    #       fourier (sin(k*x_{j}) and cos(k*x_{j})), fourier_diff 
    #       (sin(k*(x_{j}-x_{i})) and cos(k*(x_{j}-x_{i}))) and RBF (a model
    #       based on radial basis functions). These functions are shown in 
    #       table I in the main manuscript. 
    # NODE: Unit on which we are performing the reconstruction.
    #
    # Input type
    # ------------------
    # X:    Matrix{Float64}
    # K:    Integer
    # TYPE: String
    # NODE: Integer
    #
    # Output
    # ------------------
    # Expansion: Multidimensional array of size (K+1, M, N) containing the
    # evaluation of all k=0,1,...,K basis functions for all M time points and 
    # all N possible incoming connections. For power_series, (K*K+2) basis
    # functions are employed, and for fourier(_diff), 2*(K+1) are employed. 
    
    N, M = size(X)
    Expansion = zeros(K+1, M, N)
    table = Dict{Int64,Tuple{Int64,Int64,Int64}}()

    if TYPE == "polynomial"
    	for n in 1:N
    		for k in 0:K
                    Expansion[k+1, :, n] .= X[n, :].^k
                end
        end
    elseif TYPE == "polynomial_diff"
    	Xi = similar(X)
        for m in 1:M
        	Xi[:, m] .= X[:, m] .- X[NODE, m]
        end
        for n in 1:N
        	for k in 0:K
                    Expansion[k+1, :, n] .= Xi[n, :].^k
                end
        end    		
    elseif TYPE == "fourier"
    	Expansion = zeros(2*(K+1), M, N)
    	for n in 1:N
                t = 1
                for k in 0:K
                    Expansion[k+t, :, n] .= sin.(k * X[n, :])
                    Expansion[k+t+1, :, n] .= cos.(k * X[n, :])
                    t += 1
                end
        end
    elseif TYPE == "fourier_diff"
    	Expansion = zeros(2*(K+1), M, N)
        Xi = similar(X)
        for m in 1:M
        	Xi[:, m] .= X[:, m] .- X[NODE, m]
        end
        for n in 1:N
        	t = 1
                for k in 0:K
                    Expansion[k+t, :, n] .= sin.(k * Xi[n, :])
                    Expansion[k+t+1, :, n] .= cos.(k * Xi[n, :])
                    t += 1
                end
        end
    elseif TYPE == "power_series"
	Expansion = zeros((K+1)*(K+1)*(K+1), M, Int64(N*(N-1)/2))
	c = 0
        for n1 in 1:N-1
		for n2 in n1+1:N
			c += 1
			table[c] = (NODE,n1,n2)
			for k1 in 1:K+1
				for k2 in 1:K+1
					for k3 in 1:K+1
						for m in 1:M
							Expansion[(K+1)^2*(k1-1) + (K+1)*(k2-1) + k3, m, c] = X[NODE, m]^(k1-1) * X[n1, m]^(k2-1) * X[n2,m]^(k3-1)
						end
					end
				end
			end
		end
	end
    elseif TYPE == "RBF"
    	Expansion = zeros(K, M, N)
        for n in 1:N
                A = vcat(X[n, :], X[NODE, :])
                for m1 in 1:K
                    for m2 in 1:M
                        Expansion[m1, m2, n] = sqrt(2 + norm(A[:, m1] - A[:, m2])^2)
                    end
                end
        end
    else
        error("Invalid TYPE: $TYPE")
    end
    
    return Expansion,table
end

