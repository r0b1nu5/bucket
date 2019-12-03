using LinearAlgebra,SparseArrays

function res_dist(L::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int},Symmetric{Float64,Array{Float64,2}}})
	n = size(L)[1]
	
	ei = eigen(L)
	lsi = Array{Float64,1}()
	for i in 1:n
		if abs(ei.values[i]) < 1e-8
			push!(lsi,0)
		else
			push!(lsi,1/ei.values[i])
		end
	end

	Gai = ei.vectors'*diagm(0 => lsi)*ei.vectors
#	Gai = pinv(Array(L))
	
	Om = zeros(n,n)
	for i in 1:n-1
		for j in i+1:n
			Om[i,j] = Gai[i,i] + Gai[j,j] - 2*Gai[i,j]
			Om[j,i] = Om[i,j]
		end
	end
	
	return Om
end

function res_dist(L::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int},Symmetric{Float64,Array{Float64,2}}},m::Int)
	n = size(L)[1]
	
	Gaim = pinv(Array(L))^m
	
	Omm = zeros(n,n)
	for i in 1:n-1
		for j in i+1:n
			Omm[i,j] = Gaim[i,i] + Gaim[j,j] - 2*Gaim[i,j]
			Omm[j,i] = Omm[i,j]
		end
	end
	
	return Omm
end



