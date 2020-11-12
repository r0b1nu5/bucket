using SparseArrays

function L2B(L::Array{Float64,2})
	n = size(L)[1]
	
	B = Array{Float64,2}(undef,n,0)
	w = Array{Float64,1}(undef,0)
	for i in 1:n-1
		for j in i+1:n
			if L[i,j] != 0.0
				ed = zeros(n)
				ed[i] = 1.0
				ed[j] = -1.0
				B = [B ed]
				push!(w,-L[i,j])
			end
		end
	end
	
	return B,w
end

function L2B(L::SparseMatrixCSC{Float64,Int})
	n = size(L)[1]
	
	B = Array{Float64,2}(undef,n,0)
	w = Array{Float64,1}(undef,0)
	Bt = Array{Float64,2}(undef,0,n)
	for i in 1:n-1
		for j in i+1:n
			if L[i,j] != 0.0
				ed = zeros(n)
				edt = zeros(1,n)
				ed[i] = 1.0
				edt[1,i] = 1.0
				ed[j] = -1.0
				edt[1,j] = -1.0
				B = [B ed]
				Bt = [Bt;edt]
				push!(w,-L[i,j])
			end
		end
	end
	
	return sparse(B),w,sparse(Bt)
end

function L2B_dir(L::Array{Float64,2})
	n = size(L)[1]

	S = Array{Float64,2}(undef,n,0) # "source matrix", S_{ie}=1 iff i is the source of e.
	T = Array{Float64,2}(undef,n,0) # "target matrix", T_{ie}=1 iff i is the target of e.
	w = Array{Float64,1}()
	
	for i in 1:n
		for j in 1:n
			if (L[i,j] != 0.) && (i != j)
				s = zeros(n)
				t = zeros(n)
				s[i] = 1.
				t[j] = 1.
				S = [S s]
				T = [T t]
				push!(w,-L[i,j])
			end
		end
	end

	B = S - T

	return B,S,T,w
end




