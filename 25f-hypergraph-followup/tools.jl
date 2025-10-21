using LinearAlgebra

function adj2inc(AA::Matrix{Float64})
    l,w = size(AA)
    if l != w
        @info "Matrix not square! ($l,$w)"
    end

    B = zeros(l,0)
    w = Float64[]

    for i in 1:l-1
        for j in i+1:l
            if AA[i,j] > 1e-4
                x = can_bas(i,l) - can_bas(j,l)
                B = [B x]
                push!(w,AA[i,j])
            end
        end
    end

    return B
end

function adj2lap(A::Matrix{Float64})
	return diagm(0 => sum(A,dims=2)[:,1]) - A
end

function can_bas(i::Int64, n::Int64)
	v = zeros(n)
	v[i] = 1
	return v
end

function get_jac(θ::Vector{Float64},
		 A::Matrix{Float64})
	n = length(θ)
	m = size(A)[1]

	b = [2 -1 -1;
	     -1 2 -1;
	     -1 -1 2]

	J = zeros(n,n)
	for l in 1:m
		ijk = Int64.(A[l,1:3])
		a = A[l,4]
		
		J[ijk,ijk] -= a*b*diagm(0 => cos.(b*θ[ijk]))
	end

	return J
end

function get_loose_ends(A2l::Matrix{Float64}, A3l::Matrix{Float64})
	n = Int64(maximum(A3l[:,1:3]))

	le = setdiff(collect(1:n),union(Int64.(A2l[:,1]),Int64.(A3l[:,1])))

	return le
end

function hyper2edge(A::Matrix{Float64})
    m = size(A)[2]
    n = Int64(maximum(A[:,1:end-1]))

    a = zeros(n,n)

    for l in 1:m
        i,j,k = Int64.(A[l,1:end-1])
        v = A[l,end]
        a[i,j] += v
        a[j,i] += v
        a[i,k] += v
        a[k,i] += v
        a[j,k] += v
        a[k,j] += v
    end

    b = adj2inc(a)

    return b
end




