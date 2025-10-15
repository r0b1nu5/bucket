using LinearAlgebra

function get_jac(θ::Vector{Float64}, B::Matrix{Float64}, w::Vector{Float64}, α::Float64)
    n = length(θ)

    Δ = repeat(θ,1,n) - repeat(θ',n,1)

    J = -(B*diagm(0 => w)*B').*(1 .- diagm(0 => ones(n))).*cos.(Δ .- α)
    J -= diagm(0 => J*ones(n))

    return J
end

function get_rand_graph(n::Int64,p::Float64)
    A = Int64.(rand(n,n) .< p)
    for i in 1:n-1
        for j in i+1:n
            A[j,i] = A[i,j]
        end
        A[i,i] = 0
    end
    A[n,n] = 0

    return A
end

function get_incidence(A::Matrix{Int64})
    n = size(A)[1]
    B = zeros(n,0)
    for i in 1:n-1
        for j in i+1:n
            if A[i,j] == 1
                x = zeros(n)
                x[i] = 1
                x[j] = -1
                B = [B x]
            end
        end
    end

    return B
end


