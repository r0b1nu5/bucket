using DelimitedFiles, LinearAlgebra, Statistics

function ksakaguchi(θ0::Vector{Float64}, B::Matrix{Float64}, w::Vector{Float64}, ω::Vector{Float64}, α::Union{Float64,Vector{Float64}}; h::Float64=.01, max_iter::Int64=10000, tol::Float64=1e-6)
    n = length(θ0)

    θ = copy(θ0)
    Θs = copy(θ0)
    dΘs = zeros(n,0)

    iter = 0
    err = 1000.
    c = 0

    while iter < max_iter && err > tol
        iter += 1

        k1 = f_ks(θ,B,w,ω,α)
        k2 = f_ks(θ+h/2*k1,B,w,ω,α)
        k3 = f_ks(θ+h/2*k2,B,w,ω,α)
        k4 = f_ks(θ+h*k3,B,w,ω,α)

        dθ = (k1 + 2*k2 + 2*k3 + k4)/6
        θ += dθ*h

        Θs = [Θs θ]
        dΘs = [dΘs dθ]

        if mod(iter,1000) == 0
            c += 1
            writedlm("temp/thetas_$c.csv",Θs[:,1:end-1],',')
            Θs = Θs[:,end]
            writedlm("temp/dhetas_$c.csv",dΘs,',')
            dΘs = zeros(n,0)
        end

        err = maximum(dθ) - minimum(dθ)
    end

    Θ = zeros(n,0)
    dΘ = zeros(n,0)
    for i in 1:c
        Θ = [Θ readdlm("temp/thetas_$i.csv",',')]
        rm("temp/thetas_$i.csv")
        dΘ = [dΘ readdlm("temp/dhetas_$i.csv",',')]
        rm("temp/dhetas_$i.csv")
    end
    Θ = [Θ Θs[:,1:end-1]]
    dΘ = [dΘ dΘs]

    return Θ,dΘ
end

function f_ks(θ::Vector{Float64}, B::Matrix{Float64}, w::Vector{Float64}, ω::Vector{Float64}, α::Union{Float64,Vector{Float64}})
    return ω - diagm(0 => w)*B*(sin.(B'*θ .- α) .+ sin.(α))
end





