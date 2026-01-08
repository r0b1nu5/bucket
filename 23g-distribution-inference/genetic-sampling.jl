using Random, Statistics
# Extract a subsamping of the `x` and `y`, composed of `N` components, and that minimizes the correlation between the two time series, using a genetic algorithm.

function genetic_sampling(x::Vector{Float64}, y::Vector{Float64}, N::Int64, n_parents::Int64=10, n_offsprings::Int64=10, ρ::Float64=.05, n_generation::Int64=100)
    Ntot = min(length(x),length(y))
    if N > Ntot
        @info "N too large."
    end

    os = [shuffle(1:Ntot)[1:N] for i in 1:n_offsprings]
    rs = [cor(x[o],y[o]) for o in os]
    ids = sortperm(abs.(rs))[1:n_parents]
    ps = os[ids]

    Rmax = [maximum(abs.(rs[ids])),]
    Rmin = [minimum(abs.(rs[ids])),]
    gen = 0
    while gen < n_generation
        gen += 1
        @info "gen = $gen"
        os = Vector{Int64}[]

        for i in 1:n_parents
            for j in 1:n_parents
#                @info "i,j = $i,$j"
                p1 = ps[i]
                p2 = ps[j]
                
                i0 = intersect(p1,p2)
                i1 = setdiff(p1,p2)
                i2 = setdiff(p2,p1)

                r = rand(length(i1))
                i12 = sort([i0;i1[r .< .5];i2[r .>= .5]])

                for i in 1:N
                    if rand() < ρ
                        i12[i] = rand(setdiff(1:Ntot,i12))
                    end
                end

                push!(os,i12)
            end
        end

        rs = [cor(x[o],y[o]) for o in os]
        ids = sortperm(abs.(rs))[1:n_parents]
        ps = os[ids]
        push!(Rmax,maximum(abs.(rs[ids])))
        push!(Rmin,minimum(abs.(rs[ids])))
    end

    return ps, Rmax, Rmin
end












