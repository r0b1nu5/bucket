using PyPlot, LinearAlgebra

# n_year: number of iterations of the generation process
# ppyear: number of papers published in at each iteration
# ρ0: average proportion of papers by new authors at each iteration
# amin: age at which the number of publication start decreasing
# amax: age at which the authors do not publish anymore

function generate_distr(n_year::Int64, ppyear::Int64, ρ0::Float64, amin::Int64=40, amax::Int64=60)
	authors = Dict{Int64,Any}() # contains the list of all authors with their academic age and the number of papers published
	p0 = floor(Int64,ppyear/6)
	for i in 1:p0
		authors[i] = Dict{String,Int64}()
		authors[i]["age"] = 1
		authors[i]["npaper"] = 1
		authors[i+p0] = Dict{String,Int64}()
		authors[i+p0]["age"] = 1
		authors[i+p0]["npaper"] = 2
		authors[i+2*p0] = Dict{String,Int64}()
		authors[i+2*p0]["age"] = 1
		authors[i+2*p0]["npaper"] = 3
	end

	for y in 1:n_year
		@info "year: $y"

		ρ = ρ0 + (.02*rand() - .01)
#		n_new = round(Int64,ρ*ppyear)
#		n_old = ppyear - n_new
		l = length(authors)

		Ns,authks = get_Ns(authors)

		γ = gamma(Ns)

		for k in 1:length(Ns)
			Tk = get_Tk(k,Ns[k],γ,ppyear)

			ages = get_ages(authors,authks[k])
			
			if length(ages) > 0
				ids = rand_age(Tk,ages,ρ,amin,amax)
				
				for i in ids
					if i == 0
						authors[length(authors)+1] = Dict{String,Any}("age" => 1, "npaper" => 1)
					else
						authors[authks[k][i]]["npaper"] += 1
					end
				end
			end
		end

		for a in keys(authors)
			authors[a]["age"] += 1
		end
		
#		for i in 1:n_new
#			authors[l+i] = Dict{String,Any}("age" => 1, "npaper" => 1)
#		end
	end

	return [authors[a]["npaper"] for a in keys(authors)]
	
end

# returns the number of authors with each number of papers.
# Ns is a vector whose k-th component is the number of authors with k papers
# authks is a vector whose k-th component is the list of authors with k papers
function get_Ns(auths::Dict{Int64,Any})
	l = maximum([auths[i]["npaper"] for i in keys(auths)])

	Ns = zeros(Int64,l)
	authks = [Array{Int64,1}() for i in 1:l]

	for a in keys(auths)
		np = auths[a]["npaper"]
		Ns[np] += 1
		push!(authks[np],a)
	end

	return Ns,authks
end

# returns the vector the vector of ages of the authors listed in authk
function get_ages(auth::Dict{Int64,Any}, authk::Array{Int64,1})
	return [auth[authk[i]]["age"] for i in 1:length(authk)]
end

# returns the normalization factor such that sum_k gamma*k*Tk = 1
function gamma(Ns::Array{Int64,1})
	return 1 ./ sum([k*Ns[k] for k in 1:length(Ns)])
end

# returns the number of papers published by the the group of authors with k papers.
function get_Tk(k::Int64, Nk::Int64, γ::Float64, T::Int64=100)
	return round(Int64,γ*k*Nk*T)
end

function share(a::Int64, amin::Int64=40, amax::Int64=60)
	return max(min(amax-a,amax-amin),0)
end

function share(a::Array{Int64,1}, amin::Int64=40, amax::Int64=60)
	return [share(a[i],amin,amax) for i in 1:length(a)]
end

# atributes p new paper to an author considering their age. 
# with proba ρ, the paper is atributed to a new author.
function rand_age(p::Int64, ages::Array{Int64,1}, ρ::Float64, amin::Int64=40, amax::Int64=60)
	sh = share(ages,amin,amax)
	ids = Array{Int64,1}()
	for i in 1:length(ages)
		ids = [ids;i*ones(Int64,sh[i])]
	end

	n_new = floor(Int64,ρ*length(ids))+1
	ids = [zeros(Int64,n_new);ids]

	x = rand(p)
	ix = ceil.(Int64,x*length(ids))

	return ids[ix]
end

