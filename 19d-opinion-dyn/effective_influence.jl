using PyPlot, SpecialFunctions

function influence_balance(x::Float64, ϵ::Float64, p::Float64=.5, Δ::Float64=0., μ::Float64=0., σ::Float64=.2)
	left = p/2*(
		    erf((x-μ-Δ)/(sqrt(2)*σ)) - erf((x-ϵ-μ-Δ)/(sqrt(2)*σ))
		    ) 
	+ (1-p)/2*(
		   erf((x-μ+Δ)/(sqrt(2)*σ)) - erf((x-ϵ-μ+Δ)/(sqrt(2)*σ))
		   )

	right = p/2*(
		     erf((x+ϵ-μ-Δ)/(sqrt(2)*σ)) - erf((x-μ-Δ)/(sqrt(2)*σ))
		     ) 
	+ (1-p)/2*(
		   erf((x+ϵ-μ+Δ)/(sqrt(2)*σ)) - erf((x-μ+Δ)/(sqrt(2)*σ))
		   )

	return left,right
end

function influence_balance(xs::Vector{Float64}, ϵ::Float64, p::Float64=.5, Δ::Float64=0., μ::Float64=0., σ::Float64=.2)
	l = Float64[]
	r = Float64[]
	for i in 1:length(xs)
		lr = influence_balance(xs[i],ϵ,p,Δ,μ,σ)
		push!(l,lr[1])
		push!(r,lr[2])
	end

	return l,r
end

function influence_balance(x::Float64, ϵs::Vector{Float64}, p::Float64=.5, Δ::Float64=0., μ::Float64=0., σ::Float64=.2)
	l = Float64[]
	r = Float64[]
	for i in 1:length(ϵs)
		lr = influence_balance(x,ϵs[i],p,Δ,μ,σ)
		push!(l,lr[1])
		push!(r,lr[2])
	end
	
	return l,r
end

function influence_balance(xs::Vector{Float64}, ϵs::Vector{Float64}, p::Float64=.5, Δ::Float64=0., μ::Float64=0., σ::Float64=.2)
	L = zeros(length(xs),0)
	R = zeros(length(xs),0)
	for e in ϵs
		LR = influence_balance(xs,e,p,Δ,μ,σ)
		L = [L LR[1]]
		R = [R LR[2]]
	end

	return L,R
end


function influence_balance_2nd(x::Float64, ϵ::Float64, p::Float64=.5, Δ::Float64=0., μ::Float64=0., σ::Float64=.2)
	left = (
		p*(
		   .5*(x-μ-Δ)*(
			       erf((x-μ-Δ)/(sqrt(2)*σ)) - erf((x-ϵ-μ-Δ)/(sqrt(2)*σ))
			       )
		   + σ/(sqrt(2*π))*(
				    exp(-(x-μ-Δ)^2/(2*σ^2)) - exp(-(x-ϵ-μ-Δ)^2/(2*σ^2))
				    )
		   )
		+ (1-p)*(
			 .5*(x-μ+Δ)*(
				     erf((x-μ+Δ)/(sqrt(2)*σ)) - erf((x-ϵ-μ+Δ)/(sqrt(2)*σ))
				     )
			 + σ/(sqrt(2*π))*(
					  exp(-(x-μ+Δ)^2/(2*σ^2)) - exp(-(x-ϵ-μ+Δ)^2/(2*σ^2))
					  )
			 )
		)

	right = (
		 p*(
		    .5*(x-μ-Δ)*(
				erf((x-μ-Δ)/(sqrt(2)*σ)) - erf((x+ϵ-μ-Δ)/(sqrt(2)*σ))
				)
		    + σ/(sqrt(2*π))*(
				     exp(-(x-μ-Δ)^2/(2*σ^2)) - exp(-(x+ϵ-μ-Δ)^2/(2*σ^2))
				     )
		    )
		 + (1-p)*(
			  .5*(x-μ+Δ)*(
				      erf((x-μ+Δ)/(sqrt(2)*σ)) - erf((x+ϵ-μ+Δ)/(sqrt(2)*σ))
				      )
			  + σ/(sqrt(2*π))*(
					   exp(-(x-μ+Δ)^2/(2*σ^2)) - exp(-(x+ϵ-μ+Δ)^2/(2*σ^2))
					   )
			  )
		 )
	
	return left, right
end

function influence_balance_2nd(xs::Vector{Float64}, ϵ::Float64, p::Float64=.5, Δ::Float64=0., μ::Float64=0., σ::Float64=.2)
	L = Float64[]
	R = Float64[]
	for x in xs
		LR = influence_balance_2nd(x,ϵ,p,Δ,μ,σ)
		push!(L,LR[1])
		push!(R,LR[2])
	end
	
	return L,R
end

function influence_balance_2nd(x::Float64, ϵs::Vector{Float64}, p::Float64=.5, Δ::Float64=0., μ::Float64=0., σ::Float64=.2)
	L = Float64[]
	R = Float64[]
	for ϵ in ϵs
		LR = influence_balance_2nd(x,ϵ,p,Δ,μ,σ)
		push!(L,LR[1])
		push!(R,LR[2])
	end
	
	return L,R
end

function influence_balance_2nd(xs::Vector{Float64}, ϵs::Vector{Float64}, p::Float64=.5, Δ::Float64=0., μ::Float64=0., σ::Float64=.2)
	L = zeros(length(xs),0)
	R = zeros(length(xs),0)
	for ϵ in ϵs
		LR = influence_balance_2nd(xs,ϵ,p,Δ,μ,σ)
		L = [L LR[1]]
		R = [R LR[2]]
	end
	
	return L,R
end


function plot_balance(L::Matrix{Float64}, R::Matrix{Float64}, xs::Vector{Float64}, ϵs::Vector{Float64}, cm::String="RdBu", res::Int64=50)
	cmp = get_cmap(cm)

	contourf(xs,ϵs,(R-L)',res,cmap=cmp)
	
	colorbar(label="right-left")
	xlabel("x")
	ylabel("ϵ")

	return nothing
end




