using DelimitedFiles, LinearAlgebra, PyPlot

include("kerrange.jl")
include("ksakaguchi.jl")
include("tools.jl")

function iterations(f0::Array{Float64,1}, Bout::Array{Float64,2}, b::Array{Float64,2}, ω::Array{Float64,1}, u::Array{Int64,1}, Lmin::Array{Float64,1}, γ::Float64, λ::Float64=.01, T::Int64=100, plot::Bool=true)
	n,m = size(Bout)
	m2 = Int(m/2)

	hγm = h(-γ)
	hγp = h(γ)

	P = cutset_proj(b,Lmin)

	f = πω(f0,Bout,ω,(hγm,hγp))

	f1 = Array{Float64,2}(undef,m,0)
	f2 = Array{Float64,2}(undef,m,0)
	Δ1 = Array{Float64,2}(undef,m,0)
	Δ2 = Array{Float64,2}(undef,m2,0)
	ff = Array{Float64,2}(undef,m,0)
	ff = [ff f0]
	ΔΔ = Array{Float64,2}(undef,m,0)

	for t in 1:T
		@info "iter = $t"

		f1 = [f1 f]
		ff = [ff f]

		Δ = dcc(hi(f))
		Δ1 = [Δ1 Δ]
		ΔΔ = [ΔΔ Δ]

		Δ = dcc(Tu(Δ[1:m2],-Δ[m2+1:m],u,P,C,λ))
		Δ2 = [Δ2 Δ]
		ΔΔ = [ΔΔ [Δ;-Δ]]

		f = [h(Δ);h(-Δ)]
		f2 = [f2 f]
		ff = [ff f]

		f = πω(f,Bout,ω,(hγm,hγp))
	end

	if plot
		i0 = 1
		Δs = LinRange(-π/2+α,π/2-α,200)
		x = Array{Float64,1}()
		X = Array{Float64,1}()
		push!(x,f0[i0])
		push!(X,f0[i0])
		y = Array{Float64,1}()
		Y = Array{Float64,1}()
		push!(y,f0[i0+m2])
		push!(Y,f0[i0+m2])
		xx = Array{Float64,1}()
		XX = Array{Float64,1}()
		yy = Array{Float64,1}()
		YY = Array{Float64,1}()

		figure()
		
		subplot(1,2,1)
		PyPlot.plot(h(Δs),h(-Δs),"--k")
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C8",markersize=20.)
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C1",markersize=13.)
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C3",markersize=5.)
		PyPlot.plot(f0[i0],f0[i0+m2],"ok")

		subplot(1,2,2)
		PyPlot.plot([-π/2+α,π/2-α],[π/2-α,-π/2+α],"--k")
		PyPlot.plot(Δf[i0],-Δf[i0],"o",color="C8",markersize=20.)
		PyPlot.plot(Δf[i0],-Δf[i0],"o",color="C1",markersize=13.)
		PyPlot.plot(Δf[i0],-Δf[i0],"o",color="C3",markersize=5.)

		for t in 1:T
			x = [x;f1[i0,t];f2[i0,t]]
			push!(X,f1[i0,t])
			y = [y;f1[i0+m2,t];f2[i0+m2,t]]
			push!(Y,f1[i0+m2,t])

			xx = [xx;Δ1[i0,t];Δ2[i0,t]]
			push!(XX,Δ1[i0,t])
			yy = [yy;Δ1[i0+m2,t];-Δ2[i0,t]]
			push!(YY,Δ1[i0+m2,t])
		end

		subplot(1,2,1)
		PyPlot.plot(x,y,"-",color="C7")
		PyPlot.plot(X,Y,"-k")
		
		xlabel("f_e")
		ylabel("f_{\bar{e}}")

		subplot(1,2,2)
		PyPlot.plot(xx,yy,"-",color="C7")
		PyPlot.plot(XX,YY,"-k")

		xlabel("Δ_{ij}")
		ylabel("Δ_{ji}")
		
		for t in 1:T
			colo = "C$(mod(t-1,10))"

			subplot(1,2,1)
			PyPlot.plot(f1[i0,t],f1[i0+m2,t],"o",color=colo)
			PyPlot.plot(f2[i0,t],f2[i0+m2,t],"o",color=colo)

			subplot(1,2,2)
			PyPlot.plot(Δ1[i0,t],Δ1[i0+m2,t],"o",color=colo)
			PyPlot.plot(Δ2[i0,t],-Δ2[i0,t],"o",color=colo)
		end

	end
	
	return f1,f2,Δ1,Δ2
end

function cutset_proj(b::Array{Float64,2}, Lmin::Array{Float64,1})
	m2 = length(Lmin)

	I2 = diagm(0 => ones(m2))
	L2 = diagm(0 => Lmin)

	return I2 - L2*b'*pinv(b*L2*b')*b
end


function Tu(Δ1::Array{Float64,1}, Δ2::Array{Float64,1}, u::Array{Int64,1}, P::Array{Float64,2}, C::Array{Float64,2}, λ::Float64=.5)
	m2 = length(Δ1)

	Δ = dcc((Δ1 + Δ2)/2)

	return Δ - λ*diagm(0 => ones(m2))*P*(Δ - 2π*pinv(C)*u)
end

function πω(f::Array{Float64,1}, Bout::Array{Float64,2}, ω::Array{Float64,1}, hγ::Tuple{Float64,Float64})
	n,m = size(Bout)

	A = [kerrange(Bout) pinv(Bout)*ones(n)]

	D = Bout'[:,2:end]

	Ω = pinv(Bout)*ω

	Π = A*pinv(A)

	return min.(max.(Π*(f - Ω) + Ω,hγ[1]),hγ[2])
end

