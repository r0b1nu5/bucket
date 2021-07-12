using DelimitedFiles, LinearAlgebra, PyPlot

include("kerrange.jl")
include("ksakaguchi.jl")
include("tools.jl")

function iterations1(f0::Array{Float64,1}, Bout::Array{Float64,2}, b::Array{Float64,2}, ω::Array{Float64,1}, u::Array{Int64,1}, Lmin::Array{Float64,1}, γ::Float64, λ::Float64=.01, T::Int64=100, plot::Bool=true)
	n,m = size(Bout)
	m2 = Int(m/2)

	hγm = h(-γ)
	hγp = h(γ)

	P = cycle_proj(b,Lmin)

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

function iterations2(Δ0::Array{Float64,1}, Bout::Array{Float64,2}, b::Array{Float64,2}, ω::Array{Float64,1}, u::Array{Int64,1}, Lmin::Array{Float64,1}, γ::Float64, λ::Float64=.01, T::Int64=100, plot::Bool=true)
	n,m = size(Bout)
	m2 = Int(m/2)

	hγm = h(-γ)
	hγp = h(γ)

	P = cycle_proj(b,Lmin)

	Δ = Δ0

	f1 = Array{Float64,2}(undef,m,0)
	f2 = Array{Float64,2}(undef,m,0)
	Δ1 = Array{Float64,2}(undef,m,0)
	Δ2 = Array{Float64,2}(undef,m2,0)
	ff = Array{Float64,2}(undef,m,0)
	ΔΔ = Array{Float64,2}(undef,m,0)
	ΔΔ = [ΔΔ [Δ0;-Δ0]]

	for t in 1:T
		@info "iter = $t"

		Δ2 = [Δ2 Δ]
		ΔΔ = [ΔΔ [Δ;-Δ]]

		f = [h(Δ);h(-Δ)]
		f2 = [f2 f]
		ff = [ff f]

		f = πω(f,Bout,ω,(hγm,hγp))
		f1 = [f1 f]
		ff = [ff f]

		Δ = dcc(hi(f))
		Δ1 = [Δ1 Δ]
		ΔΔ = [ΔΔ Δ]

		Δ = dcc(Tu(Δ[1:m2],-Δ[m2+1:m],u,P,C,λ))
	end

	if plot
		i0 = 1
		Δs = LinRange(-π/2+α,π/2-α,200)
		x = Array{Float64,1}()
		X = Array{Float64,1}()
		y = Array{Float64,1}()
		Y = Array{Float64,1}()
		xx = Array{Float64,1}()
		XX = Array{Float64,1}()
		push!(xx,Δ0[i0])
		push!(XX,Δ0[i0])
		yy = Array{Float64,1}()
		YY = Array{Float64,1}()
		push!(yy,-Δ0[i0])
		push!(YY,-Δ0[i0])

		figure()
		
		subplot(1,2,1)
		PyPlot.plot(h(Δs),h(-Δs),"--k")
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C8",markersize=20.)
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C1",markersize=13.)
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C3",markersize=5.)

		subplot(1,2,2)
		PyPlot.plot([-π/2+α,π/2-α],[π/2-α,-π/2+α],"--k")
		PyPlot.plot(Δf[i0],-Δf[i0],"o",color="C8",markersize=20.)
		PyPlot.plot(Δf[i0],-Δf[i0],"o",color="C1",markersize=13.)
		PyPlot.plot(Δf[i0],-Δf[i0],"o",color="C3",markersize=5.)
		PyPlot.plot(Δ0[i0],-Δ0[i0],"ok")

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

function iterations3(f0::Array{Float64,1}, Bout::Array{Float64,2}, B::Array{Float64,2}, ω::Array{Float64,1}, u::Array{Int64,1}, Lmin::Array{Float64,1}, γ::Float64, λ::Float64=.01, T::Int64=100, plot::Bool=true)
	n,m = size(Bout)
	m2 = Int(m/2)

	hγm = h(-γ)
	hγp = h(γ)

	P = dir_cycle_proj(B,Lmin)

	f = πω(f0,Bout,ω,(hγm,hγp))

	f1 = Array{Float64,2}(undef,m,0)
	f2 = Array{Float64,2}(undef,m,0)
	f3 = Array{Float64,2}(undef,m,0)
	ff = Array{Float64,2}(undef,m,0)
	ff = [ff f0]

	for t in 1:T
		@info "iter = $t"

		f1 = [f1 f]
		ff = [ff f]

		f = Tu_dir(f,u,P,C,λ)
		f2 = [f2 f]
		ff = [ff f]

		Δ = dcc(hi(f[1:m2]))
		f = [h(Δ);h(-Δ)]
		f3 = [f3 f]
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

		figure()
		
#		subplot(1,2,1)
		PyPlot.plot(h(Δs),h(-Δs),"--k")
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C8",markersize=20.)
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C1",markersize=13.)
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C3",markersize=5.)
		PyPlot.plot(f0[i0],f0[i0+m2],"ok")

		for t in 1:T
			x = [x;f1[i0,t];f2[i0,t];f3[i0,t]]
			push!(X,f1[i0,t])
			y = [y;f1[i0+m2,t];f2[i0+m2,t];f3[i0+m2,t]]
			push!(Y,f1[i0+m2,t])
		end

#		subplot(1,2,1)
		PyPlot.plot(x,y,"-",color="C7")
		PyPlot.plot(X,Y,"-k")
		
		xlabel("f_e")
		ylabel("f_{\bar{e}}")
	
		for t in 1:T
			colo = "C$(mod(t-1,10))"

#			subplot(1,2,1)
			PyPlot.plot(f1[i0,t],f1[i0+m2,t],"o",color=colo)
			PyPlot.plot(f2[i0,t],f2[i0+m2,t],"x",color=colo)
			PyPlot.plot(f3[i0,t],f3[i0+m2,t],"o",color=colo)
		end

		title("u = $u")

	end
	
	return f1,f2,f3
end

function iterations4(f0::Array{Float64,1}, θf::Array{Float64,1}, Bout::Array{Float64,2}, B::Array{Float64,2}, ω::Array{Float64,1}, u::Array{Int64,1}, Lmin::Array{Float64,1}, γ::Float64, λ::Float64=.01, T::Int64=100, plot::Bool=true)
	n,m = size(Bout)
	m2 = Int(m/2)

	hγm = h(-γ)
	hγp = h(γ)

	P = dir_cycle_proj(B,Lmin)

	f = πω(f0,Bout,ω,(hγm,hγp))

	ff = Array{Float64,2}(undef,m,0)
#	ff = [ff f]

#	dmax = Array:{Float64,1}()
#	dmin = Array{Float64,1}()

	for t in 1:T
		@info "iter = $t"

		ff = [ff f]
		f = Tu_dir(f,u,P,C,λ*diagm(0 => ones(m)))
		@info "$(maximum(f)/hγp):$(minimum(f)/(-hγm))"
		
#		dx = 1 ./(sqrt.(1 .- (f .- sin(.1)).^2))
#		push!(dmax,maximum(dx))
#		push!(dmin,minimum(dx))
	end

	if plot
		i0 = 1
		Δf = B'*θf
		Δs = LinRange(-π/2+α,π/2-α,200)
		X = Array{Float64,1}()
		push!(X,f0[i0])
		Y = Array{Float64,1}()
		push!(Y,f0[i0+m2])

		figure()
		
#		subplot(1,2,1)
		PyPlot.plot(h(Δs),h(-Δs),"--k")
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C8",markersize=20.)
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C1",markersize=13.)
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C3",markersize=5.)
		PyPlot.plot(f0[i0],f0[i0+m2],"ok")

		for t in 1:T
			push!(X,ff[i0,t])
			push!(Y,ff[i0+m2,t])
		end

#		subplot(1,2,1)
		PyPlot.plot(X,Y,"-k")
		
		xlabel("f_e")
		ylabel("f_{\bar{e}}")
	
		PyPlot.plot(f0[i0],f0[i0+m2],"ok")

		for t in 1:T
			colo = "C$(mod(t-1,10))"

#			subplot(1,2,1)
			PyPlot.plot(ff[i0,t],ff[i0+m2,t],"o",color=colo)
		end

#=
		figure()
		PyPlot.plot(1:length(dmax),(dmax+dmin)./(dmax-dmin),"o")
=#

		title("u = $u")
	end	

	return ff
end

function iterations4_fail(f0::Array{Float64,1}, θf::Array{Float64,1}, Bout::Array{Float64,2}, B::Array{Float64,2}, ω::Array{Float64,1}, u::Array{Int64,1}, Lmin::Array{Float64,1}, γ::Float64, λ::Float64=.01, T::Int64=100, plot::Bool=true)
	n,m = size(Bout)
	m2 = Int(m/2)

	hγm = h(-γ)
	hγp = h(γ)

	P = dir_cycle_proj(B,Lmin)

#	f = πω(f0,Bout,ω,(hγm,hγp))
	f = f0

	ff = Array{Float64,2}(undef,m,0)
#	ff = [ff f]

#	dmax = Array:{Float64,1}()
#	dmin = Array{Float64,1}()

	for t in 1:T
		@info "iter = $t"

		ff = [ff f]
		f = Tu_dir(f,u,P,C,λ*diagm(0 => ones(m)))
		@info "$(maximum(f)/hγp):$(minimum(f)/(-hγm))"
		
#		dx = 1 ./(sqrt.(1 .- (f .- sin(.1)).^2))
#		push!(dmax,maximum(dx))
#		push!(dmin,minimum(dx))
	end

	if plot
		i0 = 1
		Δf = B'*θf
		Δs = LinRange(-π/2+α,π/2-α,200)
		X = Array{Float64,1}()
		push!(X,f0[i0])
		Y = Array{Float64,1}()
		push!(Y,f0[i0+m2])

		figure()
		
#		subplot(1,2,1)
		PyPlot.plot(h(Δs),h(-Δs),"--k")
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C8",markersize=20.)
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C1",markersize=13.)
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C3",markersize=5.)
		PyPlot.plot(f0[i0],f0[i0+m2],"ok")

		for t in 1:T
			push!(X,ff[i0,t])
			push!(Y,ff[i0+m2,t])
		end

#		subplot(1,2,1)
		PyPlot.plot(X,Y,"-k")
		
		xlabel("f_e")
		ylabel("f_{\bar{e}}")
	
		PyPlot.plot(f0[i0],f0[i0+m2],"ok")

		for t in 1:T
			colo = "C$(mod(t-1,10))"

#			subplot(1,2,1)
			PyPlot.plot(ff[i0,t],ff[i0+m2,t],"o",color=colo)
		end

#=
		figure()
		PyPlot.plot(1:length(dmax),(dmax+dmin)./(dmax-dmin),"o")
=#

		title("u = $u")
	end	

	return ff
end

function iterations5(Δ0::Array{Float64,1}, θf::Array{Float64,1}, Bout::Array{Float64,2}, B::Array{Float64,2}, C::Array{Float64,2}, ω::Array{Float64,1}, u::Array{Int64,1}, ρ1::Float64, ρ2::Float64, T::Int64=100, plot::Bool=true, verb::Bool=false)
	n,m = size(Bout)
	m2 = Int(m/2)

	hγm = h(-γ)
	hγp = h(γ)

	b = B[:,1:m2]
	P = cycle_proj(b,ones(m2))
#	X = diagm(0 => ones(n)) - ones(n,n)/n
	Cdag = pinv(C)

	Δ = Δ0
	Δs = Array{Float64,2}(undef,m2,0)

	for t in 1:T
		if verb
			@info "iter = $t"
		end

		Δs = [Δs Δ]
		
		Δ = Su(Δ,ω,b,Bout,P,Cdag,u,ρ1,ρ2)
	end

	if plot
		i0 = 1
		Δf = B'*θf
		γs = LinRange(-π/2+α,π/2-α,200)
		x = Array{Float64,1}()
		push!(x,h(Δ0[i0]))
		y = Array{Float64,1}()
		push!(y,h(-Δ0[i0]))

		figure()
		
#		subplot(1,2,1)
		PyPlot.plot(h(γs),h(-γs),"--k")
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C8",markersize=20.)
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C1",markersize=13.)
		PyPlot.plot(h(Δf[i0]),h(-Δf[i0]),"o",color="C3",markersize=5.)

		for t in 1:T
			push!(x,h(Δs[i0,t]))
			push!(y,h(-Δs[i0,t]))
		end

#		subplot(1,2,1)
		PyPlot.plot(x,y,"-k")
		
		xlabel("f_e")
		ylabel("f_{\bar{e}}")
	
		PyPlot.plot(h(Δ0[i0]),h(-Δ0[i0]),"ok")
		for t in 1:T
			colo = "C$(mod(t-1,10))"

#			subplot(1,2,1)
			PyPlot.plot(h(Δs[i0,t]),h(-Δs[i0,t]),"o",color=colo)
		end

#=
		figure()
		PyPlot.plot(1:length(dmax),(dmax+dmin)./(dmax-dmin),"o")
=#

		title("u = $u")
	end	

	return Δs
end



function check_monotonicity(ff1::Array{Float64,2}, ff2::Array{Float64,2}, P::Array{Float64,2})
	bool1 = ff1 .< ff2
	bool2 = dcc(hi(ff1)) .< dcc(hi(ff2))
	bool3 = P*dcc(ff1) .< P*dcc(ff2)

	return (bool1==bool2), (bool1==bool3), (bool2==bool3)
end
#
function cycle_proj(b::Array{Float64,2}, Lmin::Array{Float64,1})
	m2 = length(Lmin)

	I2 = diagm(0 => ones(m2))
	L2 = diagm(0 => Lmin)

	return I2 - L2*b'*pinv(b*L2*b')*b
end

function dir_cycle_proj(B::Array{Float64,2})
	n,m = size(B)

	Im = diagm(0 => ones(m))

	Bout = B.*(B .> 0.)

	return Im - B'*pinv(Bout*B')*Bout
end

function dir_cycle_proj(B::Array{Float64,2}, Lmin::Array{Float64,1})
	m = length(Lmin)

	Im = diagm(0 => ones(m))
	D = diagm(0 => Lmin)

	Bout = B.*(B .> 0.)

	return Im - D*B'*pinv(Bout*D*B')*Bout
end

		
function Su(Δ::Array{Float64,1}, ω::Array{Float64,1}, b::Array{Float64,2}, Bout::Array{Float64,2}, P::Array{Float64,2}, Cdag::Array{Float64,2}, u::Array{Int64,1}, ρ1::Float64, ρ2::Float64)
	return Δ + ρ1*pinv(b)*(ω - Bout*h([Δ;-Δ])) - ρ2*P*(Δ - 2π*Cdag*u)
end

function Tu(Δ1::Array{Float64,1}, Δ2::Array{Float64,1}, u::Array{Int64,1}, P::Array{Float64,2}, C::Array{Float64,2}, λ::Float64=1.)
	m2 = length(Δ1)

	Δ = dcc((Δ1 + Δ2)/2)

	return Δ - λ*diagm(0 => ones(m2))*P*(Δ - 2π*pinv(C)*u)
end

function Tu_dir(f::Array{Float64,1}, u::Array{Int64,1}, P::Array{Float64,2}, C::Array{Float64,2}, D::Array{Float64,2})
	m = length(f)

	Δ = dcc(hi(f))

	Cdu = pinv(C)*u

	return f - P*D*(Δ - 2π*[Cdu;-Cdu])
end

function πω(f0::Array{Float64,1}, Bout::Array{Float64,2}, ω::Array{Float64,1}, hγ::Tuple{Float64,Float64})
	n,m = size(Bout)
	
	Bod = pinv(Bout)

	a = Bod*ones(n)
	a /= norm(a)

	A = [kerrange(Bout) a]

	Ω = Bod*ω
	O = Bod*ones(n)

	Π = A*pinv(A)

	f1 = Π*(f0 - Ω)
	inf1 = max.(hγ[1] .- (f1 + Ω),0.)
	sup1 = max.((f1 + Ω) .- hγ[2],0.)
	is1 = max.(inf1,sup1)
	ma,ima = findmax(is1)

	δ = 1.
	f2 = δ*f1
	δδ = .001
	while sum(hγ[1] .< (f2 + Ω) .< hγ[2]) < m
		δ -= δδ
		f2 = δ*f1
	end

	f3 = f2 + Ω

	return f3
end

function rand_init(ω::Array{Float64,1}, Bout::Array{Float64,2}, hγ::Tuple{Float64,Float64})
	n,m = size(Bout)

	Bod = pinv(Bout)

	A = [kerrange(Bout) Bod*ones(n)]

	Ω = Bod*ω

	n,p = size(A)

	f = 2*A*(2*rand(p) .- 1) + Ω
	
	while sum(hγ[1] .< f .< hγ[2]) < m
		f = 2*A*(2*rand(p) .- 1) + Ω
	end

	return f
end
