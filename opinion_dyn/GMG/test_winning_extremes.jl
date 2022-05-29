using PyPlot, Statistics

include("tools.jl")

n = 1000
p = 3

n_iter = 1

res = 40 
ϵs = LinRange(0.,.5,res)
#ϵs = [Vector(LinRange(0.,.1,Int(res/2)));Vector(LinRange(.1,.5,Int(res/2)))]
dϵ = ϵs[2] - ϵs[1]

W = zeros(res,0)
V = zeros(Int64,n_iter,res,n)
XX0 = zeros(n_iter,n,p)
XXs = zeros(n_iter,n,p,res)

for iter in 1:n_iter
	global XX0,W,V,n,p,res,ϵs,XXs

	@info "iter: $iter"
	
	X0 = order_opinion(gen_opinion_1(n,p))
#	X0 = gen_opinion_2(n,p)
#	X0 = gen_opinion_3(n,ones(p)*(1/p))
	XX0[iter,:,:] = X0

	ws = Int64[]
	for j in 1:res
		ϵ = ϵs[j]
		L,A,D,M = Lϵ(X0,ϵ)
		X1 = M*X0
		XXs[iter,:,:,j] = X1
		R = election_result(X1)
		r,w = findmax(R)
		push!(ws,w)
		V[iter,j,:] = party_vote(X1)
	end

	W = [W ws]
end

figure()
for k in p:-1:1
	for j in 1:res
		ϵ = ϵs[j]
		ρ = sum(W[j,:] .<= k)/n_iter
		PyPlot.fill(dϵ/3*[-1,1,1,-1] .+ ϵ,[0.,0.,ρ,ρ],color="C$k")
	end
end

 #= 
figure()
for i in 1:n_iter
	for k in 1:100
		subplot(1,p,V[i,1,k])
		PyPlot.plot(ϵs,V[i,:,k],color="C$(V[i,1,k])",alpha=maximum(XX0[i,k,:]))
	end
end
# =#

 #=
Cp = zeros(n_iter,res-1,p)
Cm = zeros(n_iter,res-1,p)
C = zeros(n_iter,res-1,p)

for i in 1:n_iter
	for k in 1:n
		c = V[i,2:res,k] - V[i,1:res-1,k]
		Cp[i,:,V[i,1,k]] += Int64.(c .> 0)
		Cm[i,:,V[i,1,k]] -= Int64.(c .< 0)
		C[i,:,V[i,1,k]] += c
	end
end



figure()
mCp = zeros(p,res-1)
mCm = zeros(p,res-1)
mC = zeros(p,res-1)
#for i in 1:n_iter
	for k in 1:p
		mCp[k,:] = [mean(Cp[:,j,k]) for j in 1:res-1]
		mCm[k,:] = [mean(Cm[:,j,k]) for j in 1:res-1]
		mC[k,:] = [mean(C[:,j,k]) for j in 1:res-1]

		subplot(3,3,k)
#		PyPlot.plot(ϵs[2:res],Cp[i,:,k],color="C$k")
		PyPlot.plot(ϵs[2:res],mCp[k,:],color="C$k")
		subplot(3,3,k+3)
#		PyPlot.plot(ϵs[2:res],Cm[i,:,k],color="C$k")
		PyPlot.plot(ϵs[2:res],mCm[k,:],color="C$k")
		subplot(3,3,k+6)
#		PyPlot.plot(ϵs[2:res],C[i,:,k],color="C$k")
		PyPlot.plot(ϵs[2:res],mC[k,:],color="C$k")
	end
#end

mm = max(maximum(abs.(mCp)),maximum(abs.(mCm)),maximum(abs.(mC)))
for i in 1:9
	subplot(3,3,i)
	axis([ϵs[1],ϵs[res],-mm,mm])
end
# =#

one2two = Vector{Int64}()
one2thr = Vector{Int64}()
two2one = Vector{Int64}()
two2thr = Vector{Int64}()
thr2one = Vector{Int64}()
thr2two = Vector{Int64}()

for i in 1:n
	for j in 1:res-1
		if V[1,j,i] == 1 && V[1,j+1,i] == 2
			push!(one2two,i)
		elseif V[1,j,i] == 1 && V[1,j+1,i] == 3
			push!(one2thr,i)
		elseif V[1,j,i] == 2 && V[1,j+1,i] == 1
			push!(two2one,i)
		elseif V[1,j,i] == 2 && V[1,j+1,i] == 3
			push!(two2thr,i)
		elseif V[1,j,i] == 3 && V[1,j+1,i] == 1
			push!(thr2one,i)
		elseif V[1,j,i] == 3 && V[1,j+1,i] == 2
			push!(thr2two,i)
		end
	end
end



