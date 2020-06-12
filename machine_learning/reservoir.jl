using PyPlot,LinearAlgebra, DelimitedFiles


function reservoir_training(training_data::Tuple{Array{Float64,2},Array{Float64,2}},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64,dT::Int64=1,beta::Float64=.01)
	ut = training_data[1]
	st = training_data[2]
	
	T = size(ut)[2]
	
	@info "Computing reservoir states..."
	N = size(A)[1]
	r0 = rand(N)
	rt = reservoir_tanh(r0,ut,A,Win,a,xi)
#	rt = rt1(reservoir_tanh(r0,ut,A,Win,a,xi))
	
	nabs_r = Array{Array{Float64,2},1}()
	nabs_u = Array{Array{Float64,2},1}()
	for j in 1:size(ut)[2]
		x = nablas(rt[:,j],ut[:,j],A,Win,xi)
		push!(nabs_r,x[1])
		push!(nabs_u,x[2])
	end
	
	@info "Computing optimal post-treatment..."
	rb = sum(rt[:,1:dT:T],dims=2)./length(1:dT:T)
	rb2 = sum(rt[:,1:dT:T].^2,dims=2)./length(1:dT:T)
	sb = sum(st[:,1:dT:T],dims=2)./length(1:dT:T)

	dR = rt[:,1:dT:T] - repeat(rb,1,length(1:dT:T))
	dS = st[:,1:dT:T] - repeat(sb,1,length(1:dT:T))

	Wout = dS*transpose(dR)*inv(Symmetric(dR*transpose(dR) + beta*diagm(0 => ones(N))))
	c = -vec(Wout*rb - sb)
	
	return Wout,c,rt,nabs_r,nabs_u
end

function reservoir_prediction(us::Array{Float64,2},Wout::Array{Float64,2},c::Array{Float64,1},r0::Array{Float64,1},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64)
	T = size(us)[2]

#	@info "Predicting trajectory..."	
	rs = reservoir_tanh(r0,us,A,Win,a,xi)
#	rs = rt1(reservoir_tanh(r0,us,A,Win,a,xi))
	ss = Wout*rs + repeat(c,1,T)
	
	return ss,rs
end

function reservoir_prediction_self(s0::Array{Float64,1}, r0::Array{Float64,1}, Tp::Int64, Wout::Array{Float64,2}, c::Array{Float64,1}, A::Array{Float64,2}, Win::Array{Float64,2}, a::Float64, xi::Float64, id::Int64=0)
	n = length(s0)
	N = length(r0)

	ss = Array{Float64,2}(undef,n,0)
	ss = [ss s0]
	rs = Array{Float64,2}(undef,N,0)
	rs = [rs r0]

	for t in 1:Tp
		r = reservoir_tanh(rs[:,end],ss[:,[size(ss)[2],]],A,Win,a,xi,id)
		s = Wout*r + c
		rs = [rs r]
		ss = [ss s]
	end

	return ss,rs
end


# Uses the square of the reservoir components N/2:N for the prediction of the idx-th component of the system.

function reservoir_training2(training_data::Tuple{Array{Float64,2},Array{Float64,2}},idx::Int64,A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64,dT::Int64=1,beta::Float64=.01)
	ut = training_data[1]
	st = training_data[2]
	
	T = size(ut)[2]
	n = size(st)[1]
	
	@info "Computing reservoir states..."
	N = size(A)[1]
	r0 = rand(N)
	rt = reservoir_tanh(r0,ut,A,Win,a,xi)
#	rt = rt1(reservoir_tanh(r0,ut,A,Win,a,xi))
	
	@info "Computing optimal post-treatment..."
	rb = sum(rt[:,1:dT:T],dims=2)./length(1:dT:T)
	rb2 = sum(rt[:,1:dT:T].^2,dims=2)./length(1:dT:T)
	sb = sum(st[:,1:dT:T],dims=2)./length(1:dT:T)

	dR = rt[:,1:dT:T] - repeat(rb,1,length(1:dT:T))
	dRR = [rt[1:Int(N/2),1:dT:T] - repeat(rb[1:Int(N/2)],1,length(1:dT:T)); 
	       rt[Int(N/2)+1:N,1:dT:T].^2 - repeat(rb2[Int(N/2)+1:N],1,length(1:dT:T))]
	dS = st[:,1:dT:T] - repeat(sb,1,length(1:dT:T))

	xdi = setdiff(Array(1:n),[idx,])

	Wout_xdi = dS[xdi,:]*transpose(dR)*inv(Symmetric(dR*transpose(dR) + beta*diagm(0 => ones(size(A)[1]))))
	Wout_idx = dS[[idx,],:]*transpose(dRR)*inv(Symmetric(dRR*transpose(dRR) + beta*diagm(0 => ones(size(A)[1]))))
	Wout = [Wout_xdi[1:idx-1,:] zeros(idx-1,N);
		Wout_idx[1:Int(N/2)]' zeros(1,N) Wout_idx[1+Int(N/2):N]';
		Wout_xdi[idx+1:n,:] zeros(n-idx,N)]
	c = -vec(Wout*[rb; rb.^2] - sb)
	
	return Wout,c,rt
end

function reservoir_prediction2(us::Array{Float64,2},Wout::Array{Float64,2},c::Array{Float64,1},idx::Int64,r0::Array{Float64,1},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64)
	T = size(us)[2]

#	@info "Predicting trajectory..."	
	rs = reservoir_tanh(r0,us,A,Win,a,xi)
#	rs = rt1(reservoir_tanh(r0,us,A,Win,a,xi))
	ss = Wout*[rs; rs.^2] + repeat(c,1,T)
	
	return ss,rs
end


# Uses the square of the reservoir components N/2:N for the prediction of the idx-th component of the system. Does not unbias the times series.

function reservoir_training3(training_data::Tuple{Array{Float64,2},Array{Float64,2}},idx::Int64,A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64,dT::Int64=1,beta::Float64=.01)
	ut = training_data[1]
	sttemp = training_data[2]
	
	T = size(ut)[2]
	n = size(sttemp)[1]
	
	@info "Computing reservoir states..."
	N = size(A)[1]
	r0 = rand(N)
	rttemp = reservoir_tanh(r0,ut,A,Win,a,xi)
#	rt = rt1(reservoir_tanh(r0,ut,A,Win,a,xi))
	
	@info "Computing optimal post-treatment..."
	rt = rttemp[:,1:dT:T]
	rrt = [rttemp[1:Int(N/2),1:dT:T]; 
	       rttemp[Int(N/2)+1:N,1:dT:T].^2]
	st = sttemp[:,1:dT:T]

	xdi = setdiff(Array(1:n),[idx,])

	Wout_xdi = st[xdi,:]*transpose(rt)*inv(Symmetric(rt*transpose(rt) + beta*diagm(0 => ones(size(A)[1]))))
	Wout_idx = st[[idx,],:]*transpose(rrt)*inv(Symmetric(rrt*transpose(rrt) + beta*diagm(0 => ones(size(A)[1]))))
	Wout = [Wout_xdi[1:idx-1,:] zeros(idx-1,N);
		Wout_idx[[1,],1:Int(N/2)] zeros(1,N) Wout_idx[[1,],1+Int(N/2):N];
		Wout_xdi[idx+1:n,:] zeros(n-idx,N)]
	
	return Wout,rt
end

function reservoir_prediction3(us::Array{Float64,2},Wout::Array{Float64,2},idx::Int64,r0::Array{Float64,1},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64)
	T = size(us)[2]

#	@info "Predicting trajectory..."	
	rs = reservoir_tanh(r0,us,A,Win,a,xi)
#	rs = rt1(reservoir_tanh(r0,us,A,Win,a,xi))
	ss = Wout*[rs; rs.^2]
	
	return ss,rs
end



function reservoir_tanh(r0::Array{Float64,1},u::Array{Float64,1},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64)
	r1 = (1-a)*r0 + a*tanh.(A*r0 + Win*u .+ xi)
	
	return [r0 r1]
end

function reservoir_tanh(r0::Array{Float64,1},us::Array{Float64,2},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64,id::Int64=0)
	T = size(us)[2]
	
	rs = copy(r0)
	for t in 1:T
		if t%1000 == 0
			@info "$t/$T"
			writedlm("data1/rs_$(t)_$(id).csv",rs,',')
			rs = rs[:,end]
		end

		rs = [rs ((1-a)*rs[:,end] + a*tanh.(A*rs[:,end] + Win*us[:,t] .+ xi))]
	end

	RS = Array{Float64,2}(undef,length(r0),0)
	for t in 1000:1000:T
		RS = [RS readdlm("data1/rs_$(t)_$(id).csv",',')[:,2:end]]
		rm("data1/rs_$(t)_$(id).csv")
	end
	RS = [RS rs[:,2:end]]
	
	return RS
end

function rt1(r::Array{Float64,1})
	D = length(r)
	d = round(Int,D/2)
	
	return [r[1:d];r[d+1:end].^2]
end

function rt1(r::Array{Float64,2})
	D = size(r)[1]
	d = round(Int,D/2)
	
	return [r[1:d,:];r[d+1:end,:].^2]
end

function nablas(r::Array{Float64,1}, u::Array{Float64,1}, A::Array{Float64,2}, Win::Array{Float64,2}, xi::Float64)
	cosh2 = (cosh.(A*r + Win*u + xi*ones(length(r)))).^2
	
	nablaf_r = A./repeat(cosh2,1,length(r))
	nablaf_u = Win./repeat(cosh2,1,length(u))

	return nablaf_r,nablaf_u
end


function A_gen(n::Int64,m::Int64,rho::Float64)
	A = zeros(n,n)
	p = (2*m)/(n*(n-1))
	
	for i in 1:n-1
		for j in i:n
			if rand() < p
				A[i,j] = 2*rand()-1
			end
		end
	end
	
	A = A + transpose(A)
	A .*= rho/opnorm(Symmetric(A))
	
	return A
end
	
function Win_gen(n::Int64,N::Int64,sig::Float64)
	Win = zeros(N,n)
	
	p = floor(Int,N/n)
	q = N - n*p
	
	for i in 0:n-1
		Win[((1:p).+p*i),i+1] = sig*(2*rand(p).-1)
	end
	
	for i in 1:q
		Win[n*p+i,rand(1:n)] = sig*(2*rand()-1)
	end
	
	return Win
end


function breaktime(thr::Float64, xs::Array{Float64,2}, ss::Array{Float64,2})
	Tp = size(ss)[2]

	mas = maximum(xs,dims=2)
	mis = minimum(xs,dims=2)
	amps = mas - mis
	diff = xs[:,1:Tp] - ss
	errs = abs.(diff)./repeat(amps,1,Tp)
	merr = [maximum(errs[:,1:i]) for i in 1:Tp]
	
	test = (1:Tp).*(merr[1:end] .> thr)
	if maximum(test) > 0
		Tb = minimum(setdiff(test,[0,]))
	else
		Tb = Tp + 1
	end

	return Tb
end


