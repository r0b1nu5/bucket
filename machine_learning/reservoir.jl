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
	
	@info "Computing optimal post-treatment..."
	rb = sum(rt[:,1:dT:T],dims=2)./length(1:dT:T)
	rb2 = sum(rt[:,1:dT:T].^2,dims=2)./length(1:dT:T)
	sb = sum(st[:,1:dT:T],dims=2)./length(1:dT:T)

	dR = rt[:,1:dT:T] - repeat(rb,1,length(1:dT:T))
	dS = st[:,1:dT:T] - repeat(sb,1,length(1:dT:T))

	Wout = dS*transpose(dR)*inv(Symmetric(dR*transpose(dR) + beta*diagm(0 => ones(N))))
	c = -vec(Wout*rb - sb)
	
	return Wout,c,rt
end

function reservoir_prediction(us::Array{Float64,2},Wout::Array{Float64,2},c::Array{Float64,1},r0::Array{Float64,1},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64)
	T = size(us)[2]

#	@info "Predicting trajectory..."	
	rs = reservoir_tanh(r0,us,A,Win,a,xi)
#	rs = rt1(reservoir_tanh(r0,us,A,Win,a,xi))
	ss = Wout*rs + repeat(c,1,T)
	
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

function reservoir_tanh(r0::Array{Float64,1},us::Array{Float64,2},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64)
	T = size(us)[2]
	
	rs = copy(r0)
	for t in 1:T
		if t%1000 == 0
			@info "$t/$T"
			writedlm("data1/rs_$(t).csv",rs,',')
			rs = rs[:,end]
		end

		rs = [rs ((1-a)*rs[:,end] + a*tanh.(A*rs[:,end] + Win*us[:,t] .+ xi))]
	end

	RS = Array{Float64,2}(undef,length(r0),0)
	for t in 1000:1000:T
		RS = [RS readdlm("data1/rs_$(t).csv",',')[:,2:end]]
		rm("data1/rs_$(t).csv")
	end
	RS = [RS rs[:,2:end]]
	
	return RS
end

#= 
## ======================== DOES NOT WORK ===============================

function reservoir_training_lin(training_data::Tuple{Array{Float64,2},Array{Float64,2}},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64,beta::Float64=.01)
	ut = training_data[1]
	st = training_data[2]
	
	T = size(ut)[2]
	
	@info "Computing reservoir states..."
	r0 = rand(size(A)[1])
	rt = reservoir_lin(r0,ut,A,Win,a,xi)
#	rt = rt1(reservoir_lin(r0,ut,A,Win,a,xi))
	
	@info "Computing optimal post-treatment..."
	rb = sum(rt,dims=2)./T
	sb = sum(st,dims=2)./T
	
	dR = rt[:,2:end] - repeat(rb,1,T)
	dS = st - repeat(sb,1,T)
	
	Wout = dS*transpose(dR)*inv(dR*transpose(dR) + beta*diagm(0 => ones(size(A)[1])))
	c = -vec(Wout*rb - sb)
	
	return Wout,c,rt
end

function reservoir_prediction_lin(us::Array{Float64,2},Wout::Array{Float64,2},c::Array{Float64,1},r0::Array{Float64,1},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64)
	T = size(us)[2]

	@info "Predicting trajectory..."	
	rs = reservoir_lin(r0,us,A,Win,a,xi)
	ss = Wout*rs[:,2:end] + repeat(c,1,T)
	
	return ss,rs
end

function reservoir_lin(r0::Array{Float64,1},u::Array{Float64,1},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64)
	r1 = (1-a)*r0 + a*min(1,max(-1,(A*r0 + Win*u .+ xi)./500))
	
	return [r0 r1]
end

function reservoir_lin(r0::Array{Float64,1},us::Array{Float64,2},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64)
	T = size(us)[2]
	
	rs = copy(r0)
	for t in 1:T
		rs = [rs ((1-a)*rs[:,end] + a*min.(1,max.(-1,(A*rs[:,end] + Win*us[:,t] .+ xi)./500)))]
	end
	
	return rs
end
=#

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



