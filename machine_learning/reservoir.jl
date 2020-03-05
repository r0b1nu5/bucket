using PyPlot,LinearAlgebra, DelimitedFiles


function reservoir_training(training_data::Tuple{Array{Float64,2},Array{Float64,2}},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64,beta::Float64=.01)
	ut = training_data[1]
	st = training_data[2]
	
	T = size(ut)[2]
	
	@info "Computing reservoir states..."
	r0 = rand(size(A)[1])
	rt = reservoir_tanh(r0,ut,A,Win,a,xi)
#	rt = rt1(reservoir_tanh(r0,ut,A,Win,a,xi))
	
	@info "Computing optimal post-treatment..."
	rb = sum(rt,dims=2)./T
	sb = sum(st,dims=2)./T

	dR = rt[:,1:end] - repeat(rb,1,T)
	dS = st - repeat(sb,1,T)
	
	Wout = dS*transpose(dR)*inv(dR*transpose(dR) + beta*diagm(0 => ones(size(A)[1])))
	c = -vec(Wout*rb - sb)
	
	return Wout,c,rt
end

function reservoir_prediction(us::Array{Float64,2},Wout::Array{Float64,2},c::Array{Float64,1},r0::Array{Float64,1},A::Array{Float64,2},Win::Array{Float64,2},a::Float64,xi::Float64)
	T = size(us)[2]

	@info "Predicting trajectory..."	
	rs = reservoir_tanh(r0,us,A,Win,a,xi)
#	rs = rt1(reservoir_tanh(r0,us,A,Win,a,xi))
	ss = Wout*rs[:,2:end] + repeat(c,1,T)
	
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
	A .*= rho/maximum(abs.(A))
	
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



