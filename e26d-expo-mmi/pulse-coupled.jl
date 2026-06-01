using PyPlot, DelimitedFiles

N = 20
γ = .1 # Saturation rate of each oscillator (see Mirollo and Strogatz (1990))
S0 = .11 # Increase rate of each oscillator (see Mirollo and Strogatz (1990))
ϵ = .005 # Jump size

T = 500
dt = .01
t = 0.
ts = [t,]

A = Float64.(rand(N,N) .> .4)

x = rand(N)
xs = copy(x)

c = 0

while t < T
	global t += dt
	push!(ts,t)
	
	if mod(t,100) < .01
		@info "t = $t"

		global c += 1
		writedlm("temp/xs-$c.csv",xs[:,1:end-1],',')
		global xs = xs[:,end]
	end

	global x += (S0 .- γ*x)*dt

	while maximum(x) > 1.
		local fire = (x .> 1.)
		x[fire] = zeros(sum(fire))

		x += ϵ*A*fire
#		@info "max = $(maximum(ϵ*A*fire))"
	end

	global xs = [xs x]
end

Xs = zeros(N,0)
for i in 1:c
	global Xs = [Xs readdlm("temp/xs-$i.csv",',')]
	rm("temp/xs-$i.csv")
end
Xs = [Xs xs]

for i in 1:N
	PyPlot.plot(ts,Xs[i,:])
end





