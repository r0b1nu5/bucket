include("ks.jl")
include("tools.jl")

n = 20

# #=
B = zeros(2,2)
while abs(eigvals(B*B')[end-1]) < 1e-6
    global B = get_incidence(get_rand_graph(n,.5))
end
m = size(B)[2]
# =#

ω = rand(n); ω .-= mean(ω)
α = 3.

θs,dθs = ksakaguchi(2π*rand(n),B,ones(n),ω,α)

J = get_jac(θs[:,end],B,ones(m),α)
λs = eigvals(J)

figure("fig")
subplot(2,1,1)
PyPlot.plot(real.(λs),imag.(λs),"o")
xlabel("Re(λ)")
ylabel("Im(λ)")

subplot(2,1,2)
for i in 1:n
    PyPlot.plot(θs[i,:])
end



