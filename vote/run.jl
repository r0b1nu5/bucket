include("scripts.jl")


n1 = 1000
n2 = 1001
d = 1.
sig = .2
n_rand = 20
emin = 0.
emax = 1.
ne = 40
epss = Array(LinRange(emin,emax,ne))
fignum = 101

@info "Run 1"
x1 = [rand(Normal(-d/2,sig),n1);rand(Normal(d/2,sig),n2)]
es1 = loop_rand(n_rand,x1,epss)
ef1 = loop_fiedler(x1,epss)
em1 = loop_mini(x1,epss)

@info "Run 2"
x2 = [rand(Normal(-d/2,sig),n1);rand(Normal(d/2,sig),n2)]
es2 = loop_rand(n_rand,x2,epss)
ef2 = loop_fiedler(x2,epss)
em2 = loop_mini(x2,epss)

@info "Run 3"
x3 = [rand(Normal(-d/2,sig),n1);rand(Normal(d/2,sig),n2)]
es3 = loop_rand(n_rand,x3,epss)
ef3 = loop_fiedler(x3,epss)
em3 = loop_mini(x3,epss)

@info "Run 4"
x4 = [rand(Normal(-d/2,sig),n1);rand(Normal(d/2,sig),n2)]
es4 = loop_rand(n_rand,x4,epss)
ef4 = loop_fiedler(x4,epss)
em4 = loop_mini(x4,epss)


figure(fignum)

subplot(1,4,1)
e1 = eps_connect(x1,epss)
xima = max(maximum(abs.(es1)),maximum(abs.(ef1)),maximum(abs.(em1)))
PyPlot.plot([e1,e1],[0,1.1*xima],"--k")
plot_mean(abs.(es1),epss)
plot_fiedler(abs.(ef1),epss)
plot_mini(abs.(em1),epss)
xlabel("ε")
ylabel("ξ")

subplot(1,4,2)
e2 = eps_connect(x2,epss)
xima = max(maximum(abs.(es2)),maximum(abs.(ef2)),maximum(abs.(em2)))
PyPlot.plot([e2,e2],[0,1.1*xima],"--k")
plot_mean(abs.(es2),epss)
plot_fiedler(abs.(ef2),epss)
plot_mini(abs.(em2),epss)
xlabel("ε")
ylabel("ξ")

subplot(1,4,3)
e3 = eps_connect(x3,epss)
xima = max(maximum(abs.(es3)),maximum(abs.(ef3)),maximum(abs.(em3)))
PyPlot.plot([e3,e3],[0,1.1*xima],"--k")
plot_mean(abs.(es3),epss)
plot_fiedler(abs.(ef3),epss)
plot_mini(abs.(em3),epss)
xlabel("ε")
ylabel("ξ")

subplot(1,4,4)
e4 = eps_connect(x4,epss)
xima = max(maximum(abs.(es4)),maximum(abs.(ef4)),maximum(abs.(em4)))
PyPlot.plot([e4,e4],[0,1.1*xima],"--k")
plot_mean(abs.(es4),epss)
plot_fiedler(abs.(ef4),epss)
plot_mini(abs.(em4),epss)
xlabel("ε")
ylabel("ξ")





