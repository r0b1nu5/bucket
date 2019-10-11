include("scripts.jl")

epss = Array(LinRange(0,1,40))

x1 = [rand(Normal(-.5,.2),100);rand(Normal(.5,.2),101)]
es1 = loop_rand(20,x1,epss)
ef1 = loop_fiedler(x1,epss)

x2 = [rand(Normal(-.5,.2),100);rand(Normal(.5,.2),101)]
es2 = loop_rand(20,x2,epss)
ef2 = loop_fiedler(x2,epss)

x3 = [rand(Normal(-.5,.2),100);rand(Normal(.5,.2),101)]
es3 = loop_rand(20,x3,epss)
ef3 = loop_fiedler(x3,epss)

x4 = [rand(Normal(-.5,.2),100);rand(Normal(.5,.2),101)]
es4 = loop_rand(20,x4,epss)
ef4 = loop_fiedler(x4,epss)

x5 = [rand(Normal(-.5,.2),100);rand(Normal(.5,.2),101)]
es5 = loop_rand(20,x5,epss)
ef5 = loop_fiedler(x5,epss)

figure(100)

subplot(1,5,1)
plot_mean(abs.(es1),epss)
plot_fiedler(abs.(ef1),epss)
xlabel("ε")
ylabel("ξ")

subplot(1,5,2)
plot_mean(abs.(es2),epss)
plot_fiedler(abs.(ef2),epss)
xlabel("ε")
ylabel("ξ")

subplot(1,5,3)
plot_mean(abs.(es3),epss)
plot_fiedler(abs.(ef3),epss)
xlabel("ε")
ylabel("ξ")

subplot(1,5,4)
plot_mean(abs.(es4),epss)
plot_fiedler(abs.(ef4),epss)
xlabel("ε")
ylabel("ξ")

subplot(1,5,5)
plot_mean(abs.(es5),epss)
plot_fiedler(abs.(ef5),epss)
xlabel("ε")
ylabel("ξ")




