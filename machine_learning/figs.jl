using PyPlot, LinearAlgebra, DelimitedFiles, Statistics

Tts = Array(101:100:4001)
T0 = 1001 
Ns = Array(200:200:4000)
N0 = 1000

n = 3

Tbs = Array{Float64,1}()
Wouts = Array{Array{Float64,2},1}()
svds = Array{Float64,2}(undef,n,0)

for i in 1:length(Tts)
	global Tbs,Wouts,svds
	push!(Tbs,readdlm("data1/Tb$(i)_vs_Tt_N2000.csv",',')[1])
	push!(Wouts,readdlm("data1/Wout$(i)_vs_Tt_N2000.csv",','))
	svds = [svds svd(Wouts[end]).S]
end

figure(246)

subplot(3,1,1)
PyPlot.plot(Tts,Tbs,"-o")
xlabel("Training time [iterations]")
ylabel("Breaktime [iterations]")
title("Lorenz system, reservoir size N = $N0")

for d in [1,2,5,10]
	subplot(3,1,2)
	PyPlot.semilogy(Tts[1:end-d],norm.(Wouts[1+d:end]-Wouts[1:end-d])./norm.(Wouts[1:end-d]),"-o",label="δT = $(100*d)")
end

subplot(3,1,2)
xlabel("Training time [iterations]")
ylabel("||W(T+δT) - W(T)||_2 / ||W(T)||_2")
legend()

subplot(3,1,3)
num = ["1st","2nd","3rd"]
for i in 1:length(svds[:,1])
	PyPlot.semilogy(Tts[1:end-1],abs.(svds[i,1:end-1] - svds[i,2:end])./abs.(svds[i,1:end-1]),"-o",label=num[i]*" singular value")
end

xlabel("Training time [iterations]")
ylabel("|μ(T+100) - μ(T)|/|μ(T)|")
legend()

Tb00 = Array{Float64,1}()
Tb25 = Array{Float64,1}()
Tb50 = Array{Float64,1}()
Tb75 = Array{Float64,1}()
Tb100 = Array{Float64,1}()

for i in 1:length(Ns)
	Tb = Array{Float64,1}()
	for k in 1:50
		push!(Tb,readdlm("data1/Tb$(k).$(i)_vs_N_Tt1001.csv",',')[1])
	end
	global Tb00,Tb25,Tb50,Tb75,Tb100
	push!(Tb00,minimum(Tb))
	push!(Tb25,quantile(Tb,.25))
	push!(Tb50,median(Tb))
	push!(Tb75,quantile(Tb,.75))
	push!(Tb100,maximum(Tb))
end

figure(369)

PyPlot.plot(Ns,Tb00,"x",color="C0")
PyPlot.plot(Ns,Tb100,"x",color="C0")
PyPlot.plot(Ns,Tb25,"--",color="C0")
PyPlot.plot(Ns,Tb75,"--",color="C0")
PyPlot.plot(Ns,Tb50,"-o",color="C0")
xlabel("Reservoir size")
ylabel("Breaktime [iterations]")
title("Lorenz system, training time Tt = 4001 [iterations]")




