using PyPlot, LinearAlgebra, DelimitedFiles, Statistics

Tts = Array(51:50:2001)
T0 = 1001 
Ns = Array(200:200:2000)
N0 = 2000
thrs = [.1,.05,.02]

n = 3

Tbs1 = Array{Float64,1}()
Tbs2 = Array{Float64,1}()
Tbs3 = Array{Float64,1}()
Wouts = Array{Array{Float64,2},1}()
svds = Array{Float64,2}(undef,n,0)

for i in 1:length(Tts)
	global Tbs,Wouts,svds
	push!(Tbs1,readdlm("data1/Tb$(i)_vs_Tt_N$(N0)_thr$(thrs[1]).csv",',')[1])
	push!(Tbs2,readdlm("data1/Tb$(i)_vs_Tt_N$(N0)_thr$(thrs[2]).csv",',')[1])
	push!(Tbs3,readdlm("data1/Tb$(i)_vs_Tt_N$(N0)_thr$(thrs[3]).csv",',')[1])
	push!(Wouts,readdlm("data1/Wout$(i)_vs_Tt_N$(N0).csv",','))
	svds = [svds svd(Wouts[end]).S]
end

figure(246)

subplot(3,1,1)
PyPlot.plot(Tts,Tbs5,"-o",label="$(thrs[1]*100)% deviation")
PyPlot.plot(Tts,Tbs2,"-o",label="$(thrs[2]*100)% deviation")
PyPlot.plot(Tts,Tbs1,"-o",label="$(thrs[3]*100)% deviation")
xlabel("Training time [iterations]")
ylabel("Breaktime [iterations]")
title("Lorenz system, reservoir size N = $N0")
legend()

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

Tbs00 = Array{Array{Float64,1},1}()
Tbs25 = Array{Array{Float64,1},1}()
Tbs50 = Array{Array{Float64,1},1}()
Tbs75 = Array{Array{Float64,1},1}()
Tbs100 = Array{Array{Float64,1},1}()

for j in 1:length(thrs)
	global Tbs00,Tbs25,Tbs50,Tbs75,Tbs100
	thr = thrs[j]
	Tb00 = Array{Float64,1}()
	Tb25 = Array{Float64,1}()
	Tb50 = Array{Float64,1}()
	Tb75 = Array{Float64,1}()
	Tb100 = Array{Float64,1}()
	for i in 1:length(Ns)
		Tb = Array{Float64,1}()
		for k in 1:50
			push!(Tb,readdlm("data1/Tb$(k).$(i)_vs_N_Tt$(T0)_thr$(thr).csv",',')[1])
		end
		push!(Tb00,minimum(Tb))
		push!(Tb25,quantile(Tb,.25))
		push!(Tb50,median(Tb))
		push!(Tb75,quantile(Tb,.75))
		push!(Tb100,maximum(Tb))
	end
	push!(Tbs00,Tb00)
	push!(Tbs25,Tb25)
	push!(Tbs50,Tb50)
	push!(Tbs75,Tb75)
	push!(Tbs100,Tb100)

	figure(369)
	subplot(length(thrs),1,j)

	PyPlot.plot(Ns,Tbs00[end],"x",color="C$(j-1)")
	PyPlot.plot(Ns,Tbs100[end],"x",color="C$(j-1)")
	PyPlot.plot(Ns,Tbs25[end],"--",color="C$(j-1)")
	PyPlot.plot(Ns,Tbs75[end],"--",color="C$(j-1)")
	PyPlot.plot(Ns,Tbs50[end],"-o",color="C$(j-1)")
	
	if j == length(thrs)
		xlabel("Reservoir size")
	end
	ylabel("Breaktime [iterations]")
	title("Lorenz system, $(thr*100)% deviation, training time Tt = 4001 [iterations]")
end





