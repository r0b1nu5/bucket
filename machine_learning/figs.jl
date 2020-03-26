using PyPlot, LinearAlgebra, DelimitedFiles, Statistics

 #=
Tts = Array(51:50:2001)
T0 = 1001 
Ns = Array(200:200:2000)
N0 = 2000
# =#


thrs = [.1,.05,.02]

n = 3

Tbs1 = Array{Float64,1}()
Tbs2 = Array{Float64,1}()
Tbs3 = Array{Float64,1}()
Wouts = Array{Array{Float64,2},1}()
svds = Array{Float64,2}(undef,n,0)


 #= 

Tts = Int.(vec(readdlm("data1/last_run_Tts.csv",',')))
N0 = Int(readdlm("data1/last_run_N0.csv",',')[1])

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
PyPlot.plot(Tts,Tbs1,"-o",label="$(thrs[1]*100)% deviation")
PyPlot.plot(Tts,Tbs2,"-o",label="$(thrs[2]*100)% deviation")
PyPlot.plot(Tts,Tbs3,"-o",label="$(thrs[3]*100)% deviation")
xlabel("Training time [iterations]")
ylabel("Breaktime [iterations]")
title("Lorenz system, reservoir size N = $N0")
legend()

for d in [1,2,5,10]
	subplot(3,1,2)
	PyPlot.semilogy(Tts[1:end-d],norm.(Wouts[1+d:end]-Wouts[1:end-d])./norm.(Wouts[1:end-d]),"-o",label="δT = $(50*d)")
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

# =#

 ##= 

Ns = Int.(vec(readdlm("data1/last_run_Ns.csv",',')))
T0 = Int(readdlm("data1/last_run_Tt0.csv",',')[1])

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

	PyPlot.fill([Ns;Ns[length(Ns):-1:1]],[Tb00;Tb100[length(Tb100):-1:1]],color="C$(j-1)",alpha=.1)
	PyPlot.fill([Ns;Ns[length(Ns):-1:1]],[Tb25;Tb75[length(Tb75):-1:1]],color="C$(j-1)",alpha=.3)
	PyPlot.plot(Ns,Tbs50[end],"-o",color="C$(j-1)")
	
	if j == length(thrs)
		xlabel("Reservoir size")
	end
	ylabel("Breaktime [iterations]")
	title("Lorenz system, $(thr*100)% deviation, training time Tt = 4001 [iterations]")
end

# =#

 #= 

Ns = Int.(vec(readdlm("data1/last_run_Ns.csv",',')))
T0 = Int(readdlm("data1/last_run_Tt0.csv",',')[1])

Tbs00 = Array{Array{Float64,1},1}()
Tbs20 = Array{Array{Float64,1},1}()
Tbs40 = Array{Array{Float64,1},1}()
Tbs50 = Array{Array{Float64,1},1}()
Tbs60 = Array{Array{Float64,1},1}()
Tbs80 = Array{Array{Float64,1},1}()
Tbs100 = Array{Array{Float64,1},1}()

for j in 1:length(thrs)
	global Tbs00,Tbs25,Tbs50,Tbs75,Tbs100
	thr = thrs[j]
	Tb00 = Array{Float64,1}()
	Tb20 = Array{Float64,1}()
	Tb40 = Array{Float64,1}()
	Tb50 = Array{Float64,1}()
	Tb60 = Array{Float64,1}()
	Tb80 = Array{Float64,1}()
	Tb100 = Array{Float64,1}()
	for i in 1:length(Ns)
		Tb = Array{Float64,1}()
		for k in 1:50
			push!(Tb,readdlm("data1/Tb$(k).$(i)_vs_N_Tt$(T0)_thr$(thr).csv",',')[1])
		end
		push!(Tb00,quantile(Tb,0.))
		push!(Tb20,quantile(Tb,.2))
		push!(Tb40,quantile(Tb,.4))
		push!(Tb50,quantile(Tb,.5))
		push!(Tb60,quantile(Tb,.6))
		push!(Tb80,quantile(Tb,.8))
		push!(Tb100,quantile(Tb,1.))
	end
	push!(Tbs00,Tb00)
	push!(Tbs20,Tb20)
	push!(Tbs40,Tb40)
	push!(Tbs50,Tb50)
	push!(Tbs60,Tb60)
	push!(Tbs80,Tb80)
	push!(Tbs100,Tb100)

	figure(369)
	subplot(length(thrs),1,j)

	PyPlot.fill([Ns;Ns[length(Ns):-1:1]],[Tb00;Tb100[length(Tb100):-1:1]],color="C$(j-1)",alpha=.1)
	PyPlot.fill([Ns;Ns[length(Ns):-1:1]],[Tb20;Tb80[length(Tb80):-1:1]],color="C$(j-1)",alpha=.3)
	PyPlot.fill([Ns;Ns[length(Ns):-1:1]],[Tb40;Tb60[length(Tb60):-1:1]],color="C$(j-1)",alpha=.3)
	PyPlot.plot(Ns,Tbs50[end],"-o",color="C$(j-1)")
	
	if j == length(thrs)
		xlabel("Reservoir size")
	end
	ylabel("Breaktime [iterations]")
	title("Lorenz system, $(thr*100)% deviation, training time Tt = 4001 [iterations]")
end

# =#

