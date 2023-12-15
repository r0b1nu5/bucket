using PyPlot, DelimitedFiles, Statistics

include("eeg-tools.jl")

n = 7
type = "avg" # "avg" or "sub"
data_exist = true
Tmax = 1000
thr2 = 1.
thr3 = 10.

# Data to be loaded
subjects = list_all_subjects(109)
#subjects = ["001","002"]
#states = ["01","02"]
states = ["03","07","11"]

suffix = "xx7"

# Sensor to zone pairing
s = readdlm("eeg-data/sensors-$n.csv",',',String)
z = readdlm("eeg-data/zones-$n.csv",',',Int64)
s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

#if data_exist
#	ρ3 = readdlm("eeg-data/violin-data-"*type*"-$n.csv",',')
#else
	ρ3 = zeros(n,0)
#end

AA2 = zeros(n,n)
AAA2 = zeros(n,n,length(subjects)*length(states))
AA3 = zeros(n,n,n)
AAA3 = zeros(n,n,n,length(subjects)*length(states))
vA2 = Float64[]
vA3 = Float64[]
contr_per_subject = zeros(n,0)

c = 0
for su in subjects
	for st in states
		global c += 1
		@info "Working on S"*su*"R"*st*": run "*suffix
		if !(data_exist)
			s2signal = read_eeg("eeg-data/S"*su*"R"*st*".edf")
			asig = average_over_zones(s2signal,s2z)
			dt = 1/160
			truncat = findmin(vec(maximum(abs.([asig zeros(n)]),dims=1)) .> 1e-6)[2] - 1
#			truncat = 3001
			X0 = denoise_fourier(asig[:,1:truncat],100)
			X0 = X0[:,1:end-1]
			X = X0./mean(abs.(X0))
		
			X,ids = restrict_box_size(X,1000)
			T = size(X)[2]
		end

		A2 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A2-"*suffix*".csv",',')
		a3 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A3-"*suffix*".csv",',')
		A3 = zeros(n,n,n)
		for k in 1:n
			A3[:,:,k] = a3[:,(k-1)*n .+ (1:n)]
		end

		global vA2 = [vA2;vectorize_adj(A2)]
		global vA3 = [vA3;vectorize_adj(A3)]

		if !(data_exist)
			writedlm("eeg-data/T-S"*su*"R"*st*"-"*suffix*".csv",T,',')
			for t in 1:T
				ρ = Float64[]
				for i in 1:n
					z2 = 0.
					z3 = 0.
					for j in 1:n
						for k in 1:n
							z3 += abs(A3[i,j,k]*X[j,t]*X[k,t])
						end
						z2 += abs(A2[i,j]*X[j,t])
					end
					if z2+z3 == 0.
						push!(ρ,0.)
					else
						push!(ρ,z3/(z2+z3))
					end
				end
#				global ρ3 = [ρ3 ρ]
				writedlm("eeg-data/violin-data-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t.csv",ρ,',')
			end
			ρ = zeros(n,0)
			for t in 1:T
				ρ = [ρ vec(readdlm("eeg-data/violin-data-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t.csv",','))]
				rm("eeg-data/violin-data-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t.csv")
				if t%1000 == 0
					writedlm("eeg-data/violin-data-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T$t.csv",ρ,',')
					ρ = zeros(n,0)
				end
			end

		end

		global AA2 += (abs.(A2) .> thr2)
		global AAA2[:,:,c] = abs.(A2)
		global AA3 += (abs.(A3) .> thr3)
		global AAA3[:,:,:,c] = abs.(A3)
	end
end

MA2 = zeros(n,n)
MA3 = zeros(n,n,n)
for i in 1:n
	for j in 1:n
		MA2[i,j] = median(AAA2[i,j,:])
		for k in 1:n
			MA3[i,j,k] = median(AAA3[i,j,k,:])
		end
	end
end

q0 = zeros(n,0)
q1 = zeros(n,0)
q2 = zeros(n,0)
q3 = zeros(n,0)
q4 = zeros(n,0)
nams = String[]
for su in subjects
	for st in states
		@info "Loading S"*su*"R"*st

		T = Int64(1000*floor(readdlm("eeg-data/T-S"*su*"R"*st*"-"*suffix*".csv",',')[1]/1000))
		ρ = zeros(n,0)
		for t in 1000:1000:min(T,Tmax)
			ρ = [ρ readdlm("eeg-data/violin-data-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T$t.csv",',')]
		end
		global ρ3 = [ρ3 ρ]
		global contr_per_subject = [contr_per_subject mean(ρ,dims=2)]
		global q0 = [q0 [quantile(ρ[i,:],.00) for i in 1:n]]
		global q1 = [q1 [quantile(ρ[i,:],.25) for i in 1:n]]
		global q2 = [q2 [quantile(ρ[i,:],.50) for i in 1:n]]
		global q3 = [q3 [quantile(ρ[i,:],.75) for i in 1:n]]
		global q4 = [q4 [quantile(ρ[i,:],1.0) for i in 1:n]]
		push!(nams,"S"*su*"R"*st)
#		figure("Violins-bis",(15,6))
#		fig = plt.violinplot(ρ',showextrema=false)
#		for pc in fig["bodies"]
#			pc.set_alpha(0.05)
#		end
	end
end

if !(data_exist)
	writedlm("eeg-data/violin-data-"*type*"-$n-"*suffix*".csv",ρ3,',')
end


# Fig 1
#figure("Violins-bis")
#plt.violinplot(ρ3')

figure("Violins-ter",(10,4))
ζ = [vec(ρ3),]
#fig = plt.violinplot(vec(ρ3),positions=[1,],showextrema=false)
qs = quantile(vec(q2),[0.,.25,.5,.75,1.])
#qs = quantile(vec(q2),[.01,.25,.5,.75,.99])
#qs = quantile(vec(q2),[.05,.25,.5,.75,.95])
#qs = quantile(vec(q2),[.2,.4,.5,.8])
m0,i0 = findmin(abs.(q2 .- qs[1]))
m1,i1 = findmin(abs.(q2 .- qs[2]))
m2,i2 = findmin(abs.(q2 .- qs[3]))
m3,i3 = findmin(abs.(q2 .- qs[4]))
m4,i4 = findmin(abs.(q2 .- qs[5]))
for i in [i0,i1,i2,i3,i4] 
	push!(ζ,ρ3[i[1],(i[2]-1)*Tmax .+ (1:Tmax)])
end
#ζ = [ρ3[i0[1],(i0[2]-1)*Tmax .+ (1:Tmax)] ρ3[i1[1],(i1[2]-1)*Tmax .+ (1:Tmax)] ρ3[i2[1],(i2[2]-1)*Tmax .+ (1:Tmax)] ρ3[i3[1],(i3[2]-1)*Tmax .+ (1:Tmax)] ρ3[i4[1],(i4[2]-1)*Tmax .+ (1:Tmax)]]
fig = plt.violinplot(ζ,showextrema=false)
cmap = get_cmap("magma")
cols = [cmap(r) for r in [.1,.5,.6,.7,.8,.9]]
c = 0
for pc in fig["bodies"]
	global c += 1
	pc.set_facecolor(cols[c])
	pc.set_edgecolor("black")
	pc.set_alpha(.8)
end
cc = ["gray","black","black","black","black","black"]
for i in 1:length(ζ)
	PyPlot.plot([i,i],[minimum(ζ[i]),maximum(ζ[i])],color=cc[i])
	PyPlot.plot([i,i],quantile(ζ[i],[.25,.75]),color=cc[i],linewidth=5)
	PyPlot.plot(i,median(ζ[i]),"o",color=cc[i],markersize=10)
end
#xticks([1,2,3,4,5,6],["All","1%","25%","50%","75%","99%"])
xticks([1,2,3,4,5,6],["All","0%","25%","50%","75%","100%"])
xlabel("Percentile")
ylabel("Amount of the dynamics that is\nexplained by 3rd-order interactions")

figure("Violins",(15,6))
#for i in 1:n
#	PyPlot.plot(i*ones(length(ρ3[i,:])),ρ3[i,:],"xk",alpha=.002)
#end
fig = plt.violinplot(ρ3',showextrema=false,showmedians=false)
cols = [cmap(r) for r in LinRange(.2,.8,7)]
c = 0
for pc in fig["bodies"]
	global c += 1
#	pc.set_facecolor("#d41367")
	pc.set_facecolor(cols[c])
	pc.set_edgecolor("black")
	pc.set_alpha(.5)
end
for i in 1:n
	PyPlot.plot([i,i],[minimum(ρ3[i,:]),maximum(ρ3[i,:])],"k")
	PyPlot.plot([i,i],quantile(ρ3[i,:],[.25,.75]),"k",linewidth=5)
	PyPlot.plot(i,median(ρ3[i,:]),"ok",markersize=10)
end

xlabel("Area")
ylabel("Amount of the dynamics that is explained by 3rd-order interactions")

# Fig 2
N = length(subjects)*length(states)
M2 = maximum(AA2)
M3 = maximum(AA3)

figure("Histograms - average - $n",(15,6))
subplot(2,1,1)
vAA2 = vectorize_adj(AA2)
#PyPlot.hist(100*vec(AA2)./N,bins=100*((0:.5:M2+.5) .- .25)./N,color="#d41367")
PyPlot.hist(100*vAA2./N,bins=100*((0:1:M2+.5) .- .25)./N,color=cols[1])
PyPlot.plot([-25/N,(100*M2+25)/N],[n^2,n^2],"--k")
dM = 10*(M2/N > .5) + 5*(.05 < M2/N < .5) + 1*(M2/N < .05)
xticks(0:dM:(100*M2/N + dM - 1e-4))
#xticks(100*(0:ceil(M2/10):M2)./N)
ylabel("# of 2-edges")
subplot(2,1,2)
vAA3 = vectorize_adj(AA3)
#PyPlot.hist(100*vec(AA3)./N,bins=100*((0:.5:M3+.5) .- .25)./N,color="#d41367")
PyPlot.hist(100*vAA3./N,bins=100*((0:1:M3+.5) .- .25)./N,color=cols[1])
PyPlot.plot([-25/N,(100*M3+25)/N],[n^3,n^3],"--k")
dM = 10*(M3/N > .5) + 5*(.05 < M3/N < .5) + 1*(M3/N < .05)
xticks(0:dM:(100*M2/N + dM - 1e-4))
#xticks(100*(0:ceil(M3/10):M3)./N)
xlabel("Appears in x% of the inferred hypergraphs")
ylabel("# of 3-edges")

# Fig 3
figure("Brains, thr = $thr3",(15,8))
nr = 3
nc = 4
c = 0
nrc = nr*nc

a3temp = copy(AA3)
for i in 1:nr
	for j in 1:nc
		global c += 1
		v,idx = findmax(a3temp)
		ids = [idx[k] for k in 1:3]
		a3temp[idx] = -1000.
		subplot(nr,nc,c)
		plot_brain_3dge(ids,v/N,cols)
	end
end







