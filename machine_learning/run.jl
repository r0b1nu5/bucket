Tts = Array(101:10:4001)
N = 1000

for i in 1:length(Tts)
	global Tt = Tts[i]
	global T0 = 8002 - Tt
	include("nonlinear_prediction_rc.jl")
	writedlm("data1/Wout$(i)_vs_Tt.csv",Wout,',')
	writedlm("data1/Tb$(i)_vs_Tt.csv",Tb,',')
	@info "Tt = $Tt, i = $i done."
end


Tt = 4001
T0 = 8002 - Tt
Ns = [100,200,500,1000,1500,2000]

for i in 1:lenght(Ns)
	global N = Ns[i]
	include("nonlinear_prediction_rc.jl")
	writedlm("data1/Wout$(i)_vs_N.csv",Wout,',')
	writedlm("data1/Tb$(i)_vs_N.csv",Tb,',')
	@info "N = $N, i = $i done."
end




