using DelimitedFiles

rs = [.0,.1,.2,.3,.4,.5]

for r in rs
	P = readdlm("data-robin/P_$r.csv",',')
	Q = readdlm("data-robin/Q_$r.csv",',')
	V = readdlm("data-robin/V_$r.csv",',')

	n,m = size(P)

	dat = ["Time" "Site" "P" "Q" "V"]

	for i in 2:56
		dat = [dat;[P[2:end,1] repeat(P[[1,],i],n-1,1) P[2:end,i] Q[2:end,i] V[2:end,i]]]
	end

	writedlm("data-robin/data-$r.csv",dat,',')
end




