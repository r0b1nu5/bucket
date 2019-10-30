@everywhere using DelimitedFiles, Distributed

@everywhere dates = ["2013-01-15_00", "2013-03-10_04", "2013-04-03_02", "2013-04-03_03", "2013-04-03_07", "2013-07-30_01", "2013-07-30_04", "2013-07-30_09"]

@everywhere function load_pen(date::String)
	x = readdir("data_PEN")

	for f in x
		if f[1:13] == date
			@info "Start: "*f

			y = Float64.(readdlm("data_PEN/"*f,',')[2:end,2:end])
	
			N = Int(maximum(y[:,1])+1)
			n,m = size(y)
			did = [Int64.(setdiff([1.;Float64.((y[2:end,1] - y[1:end-1,1]) .< 0)] .* (1:n),[0,]));n+1]
			fs = Array{Float64,2}(undef,N,0)
			th = Array{Float64,2}(undef,N,0)
		
			for i in 1:length(did)-1
				F = zeros(N)
				T = zeros(N)
				F[Int.(y[did[i]:did[i+1]-1,1]) .+ 1] = y[did[i]:did[i+1]-1,6]
				T[Int.(y[did[i]:did[i+1]-1,1]) .+ 1] = y[did[i]:did[i+1]-1,5]
				fs = [fs F]
				th = [th T]
			end
	
			open("data_PEN/"*date*"_fs.csv","a") do f
				writedlm(f,fs',',')
			end
			open("data_PEN/"*date*"_th.csv","a") do f
				writedlm(f,th',',')
			end

			@info "Done: "*f
		end
	end
end

pmap(load_pen,dates)


