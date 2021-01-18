using Distributed

@everywhere using DelimitedFiles

@everywhere dates = ["2013-01-15_00", "2013-03-10_04", "2013-04-03_02", "2013-04-03_03", "2013-04-03_07", "2013-07-30_01", "2013-07-30_04", "2013-07-30_09"]


# Loads the data of Pennsylvania as a 2nxT matrix, ready to use in our algorithm.
# Missing data are interpolated.
# Angles are converted from degrees to radians.
# Last oscillator is considered as the reference node.

function load_pen(date::String)
	th = readdlm("data_pen/th_"*date*".csv",',')' * pi/180
	fs = readdlm("data_pen/fs_"*date*".csv",',')'

	N,T = size(fs)

# Remove the rows without measurements
	zt = setdiff((sum(th .!= 0,dims=2) .== 0).*(1:N),[0,])
	zf = setdiff((sum(fs .!= 0,dims=2) .== 0).*(1:N),[0,])

	z1 = sort(union(zf,zt))

# Remove the rows that are copies of other rows
	cc = Array{Int64,1}()
	for i in 1:N
		c = 1
		te = false
		while (i-c) > 0 && !te && length(intersect(z1,i)) == 0
			te = (th[i,:] == th[i-c,:]) || (fs[i,:] == fs[i-c,:])
			c += 1
		end
		if te
			push!(cc,i)
		end
	end
	
	z2 = sort(union(z1,cc))
	zc = setdiff((1:N),z2)
	
	writedlm("data_pen/pen_"*date*"_ids.csv",zc,',')

	tth = th[:,1]
	dth = th[:,2:end] - th[:,1:end-1]
	gp = dth .> pi
	lp = dth .< -pi
	tth = [th[:,1] th[:,2:end]+2pi*(lp - gp)]

	ffs = (fs .- 60) * pi/180
	
	# I cannot explain why we need this 2.5 factor, but it is needed to have the same measured and empirical derivative, i.e., to have ω_i(t) \approx (θ_i(t)-θ_i(t-1))*dt
	tth *= 2.5

	Xs = [tth[zc,:]; ffs[zc,:]]

	return Xs
end




# Extracts phases and frequencies from the raw data of Pennsylvania.

@everywhere function treat_pen(date::String)
	x = readdir("data_pen")

	for f in x
		if f[1:13] == date
			@info "Start: "*f

			y = Float64.(readdlm("data_pen/"*f,',')[2:end,2:end])
	
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
	
			open("data_pen/fs_"*date*".csv","a") do f
				writedlm(f,fs',',')
			end
			open("data_pen/th_"*date*".csv","a") do f
				writedlm(f,th',',')
			end

			@info "Done: "*f
		end
	end
end


