using DelimitedFiles

dates = ["2013-01-15_00", "2013-03-10_04", "2013-04-03_02", "2013-04-03_03", "2013-04-03_07", "2013-07-30_01", "2013-07-30_04", "2013-07-30_09"]

function load_pen(date::String)
	x = readdir("data_PEN")

	for f in x
		if f[1:13] == date
			y = readdlm("data_PEN/"*f,',')[2:end,2:end]
		end
	
		N = Int(maximum(y[:,1])+1)
		n,m = size(y)
		nt = Int(n/N)	
	
		fs = Array{Float64,2}(undef,N,0)
		th = Array{Float64,2}(undef,N,0)
	
		for i in 1:nt
			fs = [fs y[(i-1)*nt+1:i*nt,6]]
			th = [th y[(i-1)*nt+1:i*nt,5]]
		end

		open("data_PEN//"*date*"_fs.csv","a") do f
			write(f,fs')
		end
	end
	


	f = zeros(N)
	t = zeros(N)
	
	i1 = 0
	i2 = 0
	for i in 1:n
		if i%1000 == 0
			@info "$i/$n"
		end

		i1 = copy(i2)
		i2 = y[i,1]+1

		if i2 > i1
			f[i2] = y[i,6]
			t[i2] = y[i,5]
		else
			fs = [fs f]
			th = [th t]
			f = zeros(N)
			t = zeros(N)
			f[i2] = y[i,6]
			t[i2] = y[i,5]
		end
	end

	writedlm("data_PEN/"*date*"_fs.csv",fs,',')
	writedlm("data_PEN/"*date*"_th.csv",th,',')

	return fs,th
end











