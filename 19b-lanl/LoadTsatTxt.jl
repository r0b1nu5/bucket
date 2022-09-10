# The function to be called above is defined below

function LoadTsatTxt(filename::String)
	# Load .txt output file from TSAT-OA
	# Author: Bin Wang (wangbin.dianqi@gmail.com)
	# Date created: 05/28/2015

	# open the output file
	fid = open(filename,"r")
	
	# Read the first row
	nextline = readline(fid)
	title = String.(split(nextline,"' '"))
	N = length(title)

	# Read the data
	data = Array{Float64,2}(undef,N,0)
	mark(fid)
	nl = countlines(fid)
	reset(fid)
	for i in 1:nl
		nextline = split(readline(fid)," ")
		if length(nextline) == N
			temp = parse.(Float64,nextline)
		else
			temp = parse.(Float64,[nextline[1];nextline[3:end]])
		end
		data = [data temp]
	end

	close(fid)

	# Delete the duplicate data
	N,T = size(data)
	temp = Array{Int64,1}()
	for i in 2:T
		if data[1,i] == data[1,i-1]
			push!(temp,i)
		end
	end
	ids = setdiff(Array(1:T),temp)
	data = data[:,ids]

	t = data[1,:]
	title = title[2:end]
	data = data[2:end,:]

	return t, data, title
end



