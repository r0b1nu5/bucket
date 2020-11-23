using DelimitedFiles

# IEEE-300: bus_ids = (4,303), line_ids = (306,716)

function parse_cdf(file::String, name::String, bus_ids::Tuple{Int64,Int64}, line_ids::Tuple{Int64,Int64})
	xy = readdlm(file)

	x = xy[bus_ids[1]:bus_ids[2],:]
	y = xy[line_ids[1]:line_ids[2],:]

#################### BUSES #########################
	bnumber = x[:,1]
	btype = x[:,5]
	bV = x[:,6]
	ba = x[:,7]
	pl = x[:,8]
	ql = x[:,9]
	pg = x[:,10]
	qg = x[:,11]
	bvolt = x[:,12]

	P = -(pg+pl)
	Q = -(qg+ql)
	V = Float64.(bV.*bvolt)
	t = ba*pi/180

	n = length(P)

	bids =  Dict{Int64,Int64}(bnumber[i]=>i for i in 1:n)

	islack = setdiff((1:n).*(btype .== 3),[0.,])
	ipq = setdiff((1:n).*((btype .== 0)+(btype .== 1)),[0.,])
	ipv = setdiff((1:n).*(btype .== 2),[0.,])
	ids = [islack;ipq;ipv]

	
################# LINES ##########################
	lsource = y[:,1]
	ltarget = y[:,2]
	m = length(lsource)
	R = y[:,7]
	X = y[:,8]
	g = R./(R.^2 + X.^2)
	b = X./(R.^2 + X.^2)

	I = zeros(n,m)
	for i in 1:m
		I[[bids[lsource[i]],bids[ltarget[i]]],i] = [1.,-1.]
	end

	G = I*diagm(0 => g)*I'
	B = I*diagm(0 => b)*I'



############### RE-ORDER & WRITE ####################
	P = P[[ipq;ipv]]
	Q = Q[ipq]
	V = V[ipv]
	t = t[[ipq;ipv]]
	bt = Int.(btype[ids])
	G = G[ids,ids]
	B = B[ids,ids]

	writedlm("data/"*name*"_P.csb",P,',')
	writedlm("data/"*name*"_Q.csb",Q,',')
	writedlm("data/"*name*"_V.csb",V,',')
	writedlm("data/"*name*"_t.csb",t,',')
	writedlm("data/"*name*"_bt.csb",bt,',')
	writedlm("data/"*name*"_G.csb",G,',')
	writedlm("data/"*name*"_B.csb",B,',')

	return P,Q,V,t,bt,G,B
end



