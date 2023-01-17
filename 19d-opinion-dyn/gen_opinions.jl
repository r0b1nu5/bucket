# The script generates 'n' natural opinions with 'p' parties. 
# The opinions are stored as a n x p matrix in the CSV file 'file_name.csv', in the directory where the script is ran.
# The code to run the script is 
#
# "julia gen_opinions.jl p n "file_name".
#
# The code runs with Julia 1.8, requires the packages "DelimitedFiles", and the script file "simplex.jl".

using DelimitedFiles

include("simplex.jl")

p = parse(Int64,ARGS[1])
n = parse(Int64,ARGS[2])

file = ARGS[3]

as = admissible_summits(p)

X = zeros(0,p)
for i in 1:n
	global X = [X;rand_opinion(p,as)']
end

writedlm(file*".csv", X, ',')

