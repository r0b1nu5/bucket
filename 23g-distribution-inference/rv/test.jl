include("partition.jl")

list_meter = ["e_1","e_2","e_3","e_4","e_5","e_6","e_7","e_8"]
list_pe = [["e_2","e_1"],
	   ["e_7","e_5"],
	   ["e_8","e_7"]]
list_fs = [["e_3","e_4"],
	   ["e_6","e_5"],
	   ["e_2","e_7"]]

part = get_partition(list_meter,list_pe,list_fs)

c = 0
for S in part
	global c += 1
	@info "==========================="
	@info "S_$c:"
	for e in S
		@info e
	end
end
@info "==========================="




