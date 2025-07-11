from partition import get_partition  # Assumes get_partition is defined in partition.py

list_meter = ["e_1", "e_2", "e_3", "e_4", "e_5", "e_6", "e_7", "e_8"]
list_pe = [
    ["e_2", "e_1"],
    ["e_7", "e_5"],
    ["e_8", "e_7"]
]
list_fs = [
    ["e_3", "e_4"],
    ["e_6", "e_5"],
    ["e_2", "e_7"]
]

part = get_partition(list_meter, list_pe, list_fs)

print("=" * 30)
for i, S in enumerate(part, start=1):
    print(f"S_{i}:")
    for e in S:
        print(e)
    print("=" * 30)
