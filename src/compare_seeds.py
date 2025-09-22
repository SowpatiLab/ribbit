locs2 = set()
locs3 = set()

with open('/Users/avvarua/Documents/projects/genome-references/hg38/chr1_test.ribbit.seeds2', 'r') as f:
    for line in f:
        locs2.add(line.strip())

with open('/Users/avvarua/Documents/projects/genome-references/hg38/chr1_test.ribbit.seeds3', 'r') as f:
    for line in f:
        locs3.add(line.strip())

print("Number of seeds in file 2:", len(locs2))
print("Number of seeds in file 3:", len(locs3))

common = locs2 & locs3
print("Number of common seeds:", len(common))

only_in_2 = locs2 - locs3
only_in_3 = locs3 - locs2

print("Number of seeds only in file 2:", len(only_in_2))
print("Number of seeds only in file 3:", len(only_in_3))