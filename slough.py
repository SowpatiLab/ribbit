#! /usr/bin/env python


import pybedtools
from collections import defaultdict


chroms = ['F4BF9E6E99_1', 'F4BF9E6E99_2', 'F4BF9E6E99_3', 'F4BF9E6E99_4', 'F4BF9E6E99_5', 'F4BF9E6E99_6', 'F4BF9E6E99_7']
ends   = [250002018, 250000488, 250001580, 250001566, 250002582, 250000239, 249994337]
# ends   = [13500, 250000488, 250001580, 250001566, 250002582, 250000239, 249994337]

chrom_ends = { chroms[i]:ends[i] for i in range(len(chroms))}
operations = set(['=', 'X', 'I', 'D'])

window = 50000
overlap = 5000

bedfile = pybedtools.BedTool('./test.bed')
columns = list(range(2, 12))
conditions = ['collapse']*len(columns)


def calculate_purity(cigar):
    match_len = 0; align_len = 0; clen = ''
    for x in cigar:
        if x in operations:
            clen = int(clen)
            if  x == '=': match_len += clen
            align_len += clen
            clen = ''
        else: clen += x
    return match_len/align_len


def merge_overlapping(a, b, c, d, str_cigars, i, j):
    
    pos = c; clen = ''
    for c, x in enumerate(str_cigars[j]):
        if x in operations:
            clen = int(clen)
            if  x == '=' or x == 'I' or x == 'X':
                pos += clen
                if pos >= b:
                    clen = clen - (b - pos)
                    cigar = str_cigars[i] + str(clen) + str_cigars[j][c:]
                    return [calculate_purity(cigar), cigar]
            clen = ''
        else: clen += x




region_num = 0
str_num = 0
for chrom in chroms[:1]:
    for range_start in range(0, chrom_ends[chrom], window-overlap):
        range_end   = range_start + window
        interval = pybedtools.BedTool([pybedtools.Interval(chrom, range_start, range_end)])
        range_bed = bedfile.intersect(interval)
        
        if range_bed is None or (not range_bed): continue
        
        range_bed = range_bed.sort().merge(c=columns, o=conditions)
        for region in range_bed:
            region_num += 1

            chrom = region[0]
            start = region[1]
            end = region[2]
            str_starts = [int(x) for x in region[3].split(',')]
            str_ends   = [int(x) for x in region[4].split(',')]
            str_motifs   = region[5].split(',')
            str_mlens = [int(x) for x in region[6].split(',')]
            str_rlens = [int(x) for x in region[7].split(',')]
            str_units = [int(x) for x in region[8].split(',')]
            str_purities = [float(x) for x in region[9].split(',')]
            str_oris = region[10].split(',')
            str_types  = [int(x.split('-')[1]) for x in region[11].split(',')]
            str_cigars = region[12].split(',')
            
            str_ids  = ["S%09d" % (str_num+(i+1)) for i in range(len(str_starts))]
            region_id  = "R%09d" %(region_num)

            nested_relations = {}
            overlap_relations = defaultdict(set)
            global_del_strs = []

            while True:

                del_strs = []

                for i in range(len(str_starts)):
                    a = str_starts[i]
                    b = str_ends[i]
                    m = str_motifs[i]
                    p = str_purities[i]

                    overlappers = []
                    if i in global_del_strs or i in del_strs: continue

                    for j in range(len(str_starts)):
                        # do not compare with self
                        if i == j or j in global_del_strs or j in del_strs: continue

                        c = str_starts[j]
                        d = str_ends[j]
                        
                        if b < c: break

                        if c == a and b == d:
                            if p < str_purities[j] or m == str_motifs[j]:
                                del_strs.append(i); global_del_strs.append(i)
                            elif p > str_purities[j]:
                                del_strs.append(j); global_del_strs.append(j)

                        elif c <= a and b <= d: # existing location is nested within the new location.
                            if p <= str_purities[j] or m == str_motifs[j]:
                                del_strs.append(i); global_del_strs.append(i)
                        
                        elif a <= c and d <= b:
                            if str_purities[j] < p or m == str_motifs[j]:
                                del_strs.append(j); global_del_strs.append(j)
                            if str_purities[j] >= p: overlappers.append(j)

                        else:
                            if str_motifs[i] == str_motifs[j]:
                                if i < j:
                                    str_purities[i], str_cigars[i] = merge_overlapping(a, b, c, d, str_cigars, i, j)
                                else:
                                    str_purities[i], str_cigars[i] = merge_overlapping(c, d, a, b, str_cigars, j, i)
                                if str_starts[j] < str_starts[i]: str_starts[i] = 0 + str_starts[j]
                                if str_ends[i] < str_ends[j]: str_ends[i] = 0 + str_ends[j]
                                str_rlens[i] = str_ends[i] - str_starts[i]
                                str_units[i] = str_rlens[i] // str_units[i]
                                del_strs.append(j); global_del_strs.append(j)
                                continue
                            if str_purities[j] >= p:  overlappers.append(j)
                    
                    if len(overlappers) > 0:
                        coverage = 0; prev = str_starts[overlappers[0]]
                        for o in overlappers:
                            if str_starts[o] <= prev <= str_ends[o]:
                                coverage += str_ends[o] - prev
                            elif prev < str_starts[o]:
                                coverage += str_ends[o] - str_starts[o]
                            prev = str_ends[o]
                        
                        if str_rlens[i] - coverage <= 2:
                            # print(i, overlappers)
                            del_strs.append(i); global_del_strs.append(i)

                if len(del_strs) == 0: break
            
            for i in range(len(str_starts)):
                if i not in global_del_strs:
                    print(chrom, str_starts[i], str_ends[i], str_motifs[i], str_mlens[i], str_rlens[i], str_units[i],
                          str_purities[i], str_oris[i], str_types[i], str_cigars[i], sep='\t')
            str_num += len(str_starts)