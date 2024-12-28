from collections import defaultdict
from Levenshtein import distance
import pybedtools
import os
import argparse
from tqdm import tqdm

def expected_diff(parent_motif, parent_purity, child_motif, child_len, child_purity):
    if len(parent_motif) < len(child_motif):
        least_d = float('inf')
        least_m = 0
        for i in range(1, (len(parent_motif) // len(child_motif)) + 2):
            d = distance(parent_motif * i, child_motif)
            if d < least_d:
                least_d = d
                least_m = i
        parent_motif = parent_motif * least_m
    elif len(child_motif) < len(parent_motif):
        least_d = float('inf')
        least_m = 0
        for i in range(1, (len(child_motif) // len(parent_motif)) + 2):
            d = distance(child_motif * i, parent_motif)
            if d < least_d:
                least_d = d
                least_m = i
        child_motif = child_motif * least_m

    d = float('inf')
    for i in range(len(parent_motif)):
        pm = parent_motif[i:] + parent_motif[:i]
        for j in range(len(child_motif)):
            cm = child_motif[j:] + child_motif[:j]
            if distance(pm, cm) < d:
                d = distance(pm, cm)

    parent_imp = int((1 - parent_purity) * child_len)
    total_d = int(child_len // len(child_motif) * d)

    return not (child_purity > (1 - (abs(parent_imp - total_d) / child_len)))

def distance(s1, s2):
    return sum(1 for a, b in zip(s1, s2) if a != b)

def process_bed(input_bed_file, output_gff_file):
    bed = pybedtools.BedTool(input_bed_file)
    sorted_merged_bed = bed.sort().merge(c=[2,3,4,5,6,7,8], o=['collapse', 'collapse', 'collapse', 'collapse', 'collapse', 'collapse', 'collapse'])

    out = open(output_gff_file, 'w')
    print('#gff-version 3', file=out)

    region_num = 0
    str_num = 0

    for region in tqdm(sorted_merged_bed):
        region_num += 1

        chrom = region.chrom
        start = region.start
        end = region.end
        str_starts = [int(x) for x in region[3].split(',')]
        str_ends = [int(x) for x in region[4].split(',')]
        str_motifs = region[5].split(',')
        str_purities = [float(x) for x in region[6].split(',')]
        str_ids = ["S%09d" % (str_num+(i+1)) for i in range(len(str_starts))]
        str_lens = [str_ends[i] - str_starts[i] for i in range(len(str_starts))]
        seed_type = region[8].split(',')
        str_cigars = region[9].split(',')
        region_id = "R%09d" %(region_num)

        del_strs = []
        nested_relations = {}
        overlap_relations = defaultdict(set)

        for i in range(len(str_starts)):
            a = str_starts[i]
            b = str_ends[i]
            m = str_motifs[i]
            p = str_purities[i]

            for j in range(len(str_starts)):
                if i == j:
                    continue

                c = str_starts[j]
                d = str_ends[j]

                if a >= c and b <= d: # existing location is nested within the new location.
                    if p <= str_purities[j] or str_motifs[i] == str_motifs[j] or expected_diff(str_motifs[j], str_purities[j], str_motifs[i], str_lens[i], str_purities[i]):
                        del_strs.append(i) #lower purity, similar motif size, 
                    else:
                        nested_relations[i] = j
                else:
                    if c <= a and a <= d:
                        if (str_lens[i] < str_lens[j]) and ((b - d) < len(m) or ((b - d) < 3)):
                            del_strs.append(i)
                        elif (str_lens[i] == str_lens[j]) and str_purities[i] < str_purities[j]:
                            del_strs.append(i)
                        else:
                            overlap_relations[i].add(j)
                            overlap_relations[j].add(i)
                    elif c <= b and b <= d:
                        if (str_lens[i] < str_lens[j]) and (((c - a) < len(m)) or ((c - a) < 3)):
                            del_strs.append(i)
                        elif (str_lens[i] == str_lens[j]) and str_purities[i] < str_purities[j]:
                            del_strs.append(i)
                        else:
                            overlap_relations[i].add(j)
                            overlap_relations[j].add(i)

                if b < str_starts[j]:
                    break

        if len(str_starts) - len(set(del_strs)) > 1:
            print(chrom, 'ribbit', 'REGION', start, end, '.', '+', '.', f'id={region_id}', sep='\t', file=out)

            for i in range(len(str_starts)):
                if i in del_strs:
                    continue
                str_id = str_ids[i]
                parent = ''
                children = []
                overlaps = []
                if i in nested_relations:
                    parent = str_ids[nested_relations[i]]
                else:
                    parent = region_id

                for n in nested_relations:
                    if nested_relations[n] == i and n not in del_strs:
                        children.append(str_ids[n])
                children = ",".join(children)

                for n in overlap_relations[i]:
                    if n not in del_strs:
                        overlaps.append(str_ids[n])
                overlaps = ",".join(overlaps)

                print(chrom, 'ribbit', 'RSTR', str_starts[i], str_ends[i], '.', '+', '.',
                        f'id={str_id};purity={str_purities[i]};motif={str_motifs[i]};cigar={str_cigars[i]};seed_type={seed_type[i]};parent={parent};children={children};overlaps={overlaps}',
                        sep='\t', file=out)

        elif len(str_starts) - len(set(del_strs)) == 1:
            for i in range(len(str_starts)):
                str_id = str_ids[i]
                if i in del_strs:
                    continue
                print(chrom, 'ribbit', 'STR', str_starts[i], str_ends[i], '.', '+', '.',
                        f'id={str_id};purity={str_purities[i]};motif={str_motifs[i]};cigar={str_cigars[i]};seed_type={seed_type[i]}',
                        sep='\t', file=out)

        str_num += len(str_starts)

    out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a single BED file and generate a corresponding GFF file")
    parser.add_argument("input_bed_file", help="Input BED file")
    parser.add_argument("output_gff_file", help="Output GFF file")

    args = parser.parse_args()

    input_bed_file = args.input_bed_file
    output_gff_file = args.output_gff_file

    process_bed(input_bed_file, output_gff_file)
