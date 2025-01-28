from collections import defaultdict
from Levenshtein import distance
import pybedtools
import os
import argparse
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(description="Process ribbit's output BED file to merge overlapping repeats and generate GFF file")
    parser.add_argument("-bed", type=str, required=True, help="Input BED file")
    parser.add_argument("-gff", type=str, help="Output GFF file")

    args = parser.parse_args()

    if args.gff == "":
        args.gff = '.'.join((args.bed).split('.')[:-1]) + '.gff'

    return args


def expected_diff(parent_motif, parent_purity, child_motif, child_len, child_purity):
    
    if len(parent_motif) < len(child_motif):    # motif length of parent STR is shorter than motif length of nested STR
        least_d = float('inf')      # least edit distance
        parent_units = 0            # number of parent units

        # identifying length of perfect parent STR that has least edit distance with one motif of nested STR
        for i in range(1, (len(parent_motif) // len(child_motif)) + 2):
            edit_d = distance(parent_motif * i, child_motif)
            if edit_d < least_d:
                least_d = edit_d
                parent_units = i
        parent_motif = parent_motif * parent_units
    
    elif len(child_motif) < len(parent_motif):  # motif length of the nested STR is shorter than motif length of parent STR
        least_d = float('inf')      # least edit distance
        nested_units = 0            # number of nested units

        # identifying the length of nested STR that has least edit distance with one motif of parent STR
        for i in range(1, (len(child_motif) // len(parent_motif)) + 2):
            edit_d = distance(child_motif * i, parent_motif)
            if edit_d < least_d:
                least_d = edit_d
                nested_units = i
        child_motif = child_motif * nested_units

    # identifying the least edit distance between different cyclical variations of motifs of parent and nested STRs
    edit_d = float('inf')
    for i in range(len(parent_motif)):
        pm = parent_motif[i:] + parent_motif[:i]
        for j in range(len(child_motif)):
            cm = child_motif[j:] + child_motif[:j]
            if distance(pm, cm) < edit_d:
                edit_d = distance(pm, cm)

    parent_imperfections = int((1 - parent_purity) * child_len)       # number of imperfections in the parent STR based on the length of nested STR
    total_edit_d = int(child_len // len(child_motif) * edit_d)   # possible total edit distance between parent STR and nested STR for length of nested STR

    # if the purity of nested STR is greater than the edit distance between parent and nested STR
    # retain the nested STR
    return not (child_purity > (1 - (abs(parent_imperfections - total_edit_d) / child_len)))


def distance(s1, s2):
    return sum(1 for str_start_i, str_end_i in zip(s1, s2) if str_start_i != str_end_i)


def cigar_split(cigar):
    length = "";
    clens = []; ctypes = []
    
    for c in cigar:
        if c.isdigit(): length += c
        else:
            clens.append(int(length))
            ctypes.append(c)
            length = ""

    return clens, ctypes


def get_merge_coordinate(iend, jstart, clens, ctypes):
    rpos = jstart
    for i in range(len(clens)):
        clen = clens[i]; ctype = ctypes[i]
        if ctype == '=' or ctype == 'X' or ctype == 'I':
            rpos += clen
        if rpos == iend:
            return clens[i+1:], ctypes[i+1:]
        elif rpos > iend:
            if i < len(clens) - 1: [rpos-iend] + clens[i+1:], ctypes[i:]
            elif i == len(clens) - 1: return [rpos-iend], ctypes[i:]
    if rpos >= iend:
        if rpos == iend: return [], []
        elif rpos > iend: return [rpos-iend], [ctypes[i]]


def merge_repeats(iend, jstart, icigar, jcigar):
    
    iclens, ictypes = cigar_split(icigar)    
    jclens, jctypes = cigar_split(jcigar)    
    fclens, fctypes = get_merge_coordinate(iend, jstart, jclens, jctypes)
    
    clens = iclens + fclens
    ctypes = ictypes + fctypes
    merge_cigar = ''.join([f'{clens[i]}{ctypes[i]}' for i in range(len(clens))])

    alignment_length = sum(clens)
    matches = 0; mismatches = 0
    for i in range(len(ctypes)):
        clen = clens[i]; ctype = ctypes[i]
        if ctype == '=': matches += clen
        else: mismatches += clen
    
    merge_purity = matches/alignment_length
    
    return [merge_cigar, merge_purity]


def process_bed(bed_file, gff_file):
    bed = pybedtools.BedTool(bed_file)
    sorted_merged_bed = bed.sort().merge(c=list(range(2,12)), o=['collapse']*10)

    out = open(gff_file, 'w')

    # the chosen gff format is gff3
    print('#gff-version 3', file=out)

    region_num = 0
    str_num = 0

    for region in sorted_merged_bed:
        region_num += 1     # counting the number of regions

        region_fields = str(region).split('\t')
        chrom = region_fields[0]
        start = int(region[1])
        end = int(region[2])

        # splitting the attributes for individual regions in str_start_i merged region
        str_starts = [int(x) for x in region[3].split(',')]
        str_ends = [int(x) for x in region[4].split(',')]
        str_motifs = region[5].split(',')
        str_motif_lengths = region[6].split(',')
        str_lengths = [int(x) for x in region[7].split(',')]
        str_units = region[8].split(',')
        str_purities = [float(x) for x in region[9].split(',')]
        seed_orientations = region[10].split(',')
        seed_type = region[11].split(',')
        str_cigars = region[12].split(',')

        region_id = "R%09d" %(region_num)
        str_ids = ["S%09d" % (str_num+(i+1)) for i in range(len(str_starts))]

        del_strs = []
        nested_relations = {}
        overlap_relations = defaultdict(set)

        nmerge = len(str_starts)
        for i in range(nmerge):

            str_start_i = str_starts[i]
            str_end_i = str_ends[i]
            str_motif_i = str_motifs[i]
            str_purity_i = str_purities[i]

            for j in range(nmerge):
                
                if i == j: continue     # skip comparison of region with itself
                if i in del_strs or j in del_strs: continue
                
                str_start_j = str_starts[j]
                str_end_j = str_ends[j]

                if str_start_i >= str_start_j and str_end_i <= str_end_j: # location-i is nested within location-j
                    if str_motifs[i] == str_motifs[j]:
                        del_strs.append(i)
                    elif str_purity_i <= str_purities[j] or \
                       str_motifs[i] == str_motifs[j] or \
                       expected_diff(str_motifs[j], str_purities[j], str_motifs[i], str_lengths[i], str_purities[i]):
                        del_strs.append(i) #lower purity, similar motif size, 
                    else:
                        nested_relations[i] = j
                
                else:
                    if str_start_j <= str_start_i and str_start_i <= str_end_j:     # STR-j upstream of STR-i and overlapping
                        
                        # STR-i and STR-j have the same motif ~ Could happen if they are identified from different motif shifts
                        if str_motifs[i] == str_motifs[j]:
                            merge_cigar, merge_purity = merge_repeats(str_ends[j], str_starts[i], str_cigars[j], str_cigars[i])
                            del_strs.append(j)
                            str_ends[i] = str_ends[j]; str_cigars[i] = merge_cigar; str_purities[i] = merge_purity
                        
                        # STR-i is shorter than STR-j ~ Unique length of STR-i is shorter than STR-i motif length or less than 3 bp
                        elif (str_lengths[i] < str_lengths[j]) and ((str_end_i - str_end_j) < len(str_motif_i) or ((str_end_i - str_end_j) < 3)):
                            del_strs.append(i)      # filter out STR
                        
                        # STR-i length equals STR-j length ~ choose repeat with higher purity
                        elif (str_lengths[i] == str_lengths[j]) and str_purities[i] < str_purities[j]:
                            del_strs.append(i)
                        
                        # retain if both top conditions are false
                        else:
                            overlap_relations[i].add(j)
                            overlap_relations[j].add(i)
                    
                    elif str_start_j <= str_end_i and str_end_i <= str_end_j:       # STR-i upstream of STR-j and overlapping

                        if str_motifs[i] == str_motifs[j]:
                            merge_cigar, merge_purity = merge_repeats(str_ends[i], str_starts[j], str_cigars[i], str_cigars[j])
                            del_strs.append(j)
                            str_ends[i] = str_ends[j]; str_cigars[i] = merge_cigar; str_purities[i] = merge_purity
                        
                        # STR-i is shorter than STR-j ~ Unique length of STR-i is shorter than STR-i motif length or less than 3 bp
                        elif (str_lengths[i] < str_lengths[j]) and (((str_start_j - str_start_i) < len(str_motif_i)) or ((str_start_j - str_start_i) < 3)):
                            del_strs.append(i)
                        
                        # STR-i length equals STR-j length ~ choose repeat with higher purity
                        elif (str_lengths[i] == str_lengths[j]) and str_purities[i] < str_purities[j]:
                            del_strs.append(i)
                        
                        # retain if both top conditions are false
                        else:
                            overlap_relations[i].add(j)
                            overlap_relations[j].add(i)

                if str_end_i < str_starts[j]:   # if STR-j is beyond STR-i
                    break

        if nmerge - len(set(del_strs)) > 1:     # if more than 1 STRs is retained
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

        elif len(str_starts) - len(set(del_strs)) == 1:     #if only 1 STR is retained
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
    args = parse_args()

    bed_file = args.bed
    gff_file = args.gff

    process_bed(bed_file, gff_file)
