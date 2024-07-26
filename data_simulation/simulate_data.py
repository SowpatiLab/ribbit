import pandas as pd
import random
from collections import defaultdict, Counter
from tqdm import tqdm
import sys
import argparse

NUCLEOTIDES = {'A', 'C', 'T', 'G'}
MUTATION_TYPES = ['I', 'D', 'S']
PROPORTIONAL_MUTATION_TYPES = ['S']*80 + ['I']*10 + ['D']*10
BUFFER_SEQ = 'GACGTGGTCCCTACTCTCATCTTCAGAGACAAGGTTTACACTGGAAGCCTCTAGGGCAAATGGCTTTTATGATATATAGT\
GAAAAGGGACAGATCACTTAGACTGTCTTCAAAGGAGAACATAATTCTTCTGTTCATATGTCCTCTACTACTTAGGGTCT\
TTAGCAAAATCCTTTATAAGGCAAAAAACGTGCCTGTGTATCCACCTGTAGAATTTAGAGATAGTTTAAATACAGGAAGA\
ATAGCTTCTGCTATAGAGAAAGCCAACACATTTCCTTATAGTTACAAAATGTGTTCGGTAATATCTTCCCATTATATGTG\
TGTTTTATTTCAGCTTGCCTGAATGGAGAGCAAACAGCCTCAGAGGTGTCATAGGTTCTTTTAAGTCCCTTGACCATTTG\
GGGACCAGCTACTCTTTATTGGAAGGAAGATATTTAAGAGAATTCTTTGTTATTCCAAGGAAACTAAATAGTTGTAAAGG\
GACTTTTCTCCTAGGAATTAAATCTTACATAGCAACTGCATACGAATTAAAAGCAGAGTCAAAATTA'


def choose_num_units(motif_size):
    if motif_size == 2:    return random.randint(6, 100)
    elif motif_size == 3:  return random.randint(4, 100)
    elif motif_size <= 50: return random.randint(3, 100)
    else: return random.randint(2, 10)


def mutate_repeat(positions, types):
    
    mutations_info = []; mut_seq = ''; x = 0
    for pos, type in zip(positions, types):
        mut_seq += repeat_seq[x:pos]
        
        if type == 'D':
            mutations_info.append([type, str(pos), repeat_seq[pos]])
            x = pos + 1

        elif type == 'S':
            ori_nuc = repeat_seq[pos]
            sub_nuc = random.choice(list(NUCLEOTIDES - {ori_nuc}))
            mut_seq += sub_nuc
            mutations_info.append([type, str(pos), f'{ori_nuc}/{sub_nuc}'])
            x = pos + 1

        elif type == 'I':
            ins_nuc = random.choice(list(NUCLEOTIDES))
            mut_seq += ins_nuc
            mutations_info.append([type, str(pos), ins_nuc])
            x = pos

    mut_seq += repeat_seq[x:]
    
    return mut_seq, mutations_info


def parse_args():
    parser = argparse.ArgumentParser(description="Tandem Repeat Simulator")

    parser.add_argument('-l', '--num-locations', type=int, default=1000, help='number of repeat locations in the simulated dataset. default: 1000')
    parser.add_argument('-o',  '--out-prefix', type=str, default='', help='prefix to be used for generating output files. [default: id generated]')
    parser.add_argument('--min-purity', type=float, default=0.85, help='minimum total purity (0.5-1) of the repeat generated. default: 0.85')
    parser.add_argument('--max-purity', type=float, default=0.95, help='maximum total purity (0.5-1) of the repeat generated. default: 0.95')
    parser.add_argument('--motif-purity', type=float, default=0.75, help='minimum purity of each motif with the consensus motif. default: 0.85')

    parser.add_argument('-m', '--min-motif-size', type=int, default=2,   help='minimum motif length of repeats. default: 2')
    parser.add_argument('-M', '--max-motif-size', type=int, default=100, help='maximum motif length of repeats. default: 100')    

    parser.add_argument('-C',  '--max-chrom-size', type=int, default=250*(10**6), help='maximum chromosome size in the fasta. default: 2500000000')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if args.out_prefix == '':
        args.out_prefix = "%06x" % random.randint(0, 0xFFFFFFFFFF)
        args.out_prefix = args.out_prefix.upper()
    
    print(f"File prefix: {args.out_prefix}")

    proportions_file = './proportions.tsv' 
    proportions_df = pd.read_csv(proportions_file, sep='\t')
    
    motif_sizes = []
    for index, row in proportions_df.iterrows():
        count = int(row['%_proportion'] * args.num_locations / 100)  
        motif_sizes.extend([int(row['Motif_size'])] * count)

    rmotifs = defaultdict(list)
    with open('./HG38_2-100_motifs_d2d.tsv') as fh:
        for line in fh:
            motif, kmer = line.strip().split('\t')
            kmer = int(kmer)
            rmotifs[kmer].append(motif)

    # Open output files
    bed = open(f'sim_{args.out_prefix}.bed', 'w')
    faout = open(f'sim_{args.out_prefix}.fa', 'w')
    print(f'>{args.out_prefix}_1', file=faout)

    fasta_seq = ''
    position = 0
    chrom_num = 1
    seq_length = 0

    idx_digits = len(str(args.num_locations))
    ridx = 0; repeat_id = ''

    min_impurity = int(100*(1-args.max_purity))
    max_impurity = int(100*(1-args.min_purity))

    # Iterate over each locus
    for _ in tqdm(range(args.num_locations), bar_format='{l_bar}{bar:100}{r_bar}{bar:-10b}'):
        repeat_id = 'R' + format(ridx, f'0{idx_digits}d')

        # Add buffer sequence
        buffer_size = random.randint(500, 3000)
        buffer = (BUFFER_SEQ* ((buffer_size//len(BUFFER_SEQ)) + 1))[:buffer_size]
        fasta_seq += buffer
        position += len(buffer)
        seq_length += len(buffer)

        motif_size = random.choice(motif_sizes)
        runits = choose_num_units(motif_size)
        suffix_len = int((random.randint(0, 9) / 10) * motif_size)
        rlength = (motif_size * runits) + suffix_len
        if suffix_len > 0.75*motif_size: runits += 1
        rmotif = random.choice(rmotifs[motif_size])
        repeat_seq = (rmotif * (runits+1))[:rlength]
        
        # Calculate number of mutations based on percentage
        impurity = random.randint(min_impurity, max_impurity)
        num_mutations = int((impurity/100) * rlength)  # Ensure at least one mutation        

        max_motif_mutations = max(1, int((1-args.motif_purity)) * motif_size)
        max_mutations = min(num_mutations, max_motif_mutations*runits)
        motif_mutation_counter = Counter() # to keep the track of mutations introduced in the motif
        mpositions = []; mtypes = []

        # Generate mutations ensuring the distribution
        
        while len(mpositions) < max_mutations:
            mposition = random.randint(1, rlength - 1)
            if mposition in mpositions: continue
            motif_index = mposition // motif_size
            if motif_mutation_counter[motif_index] < max_motif_mutations:
                mpositions.append(mposition)
                mtypes.append(random.choice(PROPORTIONAL_MUTATION_TYPES))
                motif_mutation_counter[motif_index] += 1
        sorted_indices = sorted(list(range(len(mpositions))), key=lambda x: mpositions[x])
        mpositions_sorted = []
        mtypes_sorted = []
        for i in sorted_indices:
            mpositions_sorted.append(mpositions[i])
            mtypes_sorted.append(mtypes[i])
        # if len(mpositions) < num_mutations:
        #     print(f"Warning: Only {len(mpositions)} out of {num_mutations} mutations could be placed at R{repeat_id}", file=sys.stderr)
    
        mut_seq, mutations_info = mutate_repeat(mpositions_sorted, mtypes_sorted)
        fasta_seq  += mut_seq
        seq_length += len(mut_seq)
        
        mutations_str = ';'.join(['|'.join(a) for a in mutations_info])
        
        print(f'{args.out_prefix}_{chrom_num}', position, position + len(mut_seq), repeat_id, len(mut_seq), len(rmotif), rmotif, mutations_str, sep='\t', file=bed)
        ridx += 1
        
        while len(fasta_seq) >= 80:
            print(fasta_seq[:80], file=faout)
            fasta_seq = fasta_seq[80:]
        position += len(mut_seq)

        if seq_length >= args.max_chrom_size:
            buffer_size = random.randint(500, 3000)
            buffer = (BUFFER_SEQ* ((buffer_size//len(BUFFER_SEQ)) + 1))[:buffer_size]
            fasta_seq += buffer
            while len(fasta_seq) >= 80:
                print(fasta_seq[:80], file=faout)
                fasta_seq = fasta_seq[80:]
            if len(fasta_seq) > 0: print(fasta_seq, file=faout)
            fasta_seq = ''
            seq_length = 0
            position = 0
            chrom_num += 1
            print(f'>{args.out_prefix}_{chrom_num}', file=faout)

    buffer_size = random.randint(500, 3000)
    buffer = (BUFFER_SEQ* ((buffer_size//len(BUFFER_SEQ)) + 1))[:buffer_size]
    fasta_seq += buffer
    while len(fasta_seq) > 80:
        print(fasta_seq[:80], file=faout)
        fasta_seq = fasta_seq[80:]
    print(fasta_seq, file=faout)

    bed.close()
    faout.close()

    print(f'Simulation completed for {args.out_prefix} with {args.num_locations} loci.')
