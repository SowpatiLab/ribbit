import pandas as pd
from random import choice, randint
from collections import defaultdict
from tqdm import tqdm
import sys

nucs = {'A', 'C', 'T', 'G'}
mut_types = ['I', 'D', 'S']

non_repeat_buffer = 'CTATCTGTTAGCCTGGGAGAAAACTATAATTTGTTTTATAAGGAATTCAACCAGATAGATAGGAATTAGTCAGATTCTAACATACAAGATAGATTTTGAACAAAACAATAGGCTTGAGAGGCAGCCAGGCATTGTTGAGAAAGTTTGTTCCAAATAAATTATTTGCTCCTTTGGTATATGAACAACTTAGCTTGAAGTGCTCGTTCTCTGAGAACACCCCTCTAGCATAGAGAAGATGACTCTTCTGCTACAATCCCCAGGCGAGCTCCATGGCCGCAGAGGGTAGTAAAAGAAAAGCCCGCCCCTCTTTTGGCTCTAGGGATTGTGAGATAAAGCGGGTGGGAGAGTAAATCAGTGAAGCAATCGAACATCTCTCTCCCCATTTACTATAAGGTAAATAGTAAGCTCCAAGAAGCCGGTTAGGACTTGGGAGGTGGCCAGCTTACAGGCTGCCCTAAACTAACCTGGGTCACAGTCAAATCAACCAGAGAATGTAGTGAGCTCAGGCAGAAGAAACTCTGAGGTCATCCTTTCCTGAGGCATGCTTGGTTCCCTTGTGCACAGACAAATTCCAGAAAGAACTGGAAATCTGGTTAGGAAACTTGGCAAAGTGGGTGTCACTGTATTGGGGAGGGTGAAACAATTCTGTAGCAGGAGTTGCAAGGCAGAACCCTCGAACTAGGGCCTAATGATAGTCACAGCCAAACCTCATAAGGGTCTCTAACACCAATAGGAGCTGAGGTGAGGCAGAGCCAGCCAAGGGTCACTGCTGAGCCCTGC'

proportions_file = '/data/ccmb/anukrati/hg38/proportions.tsv' 
proportions_df = pd.read_csv(proportions_file, sep='\t')

print(proportions_df.columns)

if len(sys.argv) == 1:
    num_loci = 10 
    num_mutations = 5  
else:#
    num_loci = int(sys.argv[1])
    num_mutations = int(sys.argv[2])

id = "%06x" % randint(0, 0xFFFFFFFFFF)
id = id.upper()

motif_sizes = []
for index, row in proportions_df.iterrows():
    count = int(row['%_proportion'] * num_loci / 100)  
    motif_sizes.extend([int(row['Motif_size'])] * count)

rmotifs = defaultdict(list)
with open('/data/ccmb/anukrati/hg38/HG38_2-100_motifs_d2d.tsv') as fh:
    for line in fh:
        motif, kmer = line.strip().split('\t')
        kmer = int(kmer)
        rmotifs[kmer].append(motif)

# Open output files
bed = open(f'sim_{id}.bed', 'w')
faout = open(f'sim_{id}.fa', 'w') 
print(f'>{id}', file=faout)

fasta_seq = ''
position = 0
chr_num = 1
seq_length = 0
max_seq_length = 200 * 10**6  # 200 Mbp

# Iterate over each locus
for _ in tqdm(range(num_loci)):
    # Add buffer sequence
    buffer_size = randint(3, 12)
    buffer = non_repeat_buffer*buffer_size
    fasta_seq += buffer
    position += len(buffer)
    seq_length += len(buffer)

    motif_size = choice(motif_sizes)
    runits = 100
    if motif_size == 2:
        runits = randint(6, 100)
    elif motif_size == 3:
        runits = randint(4, 100)
    elif motif_size <= 50:
        runits = randint(3, 100)
    else:
        runits = randint(2, 10)
    
    suffix_len = int((randint(0, 9) / 10) * motif_size)
    rlength = motif_size * runits + suffix_len
    
    rmotif = choice(rmotifs[motif_size])
    
    repeat_seq = (rmotif * ((rlength // len(rmotif)) + 1))[:rlength]
    
    mtypes = [mut_types[randint(0, 2)] for _ in range(num_mutations)]
    mut_positions = [randint(1, rlength-1) for _ in range(num_mutations)]
    mutations_info = []
    mut_seq = ''

    x = 0
    for mut_position, mtype in zip(mut_positions, mtypes):
        mut_seq += repeat_seq[x:mut_position]
        
        if mtype == 'D':
            mutations_info.append([mtype, str(mut_position), repeat_seq[mut_position]])
            x = mut_position + 1
        elif mtype == 'S':
            ori_nuc = repeat_seq[mut_position]
            sub_nuc = choice(list(nucs - {ori_nuc}))
            mut_seq += sub_nuc
            mutations_info.append([mtype, str(mut_position), f'{ori_nuc}/{sub_nuc}'])
            x = mut_position + 1
        elif mtype == 'I':
            ins_nuc = choice(list(nucs))
            mut_seq += ins_nuc
            mutations_info.append([mtype, str(mut_position), ins_nuc])
            x = mut_position
    
    mut_seq += repeat_seq[x:]
    fasta_seq += mut_seq
    seq_length += len(mut_seq)
    
    mutations_str = ';'.join(['|'.join(a) for a in mutations_info])
    
    print(id, position, position + len(mut_seq), rmotif, len(mut_seq), mut_seq, mutations_str, sep='\t', file=bed)
    
    while len(fasta_seq) >= 80:
        print(fasta_seq[:80], file=faout)
        fasta_seq = fasta_seq[80:]
    
    position += len(mut_seq)

    if seq_length >= max_seq_length:
        print(fasta_seq, file=faout)
        fasta_seq = ''
        seq_length = 0
        chr_num += 1
        id = f'{id}_chr{chr_num}'
        print(f'>{id}', file=faout)

buffer_size = randint(100,200)
fasta_seq += ''.join(choice(non_repeat_buffer) for _ in range(buffer_size))
while len(fasta_seq) > 80:
    print(fasta_seq[:80], file=faout)
    fasta_seq = fasta_seq[80:]
print(fasta_seq, file=faout)

bed.close()
faout.close()

print(f'Simulation completed for {id} with {num_loci} loci.')
