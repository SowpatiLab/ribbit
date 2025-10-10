<p align=center>
    <img src="./lib/logo_rounded.png" alt="ribbit-logo" style="width:180px; border-radius: 20%"/>
</p>

<h1 align=left style="font-size: 38px; padding-left: 10px; padding-bottom: 0px">ribbit</h1>

<!-- <p style="font-size: 16px"> -->
Ribbit identifies tandem repeat (TR) regions of variable motif sizes in DNA sequences. Tandem repeats are DNA segments consisting of two or more nearly identical copies of a motif occurring contiguously. Ribbit is designed to improve motif decomposition and resolve complex repeat structures, accurately detecting TRs with motif sizes up to 100 bp.

The algorithm converts DNA sequences to 2-bit format and uses basic bit operations to deliniate potential repetitive stretches of a certain periodicity. The DNA sequence of a potential repetitive stretch is decomposed to identify a representative motif of identified peridicity. There are two different approaches to identifying the representative motif depending the on the peridicity. The sequence in the potential stretch is aligned with the perfect repeat of the identified representative motif to calculate the purity of the repeat. The potential sequence is either trimmed or dropped based on the user desired purity and minimum repeat length thresholds. A given sequence is idependently processed for potential repeats of all user-desired motif sizes. The results from each idependent search are merged and compared to resolve for nested/overlapping interpretations of a sequence as repetitive sequence of different periodicities. Overlapping tandem repeat interpretations of different motif sizes are retained/dropped based on preference for more pure stretches. 

The conversion of DNA to 2-bit stretches results in fast identification of potential repetitive stretches and allows the time for careful motif decomposition of repeat sequence. Ribbit provides a comprehensive bed file as an output and takes about 5-7 secs to resolve an MB of DNA sequence. The program can also be run on multi-threaded mode making it ideal for processing large genomes.

<h2 style="font-size: 30px; padding-left: 20px;">Table of Contents</h2>
<ol style="font-size: 16px; padding-left: 40px;">
    <li><a href="#compiling">Compiling</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#options-description">Options description</a></li>
    <li><a href="#Output">Output</a></li>
    <li><a href="#citation">Citation</a></li>
    <li><a href="#authors">Authors</a></li>
    <li><a href="#contact">Contact</a></li>
</ol>



## Compiling
<p style="font-size: 18px">
To install Ribbit, clone the repository and install the dependencies using the following commands:
</p>

### Installing dependencies

#### 1. Install boost library
```bash
sudo apt-get install boost
```

### Instruction for compiling

```bash
git clone https://github.com/SowpatiLab/ribbit.git
cd ribbit
make
```

## Usage
<p style="font-size: 18px">
    Here’s a basic usage example:
</p>

```bash
$ ./ribbit [options] -i sequence.fasta -o results.bed
```

<p style="font-size: 18px">
    To view detailed help information
</p>

```bash
./ribbit -h

Ribbit: accurate identification of tandem repeats and annotation of complex tandem repeat sequences in genomes
Version 1.0.0

Options for running the tool:
  -h [ --help ]                 Ribbit detects tandem repeat regions in DNA, accurately resolving 
                                complex repeat structures and motif sizes up to 100 bp.
  -i [ --input-file ] arg       File path for the input fasta file.
  -o [ --output-file ] arg      File path for output file. Default: {input-file}.ribbit
  -m [ --min-motif-length ] arg The minimum length of the motif of the TR loci. Default: 2
  -M [ --max-motif-length ] arg The maximum length of the motif of the TR loci. Default: 100
  -p [ --min-purity ] arg       The minimum allowed purity of repeat sequence. Purity is calculated
                                as the (matches/(matches+mismatches+indels)) in the alignment of 
                                region sequence to perfect repeat of consensus motif. Default: 0.8
  -q [ --min-motif-purity ] arg Minimum purity of each motif with consensus motif. Calculated as 
                                the average of (matches/(matches+mismatches+indels)) for each motif
                                length in the alignment of region sequence to perfect repeat of 
                                consensus motif. Default: 0.8
  -l [ --min-length ] arg       The minimum length of the repeat. Input can be an integer or a 
                                tab-separated file with two columns of motif length and the length 
                                cutoff. Default: 12 for STRs (motif length <= 6), 2*(motif length) 
                                for others.
  --min-units arg               The minimum number of units of the repeat. Input can be a integer 
                                or a tab-separated file with two columns, first is the motif size 
                                and second unit cutoff. Default: 2 for all motif sizes.
  --perfect-units arg           The minimum number of complete units with 100% match with the 
                                consensus motif in the repeat. Input can be an integer or a 
                                tab-separated file with two columns of the motif length and the 
                                unit cutoff. Default: 2
  --cigar                       Include cigar string of the alignment of the sequence with the 
                                perfect repeat of the consensus motif in the output. Default: 
                                false.
  -t [ --threads ] arg          Number of threads to be used for running. Default: 1
```

## Options description


`-p or --min-purity`

 

`-o or --output` file path of the output file. default: standard output


## Output
### Output columns

| S.No | Column           | Description                                                   |
|------|------------------|---------------------------------------------------------------|
| 1    | Chromosome       | chromosome or sequence name as specified in fasta header      |
| 2    | Repeat Start     | 0-based start position of tandem repeat in the chromosome     |
| 3    | Repeat Stop      | end position of tandem repeat in the chromosome               |
| 4    | Motif            | sequence of unit motif of tandem repeat                       |
| 5    | Purity           | purity of the repeat                                          |
| 6    | Motif length     | length of the unit motif in bases                             |
| 7    | Repeat Length    | length of the repeat based on the coordinates                 |
| 8    | Repeat Units     | number of unit motifs of the TR                               |
| 9    | Info             | additional information of the repeat region including the CIGAR string of alignment and information of the nested TRs in the region. | 


### Output file example

| chromosome | start | end  | motif | purity    | motif_length | repeat_length   | repeat_units | info |
|------------|-------|------|-------|-----------|--------------|-----------------|--------------|------|
| test       | 292   | 305  | AC    |    0.93   | 2            | 13              | 6            | I    |
| test       | 481   | 496  | AGC   |    0.82   | 3            | 15              | 5            | I    |
| test       | 827   | 843  | GCCCAGGT | 1.00   | 8            | 16              | 2            | I    |
| test       | 1017  | 1032 | TGCGGAG  | 0.93   | 7            | 15              | 2            | I    |
| test       | 1508  | 1523 | AGGC     | 0.81   | 4            | 15              | 3            | I    |
| test       | 1682  | 1833 | TGGAGGGTGGGGCCAAATGGAAGTGGGCGGGGCTGTGG | 0.90 | 38 | 151 | 3      | I    |
| test       | 1863  | 1890 | AGGGC    | 0.83   | 5            | 27              | 5            | M:1869-1890-6-0.81:AGGGGC |
| test       | 2182  | 2197 | CCGGT    | 0.93   | 5            | 15              | 3            | I    |
| test       | 2277  | 2296 | TGGCCTCC | 1.00   | 8            | 19              | 2            | I    |
</p>

## Citation
<p style="font-size: 16px">
If you found ribbit useful, we would appreciate it if you could cite our manuscript: <a href="https://doi.org/10.1101/2025.02.06.636828">Ribbit: Accurate identification and annotation of complex tandem repeat sequences in genomes</a>
</p>

## Authors
Anukrati Sharma <br>
Akshay Kumar Avvaru

## Contact
For queries or suggestions, please contact:
<br>Akshay Kumar Avvaru - avvaruakshay@gmail.com
<br>Divya Tej Sowpati - tej@ccmb.res.in
