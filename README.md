<p align=center>
    <img src="./lib/logo_rounded.png" alt="ribbit-logo" style="width:220px; border-radius: 20%"/>
</p>

<h1 align=left style="font-size: 45px; padding-left: 20px; padding-bottom: 0px">ribbit</h1>

<p style="font-size: 20px">
Ribbit is a tool to identify tandem repeats of variable motif sizes. The algorithm
converts DNA sequences to 2-bit format and uses basic bit operations to identify  tandem repeat sequences. <br>
</p>

## Usage for testing seeds
```
make
./ribbit -i genome.fa > genome.ribbit.seeds
```

## Overview of ribbit algorithm

Ribbit identifies repetitive sequences by identifying sequence regions with high density of periodically matching nucleotides.
Ribbit leverages bit operations to detect periodically matching nucleotides in a sequence. First we convert the sequence
into bit sequences with each nucleotide denoted a combination of 2 bits. The bit sequences are converted to shift XORs 
which sets positions with periodically matching nucleotides as 1s. Ribbit parses these shift XOR sequences of all the 
shifts to identify DNA tandem repeats of specific periodicity.

NOTE: Sequences longer than 10MB are chunked into 10MB sequences with an overlap of 2KB to reduce memory footprint.

### Identification of seeds of perfect repeats

Ribbit first identifies stretches of perfect seeds in all the shift XORs of the desired motif sizes. This step is required
to not miss out on perfect repeat sequences based on the rules of the binomial distribution rules used to identify
imperfect repeat sequences.

### Identification of seeds of repeat sequences with only mismatch errors

In the next step ribbit identifies stretches of potentially repetitive sequences with mismatch errors between the motifs. 
This again is a step to rescue repetitive sequences that might be missed because of the thresholds used in the anchor seed
identification.

### Identification of anchor seed repetitive sequences

The identification of anchor seeds relies on `minimumNumberofSuccesses` function on line 175 of `parse_anchor_shiftxor.cpp`

The function `minimumNumberofSuccesses` is defined in `binomial_thresholds.cpp`.