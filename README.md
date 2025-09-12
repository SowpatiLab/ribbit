<p align=center>
    <img src="./lib/logo_rounded.png" alt="ribbit-logo" style="width:220px; border-radius: 20%"/>
</p>

<h1 align=left style="font-size: 45px; padding-left: 20px; padding-bottom: 0px">ribbit</h1>

<p style="font-size: 20px">
This branch is currently to test the seeds using anchor method <br>
NOTE: Sequences longer than 10MB are chunked into 10MB sequences with an overlap of 2KB to reduce memory footprint.
</p>

## Usage for testing seeds
```
make
./ribbit -i genome.fa > genome.ribbit.seeds
```

### Output format
```
#chrom	start	stop	motif_length	seed_length	anchored_bitcount	motif_bitcount	perfect_bitcount
chr1	10000	10155	30	155	152	122	89
chr1	10000	10161	24	161	158	134	109
chr1	10019	10191	25	172	158	83	27
chr1	10001	10448	23	447	425	209	71
chr1	10451	10461	23	10	9	8	5
chr1	10069	10472	32	403	399	144	53
chr1	10001	10473	31	472	452	260	119
chr1	10000	10478	36	478	471	286	161
chr1	10059	10480	38	421	412	153	59
```