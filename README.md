<p align=center>
    <img src="./lib/logo_rounded.png" alt="ribbit-logo" style="width:220px; border-radius: 20%"/>
</p>

<h1 align=left style="font-size: 45px; padding-left: 20px; padding-bottom: 0px">ribbit</h1>

<p style="font-size: 20px">
Ribbit is a tool to identify tandem repeats of variable motif sizes. The algorithm
converts DNA sequences to 2-bit format and uses basic bit operations to identify  tandem repeat sequences. <br>
</p>

<h1 align=left style="font-size: 45px; padding-left: 20px; padding-bottom: 0px">Table of Content</h1>
<h2 style="font-size: 35px; padding-left: 20px;">Table of Contents</h2>
<ol style="font-size: 18px; padding-left: 40px;">
    <li><a href="#installation">Installation</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#inputs-and-outputs">Inputs and Outputs</a></li>
    <li><a href="#faq">FAQ</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
</ol>



<h2 style="font-size: 35px; padding-left: 20px;">Installation</h2>
<p style="font-size: 18px; padding-left: 20px;">
    To install Ribbit, clone the repository and install the dependencies using the following commands:
</p>
<pre><code>git clone https://github.com/SowpatiLab/ribbit
cd ribbit
</code></pre>

<h2 style="font-size: 35px; padding-left: 20px;">Usage</h2>
<p style="font-size: 18px; padding-left: 20px;">
    Here’s a basic usage example:
</p>
<pre><code> ./ribbit [options] -i sequence.fasta --output results.bed</code></pre>
</p>
    To view detailed help information
<pre><code> ./ribbit -h </code></pre>
    The output would be given as folllowing.
</p>
<pre style="font-size: 16px; padding-left: 20px; background-color: #f8f8f8; padding: 10px; border-radius: 5px;">
  -h [ --help ]                 Ribbit tool identifies short tandem repeats 
                                with allowed levels of impurity.
  -i [ --input-file ] arg       File path for the input fasta file.
  -o [ --output-file ] arg      File path for the output file.
  -m [ --min-motif-length ] arg The minimum length of the motif of the repeats 
                                to be identified. Default: 2
  -M [ --max-motif-length ] arg The maximum length of the motif of the repeats 
                                to be identified. Default: 100
  -p [ --purity ] arg           Threshold value for the continuous number of 
                                ones found in a seed. Default: 0.85
  -l [ --min-length ] arg       The minimum length of the repeat. Default: 12
  --min-units arg               The minimum number of units of the repeat. Can 
                                be an integer value for cutoff across all motif
                                sizes, or a tab-separated file with two columns: 
                                the first is the motif size and the second is 
                                the unit cutoff. Default: 2
  --perfect-units arg           The minimum number of complete units of the 
                                repeat. Can be an integer value for cutoff 
                                across all motif sizes, or a tab-separated file 
                                with two columns: the first is the motif size and 
                                the second is the unit cutoff. Default: 2
</pre>

<h2 style="font-size: 35px; padding-left: 20px;">Usage</h2>
<p style="font-size: 18px; padding-left: 20px;">
</p>

```-i or --input```
<div style="border: 1px solid #333; padding: 15px; background-color: #f8f8f8; border-radius: 5px;">
    <p><strong>Expects:</strong> <code>STRING</code> (to be used as filename)</p>
    <p>The input file must be a valid FASTA file.</p>
</div>

```-o or --output```
<div style="border: 1px solid #333; padding: 15px; background-color: #f8f8f8; border-radius: 5px;">
    <p><strong>Expects:</strong> <code>STRING</code> (to be used as filename)</p>
    <p>The output for ribbit is <code>.bed</code> file.</p>
    </div>
</p>

### bed file output columns

| S.No | Column           | Description                                                                                  |
|------|------------------|----------------------------------------------------------------------------------------------|
| 1    | Chromosome       | Chromosome or Sequence Name as specified by the first word in the FASTA header               |
| 2    | Repeat Start     | 0-based start position of SSR in the Chromosome                                              |
| 3    | Repeat Stop      | End position of SSR in the Chromosome                                                        |
| 4    | Repeat Class     | Class of repeat as grouped by their cyclical variations                                      |
| 5    | Repeat Length    | Total length of identified repeat in nt                                                      |
| 6    | Motif count      | Number of complete motifs in the STR                                                         |
| 7    | Purity           | Purity of STR region (perfect STR = 1)                                                       |
| 7    | Repeat Strand    | Strand of SSR based on their cyclical variation                                              |
| 8    | CIGAR            | Representing type of imperfections.                                                          | 

</p>

``` -m or --min-motif-length ```
</p>
    The minimum length of the motif of the repeats to be identified.
    
``` -M or --max-motif-length ```
</p>
    The maximum length of the motif of the repeats to be identified.

``` -p or --purity ```
    TEXT
