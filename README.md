<p align=center>
    <img src="./lib/logo_rounded.png" alt="ribbit-logo" style="width:220px; border-radius: 20%"/>
</p>

<h1 align=left style="font-size: 45px; padding-left: 20px; padding-bottom: 0px">ribbit</h1>

<p style="font-size: 20px">
Ribbit is a tool to identify tandem repeats of variable motif sizes. The algorithm
converts DNA sequences to 2-bit format and uses basic bit operations to identify  tandem repeat sequences. <br>
</p>


<h2 style="font-size: 35px; padding-left: 20px;">Table of Contents</h2>
<ol style="font-size: 18px; padding-left: 40px;">
    <li><a href="#installation">Installation</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#inputs-and-outputs">Inputs and Outputs</a></li>
    <li><a href="#license">Citation</a></li>
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

<h2 style="font-size: 35px; padding-left: 20px;">Inputs and Outputs</h2>
<p style="font-size: 18px; padding-left: 20px;">

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
</p>
    TEXT

## Bed file output example

| Chromosome   | Start   | End     | Motif                           | Motif Size | Location Size | Purity     | Strand | CIGAR  |
|--------------|---------|---------|---------------------------------|------------|---------------|-----------|--------|-----------------------------------------------------|
| Test_Seq | 90196   | 90393   | AC                              | 2       | 197           | 0.949495  | +       | 3=1X3=1X5=1D82=1X17=1X19=1X31=1I2=1X3=1X21=1I2=   |
| Test_Seq | 137451  | 137470  | CCCGCT                          | 6       | 19            | 1         | +      | 19=                                        |
| Test_Seq | 136254  | 136401  | GT                              | 2       | 147           | 0.912752  | +  |6=1X9=1D20=1D15=1X12=1X5=1X25=1X9=1X7=1X5=1X9=1X10=1X2=1X2= |
| Test_Seq | 139286  | 139306  | AGTTGCTT                        | 8       | 20            | 0.95      | +      |   8=1X11=                                   |
| Test_Seq | 3538110 | 3538168 | AATAGCAAGAGCCAGAGCTAGAGCAAAG    | 8     | 58            | 0.881356  | +      |    4=1X1=2I30=1X9=1X5=1X1=1D2=               |
| Test_Seq | 4197438 | 4197487 | CACAGCCAGCT                     | 11     | 49            | 0.959184  | +      |      26=1X12=1X9=                               |
| Test_Seq | 4858037 | 4858050 | CTCTTT                          | 6       | 13            | 0.923077  | +      |    6=1I6=                                    
| Test_Seq | 5000704 | 5000745 | TATTCGTATGCGTATTC               | 17     | 41            | 0.902439  | +      |   4=1I22=1X4=2X7=                            |

</p>

<h2 style="font-size: 35px; padding-left: 20px;">Citation</h2>
<p style="font-size: 18px; padding-left: 20px;">
    Please cite as follows :
</p>
    Ribbit: Accurate identification and annotation of imperfect tandem repeat sequences in genomes
</p>
    Akshay Kumar Avvaru, Anukrati Sharma, Divya Tej Sowpati
    Journal:
    doi:

</p>
<h2 style="font-size: 35px; padding-left: 20px;">Contact</h2>
<p style="font-size: 18px; padding-left: 20px;">
    For queries or suggestions, please contact:
    Akshay Kumar Avvaru - avvaru@ccmb.res.in
    Divya Tej Sowpati - tej@ccmb.res.in

