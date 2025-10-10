<p align=center>
    <img src="./lib/logo_rounded.png" alt="ribbit-logo" style="width:180px; border-radius: 20%"/>
</p>

<h1 align=left style="font-size: 38px; padding-left: 10px; padding-bottom: 0px">ribbit</h1>

<!-- <p style="font-size: 16px"> -->
Ribbit is a software tool to identify tandem repeat (TR) regions of variable motif sizes in DNA sequences. Tandem repeats in DNA are sequences with a two or more near exact copies of a motif occuring contigously. Ribbit is designed with a focus on improving motif decomposition, and resolving complex structures in TR sequences. Ribbit can resolve a DNA sequence for TRs of motif size upto 100 bp.

The algorithm converts DNA sequences to 2-bit format and uses basic bit operations to deliniate potential repetitive stretches of a certain periodicity. The DNA sequence of a potential repetitive stretch is decomposed to identify a representative motif of identified peridicity. There are two different approaches to identifying the representative motif depending the on the peridicity. The sequence in the potential stretch is aligned with the perfect repeat of the identified representative motif to calculate the purity of the repeat. The potential sequence is either trimmed or dropped based on the user desired purity and minimum repeat length thresholds. A given sequence is idependently processed for potential repeats of all user-desired motif sizes. The results from each idependent search are merged and compared to resolve for nested/overlapping interpretations of a sequence as repetitive sequence of different periodicities. Overlapping tandem repeat interpretations of different motif sizes are retained/dropped based on preference for more pure stretches. 

The conversion of DNA to 2-bit stretches results in fast identification of potential repetitive stretches and allows the time for careful motif decomposition of repeat sequence. Ribbit provides a comprehensive bed file as an output and takes about 5-7 secs to resolve an MB of DNA sequence. The program can also be run on multi-threaded mode making it ideal for processing large genomes.

<!-- </p> -->

<h2 style="font-size: 35px; padding-left: 20px;">Table of Contents</h2>
<ol style="font-size: 18px; padding-left: 40px;">
    <li><a href="#installation">Installation</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#inputs-and-outputs">Inputs and Outputs</a></li>
    <li><a href="#citation">Citation</a></li>
    <li><a href="#contact">Contact</a></li>
</ol>



## Installation
<p style="font-size: 18px">
To install Ribbit, clone the repository and install the dependencies using the following commands:
</p>

### Installing dependencies

#### 1. Install boost library
```
sudo apt-get install boost
```

<br>

### Compiling ribbit

```
git clone https://github.com/SowpatiLab/ribbit.git
cd ribbit
make
```
<br>

## Usage
<p style="font-size: 18px">
    Here’s a basic usage example:
</p>

```
./ribbit [options] -i sequence.fasta --output results.bed
```

<p style="font-size: 18px">
    To view detailed help information
</p>

```
./ribbit -h

Ribbit: A fast and accurate tandem repeat finder.
Version 1.0.0
Below are the running options for the tool.:
  -h [ --help ]                 Ribbit is designed to identify tandem repeats 
                                in DNA sequences with specific focus           
                                                            on annotating 
                                complex TR loci.
  -i [ --input-file ]           File path for the input fasta file.
  
  -o [ --output-file ]          File path for the input fasta file. Default: sys.stdout

  -m [ --min-motif-length ]     Minimum length of the motif of identified TR loci. Default: 2
  
  -M [ --max-motif-length ]     Maximum length of the motif of identified  TR loci. Default: 100
  
  -p [ --min-purity ]       arg The minimum allowed purity of complete repeat.  Default: 0.8
  
  -q [ --min-motif-purity ] arg Minimum match of each motif with consensus 
                                motif. Default: 0.8
  
  -l [ --min-length ]       arg The minimum length of the repeat. Default: 12
  
  --min-units               arg The minimum number of units of the repeat. Can 
                                be a integer value, for cutoff across all motif
                                sizes.                                         
                                       Tab separated file with two columns, 
                                first is the motif size and second unit cutoff.
                                Default: 2
  
  --perfect-units arg           The minimum number of complete units of the 
                                repeat. Can be a integer value, for cutoff 
                                across all motif sizes.                        
                                                        Tab separated file with
                                two columns, first is the motif size and second
                                unit cutoff. Default: 2
  
  --cigar                       Include cigar string in the output. Default is 
                                off.
  
  -t [ --threads ] arg          Number of threads to be used for running. 
                                default: 1
```

## Inputs and Outputs
<p style="font-size: 18px; padding-left: 20px;">

```-i or --input```
<div style="border: 1px solid #333; padding: 15px; border-radius: 5px; margin-bottom: 10px">
    <p><strong>Expects:</strong> <code>STRING</code> (to be used as filename)</p>
    <p>The input file must be a valid FASTA file.</p>
</div>


```-o or --output```
<div style="border: 1px solid #333; padding: 15px; border-radius: 5px;">
    <p><strong>Expects:</strong> <code>STRING</code> (to be used as filename)</p>
    <p>The output for ribbit is <code>.bed</code> file.</p>
    </div>
</p>

#### bed file output columns

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

| Chromosome   | Start   | End     | Motif                         | Purity    | Strand | CIGAR                                                       | Motif Size | Repeat Length | Motif Units   |
|--------------|---------|---------|-------------------------------|-----------|--------|-------------------------------------------------------------|------------|---------------| --------------|
| Test_Seq     | 90196   | 90393   | AC                            | 0.9494    | +      | 3=1X3=1X5=1D82=1X17=1X19=1X31=1I2=1X3=1X21=1I2=             | 2          | 197           | 98            |
| Test_Seq     | 137451  | 137470  | CCCGCT                        | 1         | +      | 19=                                                         | 6          | 19            | 3             |
| Test_Seq     | 136254  | 136401  | GT                            | 0.9127    | +      | 6=1X9=1D20=1D15=1X12=1X5=1X25=1X9=1X7=1X5=1X9=1X10=1X2=1X2= | 2          | 147           | 73            |
| Test_Seq     | 139286  | 139306  | AGTTGCTT                      | 0.95      | +      | 8=1X11=                                                     | 8          | 20            | 2             |
| Test_Seq     | 3538110 | 3538168 | AATAGCAAGAGCCAGAGCTAGAGCAAAG  | 0.8813    | +      | 4=1X1=2I30=1X9=1X5=1X1=1D2=                                 | 8          | 58            | 7             |
| Test_Seq     | 4197438 | 4197487 | CACAGCCAGCT                   | 0.9591    | +      | 26=1X12=1X9=                                                | 11         | 49            | 4             |
| Test_Seq     | 4858037 | 4858050 | CTCTTT                        | 0.9230    | +      | 6=1I6=                                                      | 6          | 13            | 2             |
| Test_Seq     | 5000704 | 5000745 | TATTCGTATGCGTATTC             | 0.9024    | +      | 4=1I22=1X4=2X7=                                             | 17         | 41            | 2             |

</p>

## Citation
<p style="font-size: 16px">
If you found ribbit useful, we would appreciate it if you could cite our manuscript: <a href="https://doi.org/10.1101/2025.02.06.636828">Ribbit: Accurate identification and annotation of complex tandem repeat sequences in genomes</a>
</p>

## Contact
<p style="font-size: 16px">
    For queries or suggestions, please contact:
    <br>Akshay Kumar Avvaru - avvaruakshay@gmail.com
    <br>Divya Tej Sowpati - tej@ccmb.res.in
</p>
