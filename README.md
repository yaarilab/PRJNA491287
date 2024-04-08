A pipeline for BCR repertoire libraries from  - Biomed IGHV


Library preperation and sequencing method:

The sequences were amplified specific Biomed IGHV primers in the framework 1 (FR1) region and isotype-specific primers for IgG, IgM and IgA.
The generated libraries were then sequenced with Illumina MiSeq.


Input files:

* Pair-end reads Sample_R1.fastq and Sample_R2.fastq 
* primers sequences - constant region specific , v specific primers.

Output file:

1. Sample.fasta - Two processed fasta files, one for Heavy chain sequence, and one for Light
2. log tab file for each steps
3. report for some of the steps


Pipeline container:

* Docker: immcantation/suite:4.3.0


Sequence processing steps:

converting the fastq files into fasta files for alignmetn
