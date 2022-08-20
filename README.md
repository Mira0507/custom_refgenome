# Building custom reference genome for single cell RNA-seq analysis using `cellranger`

2022-08-19     
Mira Sohn

This workflow demonstrates how to build custom reference genome for analyzing expression of the human transgene p53-R280K in transgenic mice. p53 (TP53) is a tumor suppressor gene having 80% identical cDNA sequences between human and mouse. This approach started from a question whether it's possible to distinguish human transgene p53 and mouse endogenous p53 in a scRNA-seq analysis.

This workflow will be a slightly modified version of the guide [Creating a Reference Package with cellranger mkref](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references) provided by 10X genomics. Here, description about about how to install `cellranger` will be skipped. Visit [Cell Ranger Installation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation) if you need to install the software.


Due to the fact that `cellranger` utilizes [`STAR`](https://github.com/alexdobin/STAR) to build reference index files and align reads to reference genome, genome fasta and annotation GTF files are needed as input, which is identical to what `STAR --runMode genomeGenerate` requires.


## 1. Downloading mouse reference genome files

There are two ways of starting: 

- downloading from original reference files (e.g. [GENCODE release M10](https://www.gencodegenes.org/mouse/release_M10.html))
- downloading from [pre-built reference packages provided by cellrange](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)

This workflow will go with the second option which is a shortcut. Note that it's required to provide your name, email address, and zip code for access. 


```bash

# Download in the directory named reference
wget link/to/refdata-gex-mm10-2020-A.tar.gz -P reference


# workingdir/reference$ ls
# refdata-gex-mm10-2020-A.tar.gz


```

Once it's downloaded, decompress it.

```bash

# Decompress
workingdir/reference$ tar -xvf refdata-gex-mm10-2020-A.tar.gz
refdata-gex-mm10-2020-A/
refdata-gex-mm10-2020-A/fasta/
refdata-gex-mm10-2020-A/fasta/genome.fa
refdata-gex-mm10-2020-A/fasta/genome.fa.fai
refdata-gex-mm10-2020-A/star/
refdata-gex-mm10-2020-A/star/chrNameLength.txt
refdata-gex-mm10-2020-A/star/SA
refdata-gex-mm10-2020-A/star/sjdbInfo.txt
refdata-gex-mm10-2020-A/star/transcriptInfo.tab
refdata-gex-mm10-2020-A/star/sjdbList.fromGTF.out.tab
refdata-gex-mm10-2020-A/star/chrLength.txt
refdata-gex-mm10-2020-A/star/SAindex
refdata-gex-mm10-2020-A/star/genomeParameters.txt
refdata-gex-mm10-2020-A/star/exonInfo.tab
refdata-gex-mm10-2020-A/star/chrName.txt
refdata-gex-mm10-2020-A/star/chrStart.txt
refdata-gex-mm10-2020-A/star/sjdbList.out.tab
refdata-gex-mm10-2020-A/star/exonGeTrInfo.tab
refdata-gex-mm10-2020-A/star/Genome
refdata-gex-mm10-2020-A/star/geneInfo.tab
refdata-gex-mm10-2020-A/genes/
refdata-gex-mm10-2020-A/genes/genes.gtf
refdata-gex-mm10-2020-A/pickle/
refdata-gex-mm10-2020-A/pickle/genes.pickle
refdata-gex-mm10-2020-A/reference.json


# Explore directory structure
workingdir/reference$ tree
.
├── refdata-gex-mm10-2020-A
│   ├── fasta
│   │   ├── genome.fa
│   │   └── genome.fa.fai
│   ├── genes
│   │   └── genes.gtf
│   ├── pickle
│   │   └── genes.pickle
│   ├── reference.json
│   └── star
│       ├── chrLength.txt
│       ├── chrNameLength.txt
│       ├── chrName.txt
│       ├── chrStart.txt
│       ├── exonGeTrInfo.tab
│       ├── exonInfo.tab
│       ├── geneInfo.tab
│       ├── Genome
│       ├── genomeParameters.txt
│       ├── SA
│       ├── SAindex
│       ├── sjdbInfo.txt
│       ├── sjdbList.fromGTF.out.tab
│       ├── sjdbList.out.tab
│       └── transcriptInfo.tab
└── refdata-gex-mm10-2020-A.tar.gz

5 directories, 21 files
```

Here, only the `refdata-gex-mm10-2020-A/fasta/genome.fa` and `refdata-gex-mm10-2020-A/genes/genes.gtf` will be used.

```bash

# Move genome.fa and genes.gtf
workingdir/reference$ cp refdata-gex-mm10-2020-A/genes/genes.gtf .
workingdir/reference$ cp refdata-gex-mm10-2020-A/fasta/genome.fa .


# Ensure the copied files are found in the corrent location
workingdir/reference$ ls
genes.gtf  genome.fa  refdata-gex-mm10-2020-A  refdata-gex-mm10-2020-A.tar.gz


```


Explore what the GTF file looks like. Details about GTF file format are introduced by the page [GFF/GTF File Format - Definition and supported options](https://useast.ensembl.org/info/website/upload/gff.html).

```bash

workingdir/reference$ head genes.gtf   # first few rows
##description: evidence-based annotation of the mouse genome (GRCm38), version M23 (Ensembl 98)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2019-09-06
chr1    HAVANA  gene    3205901 3671498 .       -       .       gene_id "ENSMUSG00000051951"; gene_version "5"; gene_type "protein_coding"; gene_name "Xkr4"; level 2; mgi_id "MGI:3528744"; havana_gene "OTTMUSG00000026353.2";
chr1    HAVANA  transcript      3205901 3216344 .       -       .       gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000162897"; transcript_version "1"; gene_type "protein_coding"; gene_name "Xkr4"; transcript_type "lncRNA"; transcript_name "Xkr4-203"; level 2; transcript_support_level "1"; mgi_id "MGI:3528744"; havana_gene "OTTMUSG00000026353.2"; havana_transcript "OTTMUST00000086625.1";
chr1    HAVANA  exon    3213609 3216344 .       -       .       gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000162897"; transcript_version "1"; gene_type "protein_coding"; gene_name "Xkr4"; transcript_type "lncRNA"; transcript_name "Xkr4-203"; exon_number 1; exon_id "ENSMUSE00000858910"; exon_version "1"; level 2; transcript_support_level "1"; mgi_id "MGI:3528744"; havana_gene "OTTMUSG00000026353.2"; havana_transcript "OTTMUST00000086625.1";
chr1    HAVANA  exon    3205901 3207317 .       -       .       gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000162897"; transcript_version "1"; gene_type "protein_coding"; gene_name "Xkr4"; transcript_type "lncRNA"; transcript_name "Xkr4-203"; exon_number 2; exon_id "ENSMUSE00000866652"; exon_version "1"; level 2; transcript_support_level "1"; mgi_id "MGI:3528744"; havana_gene "OTTMUSG00000026353.2"; havana_transcript "OTTMUST00000086625.1";
chr1    HAVANA  transcript      3206523 3215632 .       -       .       gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000159265"; transcript_version "1"; gene_type "protein_coding"; gene_name "Xkr4"; transcript_type "lncRNA"; transcript_name "Xkr4-202"; level 2; transcript_support_level "1"; mgi_id "MGI:3528744"; havana_gene "OTTMUSG00000026353.2"; havana_transcript "OTTMUST00000086624.1";


working dir/reference$ tail genes.gtf   # last few rows
JH584304.1      ENSEMBL start_codon     58614   58616   .       -       0       gene_id "ENSMUSG00000095041"; gene_version "7"; transcript_id "ENSMUST00000178343"; transcript_version "1"; gene_type "protein_coding"; gene_name "AC149090.1"; transcript_type "protein_coding"; transcript_name "AC149090.1-202"; exon_number 1; exon_id "ENSMUSE00001037709"; exon_version "1"; level 3; protein_id "ENSMUSP00000136649.1"; transcript_support_level "1"; tag "basic";
JH584304.1      ENSEMBL exon    56986   57151   .       -       .       gene_id "ENSMUSG00000095041"; gene_version "7"; transcript_id "ENSMUST00000178343"; transcript_version "1"; gene_type "protein_coding"; gene_name "AC149090.1"; transcript_type "protein_coding"; transcript_name "AC149090.1-202"; exon_number 2; exon_id "ENSMUSE00001053862"; exon_version "1"; level 3; protein_id "ENSMUSP00000136649.1"; transcript_support_level "1"; tag "basic";
JH584304.1      ENSEMBL CDS     56986   57151   .       -       1       gene_id "ENSMUSG00000095041"; gene_version "7"; transcript_id "ENSMUST00000178343"; transcript_version "1"; gene_type "protein_coding"; gene_name "AC149090.1"; transcript_type "protein_coding"; transcript_name "AC149090.1-202"; exon_number 2; exon_id "ENSMUSE00001053862"; exon_version "1"; level 3; protein_id "ENSMUSP00000136649.1"; transcript_support_level "1"; tag "basic";
JH584304.1      ENSEMBL exon    55112   55701   .       -       .       gene_id "ENSMUSG00000095041"; gene_version "7"; transcript_id "ENSMUST00000178343"; transcript_version "1"; gene_type "protein_coding"; gene_name "AC149090.1"; transcript_type "protein_coding"; transcript_name "AC149090.1-202"; exon_number 3; exon_id "ENSMUSE00000986146"; exon_version "1"; level 3; protein_id "ENSMUSP00000136649.1"; transcript_support_level "1"; tag "basic";
JH584304.1      ENSEMBL CDS     55483   55701   .       -       0       gene_id "ENSMUSG00000095041"; gene_version "7"; transcript_id "ENSMUST00000178343"; transcript_version "1"; gene_type "protein_coding"; gene_name "AC149090.1"; transcript_type "protein_coding"; transcript_name "AC149090.1-202"; exon_number 3; exon_id "ENSMUSE00000986146"; exon_version "1"; level 3; protein_id "ENSMUSP00000136649.1"; transcript_support_level "1"; tag "basic";
JH584304.1      ENSEMBL stop_codon      55480   55482   .       -       0       gene_id "ENSMUSG00000095041"; gene_version "7"; transcript_id "ENSMUST00000178343"; transcript_version "1"; gene_type "protein_coding"; gene_name "AC149090.1"; transcript_type "protein_coding"; transcript_name "AC149090.1-202"; exon_number 3; exon_id "ENSMUSE00000986146"; exon_version "1"; level 3; protein_id "ENSMUSP00000136649.1"; transcript_support_level "1"; tag "basic";
JH584304.1      ENSEMBL exon    52691   54867   .       -       .       gene_id "ENSMUSG00000095041"; gene_version "7"; transcript_id "ENSMUST00000178343"; transcript_version "1"; gene_type "protein_coding"; gene_name "AC149090.1"; transcript_type "protein_coding"; transcript_name "AC149090.1-202"; exon_number 4; exon_id "ENSMUSE00001045433"; exon_version "1"; level 3; protein_id "ENSMUSP00000136649.1"; transcript_support_level "1"; tag "basic";
JH584304.1      ENSEMBL UTR     58617   59690   .       -       .       gene_id "ENSMUSG00000095041"; gene_version "7"; transcript_id "ENSMUST00000178343"; transcript_version "1"; gene_type "protein_coding"; gene_name "AC149090.1"; transcript_type "protein_coding"; transcript_name "AC149090.1-202"; exon_number 1; exon_id "ENSMUSE00001037709"; exon_version "1"; level 3; protein_id "ENSMUSP00000136649.1"; transcript_support_level "1"; tag "basic";
JH584304.1      ENSEMBL UTR     55112   55482   .       -       .       gene_id "ENSMUSG00000095041"; gene_version "7"; transcript_id "ENSMUST00000178343"; transcript_version "1"; gene_type "protein_coding"; gene_name "AC149090.1"; transcript_type "protein_coding"; transcript_name "AC149090.1-202"; exon_number 3; exon_id "ENSMUSE00000986146"; exon_version "1"; level 3; protein_id "ENSMUSP00000136649.1"; transcript_support_level "1"; tag "basic";
JH584304.1      ENSEMBL UTR     52691   54867   .       -       .       gene_id "ENSMUSG00000095041"; gene_version "7"; transcript_id "ENSMUST00000178343"; transcript_version "1"; gene_type "protein_coding"; gene_name "AC149090.1"; transcript_type "protein_coding"; transcript_name "AC149090.1-202"; exon_number 4; exon_id "ENSMUSE00001045433"; exon_version "1"; level 3; protein_id "ENSMUSP00000136649.1"; transcript_support_level "1"; tag "basic";
```

Explore what the fasta file looks like.

```bash

workingdir/reference$ head genome.fa   # first few rows
>chr1 1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN


workingdir/reference$ tail -n 40 genome.fa   # last 40 rows
TGTGGTGTATCTGAAAGGTCACAATGTGGTCTATCTGAGAGGTTATGTTGTGGTCTATCT
GAGTGTTCCGGGGGTGTTCTATCTGAGAGGTCATGATGTGGTCTATCTGTAAGTCATGGT
GTGGTCTCTCTGAGAGGTCATGGTATTGTCTCTCTGAGAGGTCATGGTGTAGTCTATCAG
AGGGTTCAGGATGTGTTCTAGGCCAGAGAGCTTCTGAGGCGTCGCCATCTTCAGCTCCAG
ACACCTGGCCACCTTCCTAGTGAGAGCACAGGGGTCTGCCTGGCCCGAGAGGTTTCTGCC
TCAGG
>JH584295.1 JH584295.1
GGCTGAGCGGTGACATCATGGGCGGCGGGGTCCCAGACAGGAAGTGGGCGTGGCCTCCCA
CACTCACCCTGGCCCGCGGCGTCTGCCAGGTCGCTGTCCGAGATGCCGCCTGTGGGGGGG
GCGGGGCCGGGATGACTTCTTCGTGTTTTTTTGCACCAATTAGCATATTTTTTGGGGGGG
GGGGAGGCGCCAGGGCGATTTGCATAAAGCAGCAGGAACCAGCCGGGGCCGGGTCGGAAG
CGTGAGGGGGAGGGACAGGAAGTGACGTCACGAGTGGGCAGGGCTCACCGGGCTGTGATG
TCACTGCTCATTAGCATAGCCCCACCTCTTGGAGAATTAGCATCATAGGGGCTTCATTTG
CATGAATCCAGGACTCATTTGCATAAAGGCCACACCTCATTTGCATACAGCCGGGGCTCA
TTAGCGTAAGGCAGAGGGGGACGCCCACTTCCCAGGCCTCATTAGCCTAAATCACAGCAG
CGAAGGCTCCCCTCATTAGCATGCAGGCCCCGCCCCTCACTGCATGCAGGCCCTGCCCCT
CATTTGCATACAGGCCCCGCCCACCTGAGTCTCCGTGTGGCCTTGGCGGCTTCGGGTCCG
GCTGGGGCCGGTTTCCCGGCTTCTCTCCTGCGCCTGGGGGGCGGGACAGGAAGTGAGGCC
TCCCATGATGCCCTGGGGCCCCCGACCCCTGACCCTGCGGTGCTCACCTCCTCCGCCGCC
GCCGGGCAGGGCGTCCTCCAGGTCGAAGCCTCCTGGGGGCGGGGTCAGGGGGCGGGGTCA
GCGCAGGAGCAGGGGGCGGGGTCAGCGCAGGAGCAGGAAGCAGGGGGCGTGGCAGGAAGT
GTGGGGCGGAACTGGGAGTGCCGACAGGAAGTGCGTGGGGGGGGGAGGGGACACGTGAGT
GATGGACGCACAGGAAGTGGCCGGTAGTCCTACAGGAAGCAGGGACCCCGGCCCCGCCCA
CCAGACGTCACCTGACGGCTTAAAAGGGCGTTGGGGCTAAGGAAAGGGGATATCATGTTG
GGGTCTGGGGGGAGAAACAGGAAGTGTGGGCGGGGACAGGAAATTGCCATTAACCTTAAA
ACATTTTTCTAAGTTAATTTGCAGTGATTTGAATATTCATTTATGCAAATGAGAATGGCA
TGGCTGCTGCGCCTTGTCTGAGCATCACATCCCCGATTTGCATATTCATTTGCATATTCA
GCACCGCCGTGACATCACATGCCCCGCCCCTTCTCCCGTGAGGACCCCGCCCCCATCCCC
ACCCCATGTTGATGACATCACTCCCAGGCTCCCCAAGGGCGCCCAGAGCCCCGGGTATGT
AAATGACTCGATTTGCATATTCTCATCCGGGGAGGTGACGTCACCAGCCGGCCCTGTGTC
ACCCGTGATATCATTTCCTGTCACCGTGCGCCCCGCCCACCGCCGTGCGTGTCACTCACC
CTCCAGGGCGTCGCCCAGGTTGAAGTCGTCACTCGCTGAGGGGGGCGGGGACAAGATTAT
CATAAATATGCAAATTTAGCGCCTCATTTACATATACAGCATTTACATCCAATCCGGGCC
CCGCCCCCAGGCCCCACAGCTCGAGCTGATTTGAATATGCAGATGAGCAGACTCCGCCCC
CGGGCACTGCCTTTAACCCCGCCCATGAGCCTGACGCTGAGTTCCGGCTTCCAAGCCCCG
CCCCGCCTCATTTGCATACGGGAATCTCAGATTTAAATATTCAGGTCCTGGAGGTTAAAA
GTTGCTTTTATTGGACGGCCCCGCCCCCTGGGGGCCCCGCCCACCAGGCGAGAGGGGGTG
GGGAAGTCCGTCTGACGTCACCCCCGCGGCCCCGCCCCCAGCCCCGCCCCCTGCGCCCCC
AGGATGACATCACAGAGAACATTCAGGCCTCAAACCTAAATATTCAGATCTGAGGTTGGG
CTTCCCCAGCCCCGCCCCCTCACTTCCGGTCGGACCGCACTGGGCAGAAGCCCGTG

```
