# phasihunter
Welcome to phasihunter 😉

A multithreaded program for mining phasiRNA regulation pathways based on multiple reference sequences.

![image.png|425](https://sandbox-1314381151.cos.ap-nanjing.myqcloud.com/pic/202307312031023.png)

## Table of contents
- Dependencies 
- Installation
- Use case
- Copyright

## Dependencies
phasihunter is a CLI program runing on linux platform. The correction runing of phasihunter depends on some existing softwares.
- Bowtie (Langmead, et al., 2009. Genome Biol)
- Biopython (Cock, et al., 2009. Bioinformatics)
- Bedtools (Quinlan and Hall, 2010. Bioinformatics)
- Dnapi (Tsuji and Weng, 2016. PloS One)
- Trim_galore (https://github.com/FelixKrueger/TrimGalore)
- Seqkit (Shen, et al., 2016. PloS One)
- Perl5 (https://www.perl.org)
- Fasta36 (Pearson and Lipman, 1988. Proc Natl Acad Sci U S A)
- TarHunter (Ma, et al., 2018. Bioinformatics)

# Installation
1. Clone phasihunter
`git clone https://github.com/HuangLab-CBI/PhasiHunter.git .`

2. Setting enviroment variable in ~/.bashrc
`export PATH=$PATH:<phasihunter PATH>`

Now type `phasihunter -h` to check phasihunter whether installation correct.

# Use case
Parameter in < > means necessary; parameter in [ ] means optional

1. Data pre-process
```bash
phasiHunter preprocess -m <r> -i <SRR5049781.fastq.gz> -r <oryza_sativa_cdna.fa> -o [SRR5049781_cdna.map]

phasiHunter preprocess -m <r> -i <SRR5049781.fastq.gz> -r <oryza_sativa_gdna.fa> -o [SRR5049781_gdna.map]
```
2. PhasiRNA and PHAS Loci prediction
```bash 
phasiHunter phase -cm <SRR5049781_cdna.map> -c <oryza_sativa_cdna.fa> -gm <SRR5049781_gdna.map> -g <oryza_sativa_gdna.fa> -fa <SRR7851621_trimmed_format_filter.fa> -a [SRR5049781_allsiRNA.txt] -o [SRR5049781_phasiRNA.txt] -pl [21] -j [10] -pv [0.0001] -ps [15] -pr [0.4] 
```
3. PhasiRNA and PHAS Loci result integration
```bash
phasiHunter integration -io <SRR5049781_phasiRNA.txt> -ia <SRR5049781_allsiRNA.txt> -an <oryza_sativa_gdna.gff3> -o [SRR5049781_phasiRNA_dup.txt] -a [SRR5049781_allsiRNA_dup.txt] -s [SRR5049781_summary.txt] -po [SRR5049781_phas.txt] -g <y>
```
4. Print phasiRNA_cluster plot, phasiRNA.fa, PHAS.fa
```bash
phasiHunter visulization -io <SRR5049781_phasiRNA_dup.txt> -ia <SRR5049781_allsiRNA_dup.txt> -ip <SRR5049781_phas.txt> -a [SRR5049781_alignment.txt] -o [SRR5049781.phasiRNA.fa] -p [SRR5049781.PHAS.fa] -c [oryza_sativa_cdna.fa] -g [oryza_sativa_gdna.fa] -pc [y] -pg [y]
```
5. Initiator prediction and verification
```bash 
phasiHunter target -q <osa_miRNA.fa> -b <SRR5049781_PHAS.fa> -o <SRR5049781_miR.txt> -t

phasiHunter initiator -i <SRR5049781_phasiRNA_dup.txt> -j <SRR5049781_miR.txt> -ip <SRR5049781_phas.txt> -o <SRR5049781_initiator.txt>

phasiHunter deg -i <degradome_PHAS.map> -q <osa_miRNA.fa> -j <SRR5049781_initiator.txt> -t <SRR5049781_PHAS.fa> -o <SRR5049781_initiator_verified.txt> -in <y>
```

6. PhasiRNA target prediction and verification
```bash
phasiHunter target -q <SRR5049781_phasiRNA.fa> -t <oryza_sativa_cdna.fa> -o <SRR5049781_phasiRNA_target.txt>

phasiHunter deg -i <degradome_cdna.map> -q <SRR5049781_phasiRNA.fa> -j <SRR5049781_phasiRNA_target.txt> -t <oryza_sativa_cdna.fa> -o <SRR5049781_phasiRNA_target_verified.txt> -in <n>
```

# Copyright
Copyright © Crop Bioinformatics Group (CBI), College of Agricultural, Nanjing Agricultural University.

Free for academic use. For commercial use, please contact us (huangji@njau.edu.cn)

