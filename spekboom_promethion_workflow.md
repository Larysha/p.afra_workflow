
## Spekboom Genome Assembly - PromethIon Pipeline

## Overview and Bioinformatics Methods
In Sep 2023, two new data sets for the SA Spekboom accession 28718 EvJ were produced on the PromethION 24 with the SQK-LSK114 kit, using the optimised HMW DNA extraction protocol that was developed in the first phase. The first set of unfragmented reads produced ~29.32 gigabases that passed the quality check and an estimated N50 of 18.99kb, whereas the reads that underwent random fragmentation resulted in better pore activity during sequencing with 46.24 gigabases yield and estimated 22.81 N50 (based on initial fast basecalling with `Guppy`).

Fragmented - re-loaded flow cell after wash with the same library. 

## Package version
The following tools and versions are required for this analysis and can be loaded as one .sh file:

```bash
#!/bin/bash -e
module purge

module load \
chpc/BIOMODULES \
nanofilt \
nanoplot \
flye/2.9 \
minimap2/2.24 \
samtools/1.9 \
medaka/1.11.3 \
quast/5.2.0 \
busco/5.6.1 \
blast \
blobtools/1.1.1 \
jellyfish/2.2.10 \
genomescope-2 
```

All analysis was run on the CHPC (Centre for High Performance Computing through the CSIR, South Africa).

### SUP basecalling with Dorado
After data upload to the CHPC (~ 3 days) The CHPC Nvidia GPU was used to complete basecalling in 24 hours (for each set of `pod5` files), however `Dorado` can be installed and runs locally.

**To set up dorado with WSL (windows subsystem for linux):**
Visit the `nanoporetech` dorado GitHub page https://github.com/nanoporetech/dorado#installation and select the Linux version of `dorado`.

```shell
# decompress the tar file with
tar -xzvf archive.tar.gz
# you can download specific modules as follows:
dorado download --model dna_r9.4.1_e8_sup@v3.6
```

**Basecalling:**
The variable file paths are all defined in the 0.config.sh file and module `dna_r10.4.1_e8.2_400bps_sup@v4.2.0` was used for the `SUP` assembly

```shell
#!/bin/bash -e
dorado basecaller ${model} ${raw_pod5} > ${basecalled_dir}sup.bam
# the --resume-from flag can be used if basecalling does not complete or gets interupted
```

Basecalling output was two unmapped `BAM` files - `sup-unfragmented.bam` and `sup-fragmented.bam`. Since most QC tools require `fastq` files, these were converted with `Samtools`:

```shell
samtools bam2fq ${bam_file}.bam > ${output_file}.fastq
# which can be compressed with:
gzip -c ${fastq_file}.fastq > ${output_file}.fastq.gz
```

## Read QC with NanoPlot and Chopper
To get an idea of the read quality before and after filtering, QC was run before and after `Chopper`:
```shell
NanoPlot -t 8 --color blue --fastq ${basecalled}sup_reads.fastq.gz -o ${read_QC}nanoplot
```

### Unfragmented QC:

| metric | raw stats | filtered stats|
| -- | -- | -- |
|n | 7 198 152| 4 919 667 |
| N50 | 19 718 | 20 724|
|mean length | 4 080.8| 5 345.7|
|mean quality | 16.9 | 18.5|
|max  | 341 637 | 341 637 |
| sum bases | 29 374 234 250 | 26 298 989 984 |

### Fragmented QC

| metric | raw stats | filtered stats|
| -- | -- | -- |
|n | 10 599 346| 6 781 778 |
| N50 | 23 369 | 23 972 |
|mean length |4 370.9| 6 117.1 |
|mean quality | 16.6 | 18.4 |
|max  | 327 351 | 327 351  |
| sum bases | 46 328 730 327 | 41 485 041 224 |


From these results, chopper was used to retain all read with Q15 or greater, which translates to an error rate between 3.2% - 7.1% (69.3% of the unfragmented and 65.1% of the fragmented reads retained). Based on the `NanoPlot` graph, we retained all our longest reads and a min read length was set to 1000 bps as follows:

```bash
zcat \
    ${raw_reads} |\
    NanoFilt -l 1000 -q 15 --headcrop 10 --tailcrop 10 \
    > 'pri_filt.fastq'
```
## Jellyfish and Genome Scope

K-mer analysis can give an estimate of the ploidy, heterozygosity and haploid genome length. The k-mer size is set to 21 and since the Q score threshold has been set to 15 (error rate of ~ 2-3%), which is less than an error every 30 bases (greater than the k-mer size). For this reason, it was considered appropriate to use jellyfish on the PromethION data, even though high quality short read data is recommended.

```bash
jellyfish count \
    -C -m 21 -t 10 \
    -s ${params.approx_size} \
    $filtered_reads \
    -o reads-kmers.jf

jellyfish histo \
    -t 10 reads-kmers.jf \
    > reads-kmers.histo

genomescope.R \
    -i $jelly \
    -o genomescope \
    -k 21

```
We used MGI and Illumina WGS data to perform the k-mer analysis, with the following GenomeScope2 output:

GenomeScope version 2.0 \
input file = reads-kmers.histo \
output directory = genomescope2\
p = 4
k = 21

|property     |                 min         |      max |
|--|--|--|
|Homozygous (aaaa)            |  96.9606%      |     97.8931%
|Heterozygous (not aaaa)      |  2.10693%      |     3.03942%
|aaab       |   1.68976%      |  1.72239%
|aabb       |   0.417176%     |  0.815141%
|aabc       |  0%             |  0.35173%
|abcd       | 0%              |  0.150158%
|Genome Haploid Length        |  283,159,051 bp   |  283,476,057 bp
|Genome Repeat Length         |  70,274,595 bp    |  70,353,269 bp
|Genome Unique Length         |  212,884,456 bp   |  213,122,788 bp
|Model Fit                    |  78.4124%         |  95.7216%
|Read Error Rate              |  0.429873%        |  0.429873%



## Genome assembly with Flye 2.9
genome assembly was run using both unfragmented and fragmented reads combined: this assembly was built with the following parameters, with the `--keep-haplotypes` flag. Since plant genome assemblies are typically complex with the potential for ploidy and highly repetitive elements - a collapsed approach can lead to missassbly of these difficult to resolve features. The "un-collpsed" approach leads to a less contiguous assembly but with less chance of misassembly - these "duplicated regions" can then be removed in post-processing steps.

```bash
flye \
--nano-hq ${filt_reads}*.fastq \
-g 0.8G \
-o ${assembly} \
--read-error 0.03 \
--asm-coverage 40 \
--threads 56 \
--keep-haplotypes
```

without the `--keep-haplotigs` assembly stats were as follows

|metric|measure|
|--|--|
|.fasta file size| 425MB|
| total length (bps)| 437849582|
| no. of fragments| 4569|
|fragment N50 | 271538|
|largest fragment| 2839499|
|scaffolds| 0|
| mean coverage  | 139|

With this flag (uncollapsed assembly), the stats were:

|metric|measure|
|--|--|
|.fasta file size| 442MB|
| total length (bps)| 455517192|
| no. of fragments| 8456|
|fragment N50 | 139504|
|largest fragment| 3464806|
|scaffolds| 0|
| mean coverage  | 127|

## Long read polishing with Medaka - to be edited 
`Racon` is a tool that corrects raw contigs from rapid assembly (typically the output of `Flye`) to generate a genomic consensus and error correction (especially useful for indel errors in the assembly). `Medaka`, which is the ONT tool that produces the final consensus, requires an input assembly that is pre-processed with `Racon`. `Racon`, requires a `SAM` input file that contains the overlap between the raw, filtered reads mapped back to the draft assembly. To achieve this, we first use `minimap2`:

At the time the assembly was run, `Medaka` had not yet released a module specific to the 10.4 chemistry or `SUP` basecalling with `Dorado`. For this reason, the next highest module was chosen for PromethION data - `r103_prom_high_g360`.

```bash
medaka_consensus \
${chopped_reads}*.fastq \
-d ${haploid}assembly.fasta \
-o ${polish} \
-t 52 \
-m r103_prom_high_g360
```

## Assembly QC with QUAST, BUSCO and BlobTools

### QUAST
`QUAST` gives an evaulation of all key assembly metrics to assess the output:

```bash
# run quast
python3 \
${quast}quast.py \
-o ${out}quast \
--gene-finding \
--eukaryote \
${assembly}assembly.fasta
```

### BUSCO
a widely used metric based on evolutionarily-informed expectations of gene content of near-universal single-copy orthologs (core genes in a particular lineage). A measure of the complete, core gene content that is represented in the assembly. 

```bash
# run busco with the viridiplantae lineage
busco \
-i ${assembly}assembly.fasta \
-o busco_viridiplantae \
-l ${lineages}viridiplantae_odb10 \
-m genome \
--offline
```

## Assembly Contamination and Small Contigs Filtering

### BLAST and BlobTools
`BlobTools` provides a modular solution "for visualisation, quality control and taxonomic partitioning of genome datasets" and requires a `BAM` file of the raw reads mapped to the assembly, as well as the assembly contigs BLASTED against a nucleotide database as input:

Start by aligning reads back to the genome with minimap2
```bash
minimap2 \
-ax map-ont \
${assembly}assembly.fasta \
${reads}*.fastq | \
samtools sort > aln.bam

samtools index aln.bam
```

BLAST-ing a whole assembly can be a resource intensive exercise - with 400GB RAM and 48 hours walltime, the `BLAST` still does not complete (~60% of the assembly only). So, `BLAST` is run in parallel, by splitting up the assembly contigs and then stitching the results back together, as follows:

```shell
# to divide the assembly file by contigs, first count the contigs:
grep -o '^>[^ ]*' ${assem}pri_asm.fasta | cut -c 2- > ${out}contigs.txt

# count the number of contigs
contig_count=$(wc -l < ${out}contigs.txt)
echo "The number of contigs in contigs.txt is: $contig_count" > ${out}num_contig.txt

# split the file into parts (I have chosen 12 - adjust as needed)
split -l $((contig_count/12)) -d -a 2 ${out}contigs.txt ${out}contigs_

# seqtk subseq command is used to extract sequences from a FASTA file based on a list of sequence names
for file in ${out}contigs_*; do
    index=$(echo "$file" | cut -d'_' -f2)  # Extract the numeric part from the filename / index
    seqtk subseq ${assem}pri_asm.fasta $file > ${out}contigs_${index}.fasta
done
```

`BLAST` to NCBI nucleotide DB - repeat for each sub-file

```shell
# find the blast database in the CHPC environment
which blastn
env | grep DB

blastn \
-db nt \
-query ${assembly}contigs_00.fasta \
-max_target_seqs 10 \
-max_hsps 1 \
-evalue 1e-25 \
-outfmt "6 qseqid staxids bitscore std" \
-out ${out}blast-nt_00.out 
```

The stitch the BLAST outputs back together again with:
```shell
cat ${blast}blast-nt_*.out > ${blast}blast-nt.out
```

Now, BlobTools can run with:

```shell
blobtools create \
-i ${assembly}pri_asm.fasta \
-b ${map}aln.bam \
-t ${blast}blast-nt.out \
-o ${out}contam \
--nodes ${TAXNODES}nodes.dmp \
--names ${TAXNODES}names.dmp

blobtools plot -i contam.blobDB.json -o blobP
blobtools view -i contam.blobDB.json -o blobV
```
Once BlobTools has run, any contaminants were removed from the assembly with:

```shell
grep -v '^#' ${blob}blobout.contam.blobDB.table.txt |\
grep -vE 'Chlorophyta|Artverviricota' |\
awk '{print $1}' > keep.lst

grep -v '^#' ${blob}blobout.contam.blobDB.table.txt |\
grep -E 'Chlorophyta|Artverviricota' |\
awk '{print $1}' > contaminant.lst


seqtk subseq ${genome}assembly.fasta keep.lst > pri_asm.fasta
seqtk subseq ${genome}assembly.fasta contaminant.lst > contamination.fasta
```

In this case only 2 small contaminant sequences were found from the taxanomic groups Chlorophyta (green algae) and Artverviricota (small RT-virus).



## Assembly post-processing: Purge Haplotigs
When regions of the parental haplotypes differ the assembler will produce two sequences that represent the same genetic region. As we are building a pseudo-haploid genome we need to identify the sequences that make up these doubly represented regions and remove one. These duplicates are particularly likely in a polyploid genome, where different versions of the same gene may be represented in the assembly many times, especially if there is high heterozygosity (also seen in higher duplicated `BUSCO` scores). Haplotig removal is performed to improve genome contiguity, as such only the longest sequence representative of a region is kept. We can also potentially identify plasmid sequences in our assembly as they will have a very high or very low read coverage.


```bash
minimap2 -x map-ont ${data}pri_filt.fastq ${assembly}pri_asm.fasta |\
 gzip -c - > ${out}asm_aln.paf.gz
 
# generate BP.base.cov and PB.stat files
${exe}pbcstat ${out}*.paf.gz
${exe}calcuts ${out}PB.stat > cutoffs 2>${out}calcults.log

# alternatively, you can set your cutoff manually by observing the coverage histogram
python3 ${purge_dups}scripts/hist_plot.py -c cutoffs PB.stat PB.cov.png
# and then assigning the lower, middle and upped thresholds as follows:
calcuts -l 5 -m 85 -u 200 PB.stat > cutoffs_manual

# create split fasta file
${exe}split_fa ${assembly}polished.fasta > ${out}asm.split.fasta
# map back to it
minimap2 -x map-ont -DP ${out}asm.split.fasta ${out}asm.split.fasta |\
 qzip -c - > ${out}asm_split_self.paf.gz

# purge the duplicates based on manual cutoffs
${exe}purge_dups -2 -T cutoffs_manual -c PB.base.cov ${out}asm_split_self.paf.gz \
> dups.bed 2> purge_dups.log
${exe}get_seqs 
```

## Short read polishing with NextPolish and POLCA: A comparitive analysis

We have access to both WGS MGI as well as Illumina short read sequencing data for our accesssion. These two technologies were used to generate the high quality (Q> 30) reads for assembly polishing, which were benchmarked based on two polishing tools.

`POLCA` (POLishing by Calling Alternatives) records the number of changes made between the final polished and the original assembly in a `VCF` file, making it a useful tool for a comparison study. Since `POLCA` polishes the assembly by aligning the short reads to the genome and calling variants - this `VCF` file contains the “variants” that were considered errors and subsequently polished. A variant is considered an error “if the count of the alternative allele observations is greater than 1 and at least twice the count of the reference allele.” The error variant is replaced with the highest scoring alternative allele from the Illumina reads.

```shell
module load \
chpc/BIOMODULES \
MaSuRCA/4.0.1

polca.sh \
-a ${genome} \
-r ${reads} \ # reads is a text files listing the full paths to the short reads
-t 56
```

`NextPolish` unfortunately does not output a file of the errors that were corrected - but is one of the best polishing tools available. It’s output was measured by aligning the polished assembly back to the original assembly and measuring the differences, with `DNAdiff`

```shell
nextPolish ${output}run.cfg

# where run.cfg contains:
task = best
genome = /mnt/lustre/users/lrothmann/spekboom/assembly/flye_asm.fasta
sgs_fofn = /mnt/lustre/users/lrothmann/mgi/polish/sgs.fofn
workdir = ./01_nextpolish
rewrite = yes
deltmp = yes
rerun = 4

# and sgs.fofn contains:
/mnt/lustre/users/lrothmann/MGI/raw_reads/V350190840_L01_13_1.fq.gz
/mnt/lustre/users/lrothmann/MGI/raw_reads/V350190840_L01_13_2.fq.gz
/mnt/lustre/users/lrothmann/MGI/raw_reads/V350190840_L01_14_1.fq.gz
/mnt/lustre/users/lrothmann/MGI/raw_reads/V350190840_L01_14_2.fq.gz
/mnt/lustre/users/lrothmann/MGI/raw_reads/V350190840_L01_15_1.fq.gz
/mnt/lustre/users/lrothmann/MGI/raw_reads/V350190840_L01_15_2.fq.gz
/mnt/lustre/users/lrothmann/MGI/raw_reads/V350190840_L01_16_1.fq.gz
/mnt/lustre/users/lrothmann/MGI/raw_reads/V350190840_L01_16_2.fq.gz
```

`DNAdiff` can then be run with:
```shell
mummer-4_PR30

dnadiff \
${data}flye_asm.fasta \
${data}nextpolished_mgi.fasta
```

## Measure the collapse of the genome with HapPy
This tool assesses the 'haploidy' _H_ of a given assembly. _H_ is defined as the fraction of the bases of the genome that are in the collapsed peak _C_. This metrics is calculated as _H_=_C_/(_C_+_U_/2), where _C_ is the size (area) of the collapsed peak and _U_ the size of the uncollapsed peak in the per-base coverage histogram of the assembly. 

***Here 'haploidy' is defined as the measure of properly collapsed haplotypes (the aim of Purge Haplotigs) - for Spekboom this may be less than n if a lot of collapse takes place*** 


## Repeat annotation with EDTA and RepeatMasker
https://github.com/oushujun/EDTA 
https://github.com/rmhubley/RepeatMasker

## Gene annotation with Augustus 
https://github.com/Gaius-Augustus/Augustus