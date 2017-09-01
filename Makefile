#!/usr/bin/make -Rrf
# Assemble Drosophila melanogaster using ABySS

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all

# Reference genome
ref=dmelanogaster

# Number of threads
t=16

# Genome size including Ns
G=143725995

# Parallel gzip with pigz
gzip=pigz -p$t

# Report run time and memory usage.
time=env time -v -o $@.time
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J

all: reads k32 k48 k64 k80 k96

sra: SRR3663859.sra SRR3663860.sra

reads: dmelanogaster.pe.fq.gz dmelanogaster.mp.fq.gz

fastqc: dmelanogaster.pe.fastqc dmelanogaster.mp.fastqc

k32: abyss/k32/dmelanogaster-scaffolds.fac.tsv \
	abyss/k32/dmelanogaster-scaffolds.bwa.samtobreak.tsv

k48: abyss/k48/dmelanogaster-scaffolds.fac.tsv \
	abyss/k48/dmelanogaster-scaffolds.bwa.samtobreak.tsv

k64: abyss/k64/dmelanogaster-scaffolds.fac.tsv \
	abyss/k64/dmelanogaster-scaffolds.bwa.samtobreak.tsv

k80: abyss/k80/dmelanogaster-scaffolds.fac.tsv \
	abyss/k80/dmelanogaster-scaffolds.bwa.samtobreak.tsv

k96: abyss/k96/dmelanogaster-scaffolds.fac.tsv \
	abyss/k96/dmelanogaster-scaffolds.bwa.samtobreak.tsv

nxtrim: dmelanogaster.mp.nxtrim.fq.gz

nxtrim_k32: nxtrim/abyss/k32/dmelanogaster-scaffolds.fac.tsv \
	nxtrim/abyss/k32/dmelanogaster-scaffolds.bwa.samtobreak.tsv

ifndef k
abyss/k%/dmelanogaster-scaffolds.fac.tsv:
	mkdir -p $(@D)
	$(time) $(MAKE) k=$* $@ 2>&1 | tee $@.log

abyss/k%/dmelanogaster-scaffolds.bwa.samtobreak.tsv:
	$(time) $(MAKE) k=$* $@ 2>&1 | tee $@.log

nxtrim/abyss/k%/dmelanogaster-scaffolds.fa:
	mkdir -p $(@D)
	$(time) $(MAKE) k=$* $@ 2>&1 | tee $@.log
endif

# Download data from Ensembl.
dmelanogaster.fa:
	curl ftp://ftp.ensembl.org/pub/release-90/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.gz \
		| gunzip -c | seqtk seq >$@

# Download the reference genome.
SRR366%.sra:
	curl -o $@ ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR366/SRR366$*/$@

# sratoolkit

# Convert SRA to FASTQ format with fastq-dump.
SRR%.fq.gz: SRR%.sra
	fastq-dump -Z --split-spot $< | $(gzip) >$@

# Download the FASTQ data with fastq-dump.
SRR%.fq.gz:
	fastq-dump -Z --split-spot SRR$* | $(gzip) >$@

# Symlink the paired-end data.
dmelanogaster.pe.fq.gz: SRR3663859.fq.gz
	ln -sf $< $@

# Symlink the mate-pair data.
dmelanogaster.mp.fq.gz: SRR3663860.fq.gz
	ln -sf $< $@

# FastQC

# Inspect the quality of the reads using FastQC.
%.fastqc.html: %.fq.gz
	fastqc -t $t $<
	mv $*_fastqc.html $*.fastqc.html
	mv $*_fastqc.zip $*.fastqc.zip

# NxTrim

# Trim mate-pair reads using NxTrim.
%.nxtrim.fq.gz: %.fq.gz
	nxtrim --stdout --justmp --rf -1 <(seqtk seq -1 $<) -2 <(seqtk seq -2 $<) | $(gzip) >$@

# Symlink the paired-end reads.
nxtrim/%.pe.fq.gz: %.pe.fq.gz
	mkdir -p $(@D)
	ln -sf ../$< $@

# Symlink the trimmed mate-pair reads.
nxtrim/%.mp.fq.gz: %.mp.nxtrim.fq.gz
	mkdir -p $(@D)
	ln -sf ../$< $@

# samtools

# Index a FASTA file.
%.fa.fai: %.fa
	samtools faidx $<

# Sort a SAM file and produce a sorted BAM file.
%.sort.bam: %.sam
	samtools sort -@$t -o $@ $<

# Index a BAM file.
%.bam.bai: %.bam
	samtools index $<

# BWA

# Index the target genome.
%.fa.bwt: %.fa
	bwa index $<

# Align sequences to the target genome.
%.bwa.sam: %.fa $(ref).fa.bwt
	bwa mem -t$t -xintractg $(ref).fa $< >$@

# ABySS

# Assemble paired-end and mate-pair reads using ABySS.
abyss/k$k/%-scaffolds.fa: %.pe.fq.gz %.mp.fq.gz
	test ! -e $@
	mkdir -p $(@D)
	$(time) abyss-pe -C $(@D) mpirun=mpirun np=$t G=$G v=-v name=$* k=$k lib=pe1 mp=mp1 pe1=../../$*.pe.fq.gz mp1=../../$*.mp.fq.gz 2>&1 | tee $@.log

# Assemble paired-end and trimmed mate-pair reads using ABySS.
%/abyss/k$k/dmelanogaster-scaffolds.fa: %/dmelanogaster.pe.fq.gz %/dmelanogaster.mp.fq.gz
	test ! -e $@
	mkdir -p $(@D)
	$(time) abyss-pe -C $(@D) mpirun=mpirun np=$t G=$G v=-v name=dmelanogaster k=$k lib=pe1 mp=mp1 pe1=../../dmelanogaster.pe.fq.gz mp1=../../dmelanogaster.mp.fq.gz 2>&1 | tee $@.log

# Calculate assembly contiguity stats using abyss-fac.
%.fac.tsv: %.fa
	abyss-fac -G$G -t1000 $< >$@

# Calculate assembly contiguity and correctness metrics with abyss-samtobreak.
%.samtobreak.txt: %.sam
	(echo '==> $< <=='; abyss-samtobreak -G$G -l1000 $<) >$@

# Convert samtobreak.txt to TSV.
%.samtobreak.tsv: %.samtobreak.txt
	abyss-samtobreak-to-tsv $< >$@
