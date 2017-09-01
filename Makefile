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

all: reads nxtrim \
	k32 k48 k64 k80 k96 \
	nxtrim-k32 nxtrim-k48 nxtrim-k64 nxtrim-k80 nxtrim-k96 \
	notebook

sra: SRR3663859.sra SRR3663860.sra

reads: dmelanogaster.pe.fq.gz dmelanogaster.mp.fq.gz

k32 k48 k64 k80 k96: k%: \
	abyss/k%/dmelanogaster.scaffolds.fac.tsv \
	abyss/k%/dmelanogaster.scaftigs.fac.tsv \
	abyss/k%/dmelanogaster.scaftigs.bwa.samtobreak.tsv

nxtrim: dmelanogaster.mp.nxtrim.fq.gz

nxtrim-k32 nxtrim-k48 nxtrim-k64 nxtrim-k80 nxtrim-k96: nxtrim-k%: \
	nxtrim/abyss/k%/dmelanogaster.scaffolds.fac.tsv \
	nxtrim/abyss/k%/dmelanogaster.scaftigs.fac.tsv \
	nxtrim/abyss/k%/dmelanogaster.scaftigs.bwa.samtobreak.tsv

%.samtobreak.tsv: \
		abyss/k32/%.scaftigs.bwa.samtobreak.tsv \
		abyss/k48/%.scaftigs.bwa.samtobreak.tsv \
		abyss/k64/%.scaftigs.bwa.samtobreak.tsv \
		abyss/k80/%.scaftigs.bwa.samtobreak.tsv \
		abyss/k96/%.scaftigs.bwa.samtobreak.tsv \
		nxtrim/abyss/k32/%.scaftigs.bwa.samtobreak.tsv \
		nxtrim/abyss/k48/%.scaftigs.bwa.samtobreak.tsv \
		nxtrim/abyss/k64/%.scaftigs.bwa.samtobreak.tsv \
		nxtrim/abyss/k80/%.scaftigs.bwa.samtobreak.tsv \
		nxtrim/abyss/k96/%.scaftigs.bwa.samtobreak.tsv
	mlr --tsvlite cat $^ >$@

fastqc: \
	dmelanogaster.pe.fastqc.html \
	dmelanogaster.mp.fastqc.html \
	dmelanogaster.mp.nxtrim.fastqc.html

notebook: dmelanogaster.samtobreak.nb.html

ifndef k
abyss/k%/dmelanogaster.scaffolds.fa:
	mkdir -p $(@D)
	$(time) $(MAKE) k=$* $@ 2>&1 | tee $@.log

nxtrim/abyss/k%/dmelanogaster.scaffolds.fa:
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
%.bwa.sam.gz: %.fa $(ref).fa.bwt
	bwa mem -t$t -xintractg $(ref).fa $< | $(gzip) >$@

# seqtk

# Break scaffolds into scaftigs using seqtk.
%.scaftigs.fa: %.scaffolds.fa
	seqtk cutN -n1 $< | seqtk seq >$@

# ABySS

# Assemble paired-end and mate-pair reads using ABySS.
abyss/k$k/%-scaffolds.fa: %.pe.fq.gz %.mp.fq.gz
	test ! -e $@
	mkdir -p $(@D)
	$(time) abyss-pe -C $(@D) mpirun=mpirun np=$t G=$G v=-v name=$* k=$k lib=pe1 mp=mp1 mp1_de=--rf pe1=../../$*.pe.fq.gz mp1=../../$*.mp.fq.gz 2>&1 | tee $@.log

# Assemble paired-end and trimmed mate-pair reads using ABySS.
%/abyss/k$k/dmelanogaster-scaffolds.fa: %/dmelanogaster.pe.fq.gz %/dmelanogaster.mp.fq.gz
	test ! -e $@
	mkdir -p $(@D)
	$(time) abyss-pe -C $(@D) mpirun=mpirun np=$t G=$G v=-v name=dmelanogaster k=$k lib=pe1 mp=mp1 pe1=../../dmelanogaster.pe.fq.gz mp1=../../dmelanogaster.mp.fq.gz 2>&1 | tee $@.log

# Symlink .scaffolds.fa
%.scaffolds.fa: %-scaffolds.fa
	ln -sf $(<F) $@

# Calculate assembly contiguity stats using abyss-fac.
%.fac.tsv: %.fa
	abyss-fac -G$G -t1000 $< >$@

# Calculate assembly contiguity and correctness metrics with abyss-samtobreak.
%.samtobreak.txt: %.sam.gz
	(echo "File: $<"; gunzip -c $< | abyss-samtobreak -G$G -l1000) >$@

# Convert samtobreak.txt to TSV using Miller.
%.samtobreak.tsv: %.samtobreak.txt
	mlr --ixtab --ips ': ' --otsvlite --from $< \
		then rename 'Number of unmapped contigs,Unmapped_contigs' \
		then rename 'Total length of unmapped contigs,Unmapped_bases' \
		then rename 'Mapped contig bases,Mapped_bases' \
		then rename 'Mapped NG50,Contig_NGA50' \
		then rename 'Number of Q10 break points longer than 500 bp,Contig_breakpoints' \
		then rename 'Scaffold NG50,Scaffold_NG50' \
		then rename 'Aligned scaffold NG50,Scaffold_NGA50' \
		then rename 'Number of Q10 scaffold breakpoints longer than 500 bp,Scaffold_breakpoints' \
		then cut -r -x -f ' ' \
		then put '$$Total_breakpoints = $$Contig_breakpoints + $$Scaffold_breakpoints' \
		>$@

# RMarkdown

# Generate a report of assembly metrics using RMarkdown.
%.samtobreak.nb.html: %.samtobreak.tsv assembly-metrics.rmd
	Rscript -e 'rmarkdown::render("assembly-metrics.rmd", "html_notebook", "$*.samtobreak.nb.html", params = list(input_tsv="$<"))'
