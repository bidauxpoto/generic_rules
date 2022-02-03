NGS_LIBRARY_PROPERTY_STRANDED ?= __NGS_LIBRARY_PROPERTY_STRANDED
# 1 -> forward
# 2 -> reverse
# 0 -> unstranded

NGS_LIBRARY_PROPERTY_PAIRED ?= N
# Y|N

######################
#
# System configuration
#

SCRATCH?=__SCRATCH__
TMP_DIR?=$(SCRATCH)/tmp/
NCPUS?=8

######################
#
#   Reference data location
#

REFERENCE_GENOME_FASTA?=__REFERENCE_GENOME_FASTA__
REFERENCE_FEATURECOUNTS_GTF?=__REFERENCE_FEATURECOUNTS_GTF__

BLASTDB?=__BLASTDB__ 
# directory containing ncbi blast formatted sequence database, see https://www.ncbi.nlm.nih.gov/books/NBK62345/ and https://ftp.ncbi.nlm.nih.gov/blast/db/

REFERENCE_RSEQC_HOUSEKEEPING_BED?=__REFERENCE_RSEQC_HOUSEKEEPING_BED__ 
# https://sourceforge.net/projects/rseqc/files/BED/

RIBOSOMAL_BED_FILE?=__RIBOSOMAL_BED_FILE__
######################
#
#   Parameters
#


FEATURECOUNTS_PARAM ?= -t exon -g gene_name -s $(NGS_LIBRARY_PROPERTY_STRANDED)
TRIM_GALORE_PARAM ?= --stringency=3
BIGWIG_BIN_SIZE?=10

######################
#
#    Deduplication
#

%.mark_dup.bam: %.bam
	picard MarkDuplicates -I $< -M $@.log -O $@

%.no_dup.bam: %.mark_dup.bam
	samtools view -F 1024 $< > $@

%.umi_dedup.bam: %.bam %.bam.bai
	umi_tools dedup --extract-umi-method read_id --umi-separator $(UMI_SEPARATOR) -I $< -S $@ --output-stats $@ --method unique


#######################
#
#   Sequence Filtering
#

%.dust.fa.gz: %.fa.gz
	zcat $< | dustmasker -outfmt fasta | gzip > $@

#######################
#
#   Format Conversion
#

%.bedGraph: %.bw
	bigWigToBedGraph $< $@

%.bed: %.bam
	bedtools bamtobed -splitD < $< | bsort -k1,1V -k2,2n > $@
  
%.$(BIGWIG_BIN_SIZE).bw: %.bam %.bam.bai
	bamCoverage --binSize=$(BIGWIG_BIN_SIZE) -b $< -o $@ --numberOfProcessors=$(NCPUS)
%.bw: %.bam %.bam.bai
	bamCoverage --binSize=$(BIGWIG_BIN_SIZE) -b $< -o $@ --numberOfProcessors=$(NCPUS)
%.norm.$(BIGWIG_BIN_SIZE).bw: %.bam %.bam.bai
	bamCoverage --binSize=$(BIGWIG_BIN_SIZE) --normalizeUsing=CPM -b $< -o $@ --numberOfProcessors=$(NCPUS)
%.norm.bw: %.bam %.bam.bai
	bamCoverage --binSize=$(BIGWIG_BIN_SIZE) --normalizeUsing=CPM -b $< -o $@ --numberOfProcessors=$(NCPUS)

%.fa.gz: %.fastq.gz
	zcat $< | fastq2tab | enumerate_rows | cut -f 1,3 | tab2fasta -s | gzip >$@
 
 %.cram: %.bam
	samtools view -@ $(NCPUS) -T $(REFERENCE_GENOME_FASTA) -C -o $@ $<

######################
#
#   RNA-seq
#

ifeq ($(NGS_LIBRARY_PROPERTY_PAIRED),N)

%.trimmed.fastq.gz: %.fastq.gz
  docker run -u `id -u`:$(DOCKER_GRP) --rm -v $(BIOINFO_ROOT):$(BIOINFO_ROOT) -v $(PRJ_ROOT):$(PRJ_ROOT) -v $(SCRATCH):$(SCRATCH) quay.io/biocontainers/trim-galore:0.6.5--0 bash -c "cd $(PWD); trim_galore -j $(NCPUS) $< -o `dirname $@` $(TRIM_GALORE_PARAM); mv `dirname $@`/$*_trimmed.fq.gz `dirname $@`/$*.trimmed.fastq.gz";\

else

%_R1.trimmed.fastq.gz: %_R1.fastq.gz %_R2.fastq.gz
	CMD="\
		cd $(PWD); trim_galore --paired -j $(NCPUS) $^ -o `dirname $@` $(TRIM_GALORE_PARAM);\
		mv $*_R1_val_1.fq.gz $*_R1.trimmed.fastq.gz;\
		mv $*_R2_val_2.fq.gz $*_R2.trimmed.fastq.gz\
	";\
	docker run -u `id -u`:$(DOCKER_GRP) --rm -v $(BIOINFO_ROOT):$(BIOINFO_ROOT) -v $(PRJ_ROOT):$(PRJ_ROOT) -v $(SCRATCH):$(SCRATCH) quay.io/biocontainers/trim-galore:0.6.5--0 bash -c "$$CMD";\

%_R2.trimmed.fastq.gz: %_R1.trimmed.fastq.gz
	@echo alreay done

endif


%_fastqc.html: %.fastq.gz
	docker run -u `id -u`:`id -g` --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) biocontainers/fastqc:v0.11.5_cv3  bash -c "cd $(PWD); fastqc -o `dirname $@` -t $(NCPUS) $<"

%.trimmed.fastq.gz: %.fastq.gz
	docker run -u `id -u`:`id -g` --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(PRJ_ROOT):$(PRJ_ROOT) quay.io/biocontainers/trim-galore:0.6.5--0  bash -c "cd $(PWD); trim_galore -j 6 $< -o . $(TRIM_GALORE_PARAM); mv $*_trimmed.fq.gz $@"


%.bam.featurecounts.count: $(REFERENCE_FEATURECOUNTS_GTF) %.bam
	featureCounts -o $@ -a $< $(FEATURECOUNTS_PARAM) --tmpDir $(TMP_DIR) -T $(NCPUS) $^2

%.read_distribution.txt: %.bam $(REFERENCE_RSEQC_HOUSEKEEPING_BED)
	read_distribution.py  -i $< -r $^2 > $@

%.junctionSaturation_plot.r: %.bam $(REFERENCE_RSEQC_HOUSEKEEPING_BED)
	junction_saturation.py -i $< -r $^2 -o rseqc/$*

%.read_duplication.xls: %.bam %.bam.bai
	read_duplication.py -i $< -o $*

%.geneBodyCoverage.txt: %.bam $(GENCODE_DIR)/rseqc.HouseKeepingGenes.bed.gz %.bam.bai
	docker run -u `id -u`:`id -g` --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) quay.io/biocontainers/rseqc:4.0.0--py38h0213d0e_0  bash -c "cd $(PWD); geneBody_coverage.py -i $< -r <(zcat $^2) -o rseqc/$*"
	sed -i 's|Aligned.sortedByCoord.out|$*|' $@

rseqc/%.inner_distance.txt: STAR/%.STAR/Aligned.sortedByCoord.out.bam
	mkdir -p $$(dirname $@)
	inner_distance.py -i $< -o rseqc/$* -r $(RSEQC_REF_BED

%.ribo.ex.bam: %.bam $(RIBOSOMAL_BED_FILE) %.bam.bai
	docker run -u `id -u`:`id -g` --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) quay.io/biocontainers/rseqc:4.0.0--py38h0213d0e_0  bash -c "cd $(PWD); split_bam.py -i $<  -r $^2 -o $*.ribo > $@.summary"
%.ribo.in.bam: %.ribo.ex.bam
	@echo done
%.ribo.ex.bam.summary: %.ribo.ex.bam
	@echo done

%.bam.read_count: %.bam
	samtools view $< | cut -f -1 | count > $@ 
%.uniq_map.bam: %.bam
	samtools view -h -F 260 $< | bawk '$$1~/^@/ || $$12=="NH:i:1"' | samtools view -O bam > $@   
  #-F 260  only output reads that are not unmapped (flag 4 is not set) and only primary alignment (flag 256 is not set) http://www.metagenomics.wiki/tools/samtools/number-of-reads-in-bam-file*
%.uniq_map.bam.read_count: %.bam
	samtools view    -F 260 $< | bawk '$$1~/^@/ || $$12=="NH:i:1"' \      
	| cut -f -1 | count > $@ 
	#-F 260  only output reads that are not unmapped (flag 4 is not set) and only primary alignment (flag 256 is not set) http://www.metagenomics.wiki/tools/samtools/number-of-reads-in-bam-file*
%.multi_map.bam.read_count: %.bam
	samtools view    -F 4  $< | bawk '$$1~/^@/ || $$12!="NH:i:1"' \
	| cut -f -1 | count > $@
	#*-F 4 only output reads that not unmapped (flag 4 is not set)*
%.unmap.bam: %.bam
	samtools view -f 4 $< -b > $@
  
######################
#
#   NCBI BLAST
#

%.megablast-NT.gz: %.fa.gz
	export BLASTDB=$(BLASTDB);\
	blastn -task megablast -word_size 16 -evalue 0.01 -query <(zcat $<) -num_threads $(NCPUS) -db nt -lcase_masking \
		-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms' \
	| gzip > $@
%.dc_megablast-NT.gz: %.fa.gz
	export BLASTDB=$(BLASTDB);\
	blastn -task dc_megablast -evalue 0.01 -query <(zcat $<) -num_threads $(NCPUS) -db nt \
		-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms' \
	| gzip > $@
%.blastn-NT.gz: %.fa.gz
	export BLASTDB=$(BLASTDB);\
	blastn -task blastn -evalue 0.01 -query <(zcat $<) -num_threads $(NCPUS) -db nt  \
		-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms' \
	| gzip > $@

######################
#
#   File indexing
#
  
  
%.bam.bai: %.bam
	samtools index $<
