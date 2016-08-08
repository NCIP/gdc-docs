# mRNA Analysis Pipeline

## Introduction
The GDC mRNA quantification analysis pipeline measures gene and exon level expression in [HT-Seq](http://www-huber.embl.de/HTSeq/doc/overview.html) raw mapping count, Fragments per Kilobase of transcript per Million mapped reads (FPKM) and FPKM-UQ (upper quartile normalization).  These values are generated through this pipeline by first aligning reads to the GRCh38 reference genome and then by quantifying the mapped reads.

For more information see (is there some documentation we should link to?)

## Data Processing Steps

### Alignment Workflow
The mRNA Analysis pipeline begins with the [Alignment Workflow](/Data_Dictionary/viewer/#?view=table-definition-view&id=alignment_workflow), which is performed using a 2-pass method with [STAR](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf). STAR aligns each [read group](/Data_Dictionary/viewer/#?view=table-definition-view&id=read_group) separately and then merges the resulting alignments into one. Following the methods used by the International Cancer Genome Consortium [ICGC](https://icgc.org/) ([github](https://github.com/akahles/icgc_rnaseq_align)), the 2-pass method includes a splice junction detection step, which is used to generate the final alignment. This workflow outputs a BAM file, which contains both aligned and unaligned reads. Quality assessment is performed pre-alignment with [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and post-alignment with [RNA-SeQC](https://www.broadinstitute.org/cancer/cga/rna-seqc) and [Picard Tools](http://broadinstitute.github.io/picard/).

![mRNA Align](NEW_SMALLER_PHOTO_URL)

| I/O | Entity | Format |
|---|---|---|
| Input | [Submitted Unaligned Reads](/Data_Dictionary/viewer/#?view=table-definition-view&id=submitted_unaligned_reads) or [Submitted Aligned Reads](/Data_Dictionary/viewer/#?view=table-definition-view&id=submitted_aligned_reads) |  FASTQ or BAM |
| Output | [Aligned Reads](/Data_Dictionary/viewer/#?view=table-definition-view&id=aligned_reads) | BAM  |


### mRNA Expression Workflow
Following alignment, BAM files are processed through the [RNA Expression Workflow](/Data_Dictionary/viewer/#?view=table-definition-view&id=rna_expression_workflow).

First the BAM files are filtered for reads that were aligned to protein-coding genes using the [samtools](http://samtools.sourceforge.net) view function. The reads mapped to each protein-coding gene are enumerated using HT-Seq count. The number of reads mapped to each gene (raw count), gene length and the total number of reads mapped to protein-coding genes are used to calculate normalized expression levels using FPKM and FPKM-UQ (upper quartile normalization).

![mRNA Quant](NEW_SMALLER_PHOTO_URL)

| I/O | Entity | Format |
|---|---|---|
| Input | [Aligned Reads](/Data_Dictionary/viewer/#?view=table-definition-view&id=aligned_reads) |  BAM |
| Output | [Gene Expression (count/FPKM/FPKM-UQ)](/Data_Dictionary/viewer/#?view=table-definition-view&id=gene_expression) | TXT  |
