# miRNA Analysis Pipeline

## Introduction
The GDC miRNA quantification analysis makes use of a modified version of the profiling pipeline that the British Columbia Genome Sciences Centre developed.  The pipeline generates TCGA-formatted miRNAseq data.  The first step is read alignment. The tool then compares the individual reads to sequence feature annotations in  miRBase and UCSC. Of note, however, the tool only annotates those reads that have an exact match with known miRNAs in miRBase and should therefore not be considered for novel miRNA identification or mismatched alignments.

For more information see [BCGSC's GitHub](https://github.com/bcgsc/mirna) or the [original publication](http://nar.oxfordjournals.org/content/early/2015/08/13/nar.gkv808.full).


## GDC Analysis Steps

### Alignment Workflow
The miRNA pipeline begins with the <a href="/Data_Dictionary/viewer/#?view=table-definition-view&id=alignment_workflow">Alignment Workflow</a>, which in the case of miRNA uses BWA-aln.  This outputs one BAM file for each read group in the input.

Input - [Submitted Unaligned Reads](/Data_Dictionary/viewer/#?view=table-definition-view&id=submitted_unaligned_reads) or [Submitted Aligned Reads](/Data_Dictionary/viewer/#?view=table-definition-view&id=submitted_aligned_reads) (format FASTQ or BAM)
</br>
Output - [Aligned Reads](/Data_Dictionary/viewer/#?view=table-definition-view&id=aligned_reads) (format BAM)

### miRNA Expression Workflow
Following alignment, BAM files are processed through the [miRNA Expression Workflow](/Data_Dictionary/viewer/#?view=table-definition-view&id=mirna_expression_workflow).

The outputs of the miRNA profiling pipeline report raw read counts and counts normalized to reads per million mapped reads (RPM) in two separate files mirnas.quantification.txt and isoforms.quantification.txt. The former contains summed expression for all reads aligned to known miRNAs in the miRBase reference. If there are multiple alignments to different miRNAs or different regions of the same miRNA, the read is flagged as cross-mapped and every miRNA annotation is preserved. The latter contains observed isoforms.

Input - [Aligned Reads](/Data_Dictionary/viewer/#?view=table-definition-view&id=aligned_reads) (format BAM)
</br>
Output - [miRNA Expression](/Data_Dictionary/viewer/#?view=table-definition-view&id=mirna_expression) (format TXT)
