# GDC Visualization Suite, User Acceptance Testing, April 2017

__URL: https://gdc-portal-staging.datacommons.io__

This guide describes the new data visualization features that can be accessed at the GDC Data Portal. Visualization features include detailed descriptions of mutations, genes, and their frequency; graphical representation of mutation positions, and dynamic survival analysis plots. Additionally, new API endpoints are available for users to programmatically retrieve the data used to generate the visualization features.

The mutation-based visualization features are derived from open-access MAF files that were produced by GDC variant-calling pipelines. The format of these MAF files was developed by and for the GDC and is outlined in the MAF Format documentation.

---

## GDC Portal Home Page

The GDC Portal home page is the entry way to accessing data in the Genomic Data Commons.

[![GDC Home Page](images/GDC-Home-Page.png)](images/GDC-Home-Page.png "Click to see the full image.")

The Genomic Data Commons can be queried by project or by primary site (i.e Brain). Data can be narrowed down in a few ways listed below.

* __Projects__: The projects link directs users to the Project List Page, which gives an overall summary of project-level information.

* __Repository__: The repository link allows users to see the data files available at the GDC and apply file/case filters to narrow down their search.

* __Human Outline__: The home page displays a human anatomical outline that can be used to refine their search. Choosing an associated organ will direct the user to a listing of all projects associated with that primary site. For example, clicking on the human brain will show only cases and projects associated with brain cancer (TCGA-GBM and TCGA-LGG).  The number of cases associated with each primary site is also displayed here and separated by project.

---

# Project List Page

The project list page displays statistics about the projects that are available at the GDC.

[![Project List Graphs](images/GDC-PL-Graphs.png)](images/GDC-PL-Graphs.png "Click to see the full image.")

## Top Mutated Genes in Selected Projects

This dynamically generated bar graph shows the ten genes with the most mutations across all projects. The bars represent the frequency of each mutation and is broken down into different colored segments by project and disease type. The graphic is updated as filters are applied for projects, programs, disease types, and data categories available in the project.

Hovering the cursor over each bar will display information about the number of cases affected by the disease type and clicking on each bar will bring the user to the Gene Summary Page page for the gene associated with the mutation.

## Case Distribution per Project

A pie graph displays the relative number of cases for each primary site (inner circle), which is further divided by project (outer circle). Hovering the cursor over each portion of the graph will display the primary site or project with the number of associated cases. Filtering projects at the left panel will update the pie chart.

## Project Table

[![Project Table](images/GDC-PL-Table.png)](images/GDC-PL-Table.png "Click to see the full image.")

General information about each project is displayed below the graphs. The numbers of available cases per data category for each project are displayed as links. Clicking on a specific project will bring a user to that project's Project Detail Page, which provides an overview of all cases, files and annotations available for the project.

---
# Project Detail Page

Each project has a detail page that provides an overview of all available cases, files, and annotations available. Clicking on the numbers in the summary table will display the corresponding data.

[![Project Detail Page](images/gdc-project-entity-page.png)](images/gdc-project-entity-page.png "Click to see the full image.")

Three download buttons in the top right corner of the screen allow the user to download the entire project dataset, along with the associated project metadata:

* __Download Manifest__: Downloads a manifest for all data files available in the project. The manifest can be used with the GDC Data Transfer Tool to download the files.
* __Download Clinical__: Downloads clinical metadata about all cases in the project.
* __Download Biospecimen__: Downloads biospecimen metadata associated with all cases in the project.

## Most Frequently Mutated Genes

The project detail page also reports the genes that have somatic mutations in the greatest numbers of cases in a graphical and tabular format.

[![Frequently Mutated Genes](images/gdc-frequently-mutated-genes.png)](images/gdc-frequently-mutated-genes.png "Click to see the full image.")

The top of this section contains a survival plot of all the cases within the specified project, in additiona to a bar graph of the most frequently mutated genes. Hovering over each bar in the plot will display information about the number of cases affected. Users may choose to download the underlying data in JSON or TSV format or an image of the graph in SVG or PNG format by clicking the `download` icon at the top of each graph.

Below these graphs is a tabular view of the genes affected, which includes the following information:

* __Symbol:__ The gene symbol, which links to the Gene Summary Page
* __Name:__ Full name of the gene
* __Cytoband:__ The location of the mutation on the chromosome in terms of Giemsa-stained samples.
* __# Affected Cases in Project:__ The number of cases within the project that contain a mutation on this gene
* __# Affected Cases Across all Projects:__ The number of cases within all the projects in the GDC that contain a mutation on this gene.  Clicking the red arrow will display the cases broken down by project
* __# Mutations:__ The number of SSMs (simple somatic mutations) detected in that gene
* __Survival Analysis:__ An icon that, when clicked, will plot the survival rate between cases in the project with mutated and non-mutated forms of the gene

### Survival Analysis

Survival analysis is used to analyze the occurrence of event data over time.  In the GDC, survival analysis is performed on the mortality of the cases. Survival analysis requires:

*  Data on the time to a particular event (days to death or last follow up)
    * Fields:  __diagnoses.days_to_death__ and __diagnoses.days_to_last_follow_up__
*  Information on whether the event has occurred (alive/deceased)
    * Fields:  __diagnoses.vital_status__
*  Data split into different categories or groups (i.e. gender, etc.)
    * Fields:  __demographic.gender__

The survival analysis in the GDC uses a Kaplan-Meier estimator:

[![Kaplan-Meier Estimator](images/gdc-kaplan-meier-estimator.png)](images/gdc-kaplan-meier-estimator "Click to see the full image.")

Where:

 * S(t<sub>i</sub>) is the estimated survival probability for any particular one of the t time periods
 * n<sub>i</sub> is the number of subjects at risk at the beginning of time period t<sub>i</sub>
 * and d<sub>i</sub> is the number of subjects who die during time period t<sub>i</sub>

The table below is an example data set to calculate survival for a set of seven cases:

[![Sample Survival Analysis Table](images/gdc-sample-survival-table.png)](images/gdc-sample-survival-table.png "Click to see the full image.")

The calculated cumulated survival probability can be plotted against the interval to obtain a survival plot like the one shown below.

[![Sample Survival Analysis Plot](images/gdc-survival-plot.png)](images/gdc-survival-plot.png "Click to see the full image.")


## OncoGrid

The project detail page includes an OncoGrid plot of the cases with the most mutations, for the top 50 mutated genes affected by high impact mutations. Genes displayed on the left of the grid (Y-axis) correspond to individual cases on the bottom of the grid (X-axis).

[![Oncogrid](images/gdc-oncogrid.png)](images/gdc-oncogrid.png "Click to see the full image.")

The grid is color-coded with a legend at the top left which describes what type of mutation consequence is observed for each gene/case combination. Clinical information and the available data for each case are available at the bottom of the grid.

The right side of the grid displays additional information about the genes:

* __Gene Sets:__ Describes whether a gene is part of the [Gene Census](http://cancer.sanger.ac.uk/census/).  (The cancer Gene Census is an ongoing effort to catalogue those genes for which mutations have been causally implicated in cancer)
* __GDC:__ Heat-map of all cases in the GDC affected with a mutation in this gene

### OncoGrid Options

To facilitate readability and comparisons, drag-and-drop can be used to reorder the gene rows.  Double clicking a row in the "# Cases Affected" bar at the right side of the graphic will bring the user to the respective Gene Summary Page page. Hovering over a cell will display information about the mutation such as its ID, affected case, and biological consequence. Clicking on the cell will bring the user to the respective Mutation Summary page.

A tool bar at the top right of the graphic allows the user to export the data as a JSON object, PNG image, or SVG image.  Seven buttons are available in this toolbar:

* __Download:__ Users can choose to export the contents either to a static image file (PNG or SVG format) or the underlying data in JSON format
* __Refresh:__ Sets all OncoGrid rows, columns, and zoom levels back to their initial positions
* __Cluster Data:__ Clusters the rows and columns to place mutated genes with the same cases and cases with the same mutated genes together
* __Toggle Heatmap:__ The view can be toggled between cells representing mutation consequences or number of mutations in each gene
* __Toggle Gridlines:__ Turn the gridlines on and off
* __Toggle Crosshairs:__ Turns crosshairs on, so that users can zoom into specific sections of the OncoGrid
* __Fullscreen:__ Turns Fullscreen mode on/off

## Most Frequent Mutations
The project detail page also displays the 20 most frequent mutations in the project as a bar graph that indicates the number of cases that are affected by each mutation. Hovering over each bar in the plot will show information about the number of cases affected.

[![Most Frequent Mutations](images/gdc-frequent-mutations.png)](images/gdc-frequent-mutations.png "Click to see the full image.")

A table is displayed below that lists information about each mutation:

* __ID:__ A UUID for the mutation assigned by the GDC, when clicked will bring a user to the Mutation Summary Page
* __DNA Change:__ The chromosome and starting coordinates of the mutation are displayed along with the nucleotide differences between the reference and tumor allele
* __Type:__ A general classification of the mutation
* __Consequences:__ The effects the mutation has on the gene coding for a protein (i.e. synonymous, missense, non-coding transcript).  A link to the Gene Summary Page for the gene affected by the mutation is included
* __# Affected Cases in Project:__ The number of affected cases in the project
* __# Affected Cases in Across all Projects:__ The number of affected cases, expressed as number across all projects. Choosing the arrow next to the percentage will display a breakdown of each affected project
* __Impact:__ A subjective classification of the severity of the variant consequence. The categories are:
    * __HIGH__: The variant is assumed to have high (disruptive) impact in the protein, probably causing protein truncation, loss of function, or triggering nonsense mediated decay
    * __MODERATE__: A non-disruptive variant that might change protein effectiveness
    * __LOW__: Assumed to be mostly harmless or unlikely to change protein behavior
* __Survival Analysis:__ An icon that when clicked, will plot the survival rate between the gene's mutated and non-mutated cases

## Most Affected Cases

The final section of the project detail page is a display of the top 20 cases in a specified project, with the greatest number of affected genes.

[![Most Affected Cases](images/gdc-most-affected-cases.png)](images/gdc-most-affected-cases.png "Click to see the full image.")

Below the bar graph is a table contains information about these cases:

* __UUID:__ The UUID of the case, which links to the Case Summary Page
* __Site:__  The anatomical location of the site affected
* __Gender:__ Text designations that identify gender. Gender is described as the assemblage of properties that distinguish people on the basis of their societal roles
* __Age at Diagnosis:__ Age at the time of diagnosis expressed in number of days since birth
* __Stage:__ The extent of a cancer in the body. Staging is usually based on the size of the tumor, whether lymph nodes contain cancer, and whether the cancer has spread from the original site to other parts of the body. The accepted values for tumor_stage depend on the tumor site, type, and accepted staging system
* __Survival:__ The number of days until death
* __Last Follow Up:__ Time interval from the date of last follow up to the date of initial pathologic diagnosis, represented as a calculated number of days
* __Available Files per Data Category:__ Five columns displaying the number of files available in each of the five data categories.  These link to the files for the specific case.
* __Genes:__ The number of genes affected by mutations for the case

---
# Case Summary Page

The Case Summary Page displays case details including the project and disease information, data files that are available for that case, and the experimental strategies employed. A button in the top-right corner of the page allows the user to add all files associated with the case to the file cart.

[![Case Page](images/gdc-case-entity-page.png)](images/gdc-case-entity-page.png "Click to see the full image.")

## Clinical and Biospecimen Information

The page also provides clinical and biospecimen information about that case. Links to export clinical and biospecimen information in JSON format are provided.

[![Case Page, Clinical and Biospecimen](images/gdc-case-clinical-biospecimen.png)](images/gdc-case-clinical-biospecimen.png "Click to see the full image.")

For clinical records that support multiple records of the same type (Diagnoses, Family Histories, or Exposures), a UUID of the record is provided on the left hand side of the corresponding tab, allowing the user to select the entry of interest.

### Biospecimen Search

A search filter just below the biospecimen section can be used to find and filter biospecimen data. The wildcard search will highlight entities in the tree that match the characters typed. This will search both the case submitter ID, as well as the additional metadata for each entity. For example, searching 'Primary Tumor' will highlight samples that match that type.

[![Biospecimen Search](images/gdc-case-biospecimen-search.png)](images/gdc-case-biospecimen-search.png "Click to see the full image.")

---
# Gene Summary Page

The Gene Summary Page describes each gene with mutation data featured at the GDC and provides results related to the analyses that are performed on these genes.

## Summary

The summary section of the gene page contains the following information:

[![Gene Summary](images/GDC-Gene-Summary.png)](images/GDC-Gene-Summary.png "Click to see the full image.")

* __Symbol:__ The gene symbol
* __Name:__ Full name of the gene
* __Synonyms:__ Synonyms of the gene name or symbol, if available
* __Type:__ A broad classification of the gene
* __Location:__ The chromosome on which the gene is located and its coordinates
* __Strand:__ If the gene is located on the forward (+) or reverse (-) strand
* __Description:__ A description of gene function and downstream consequences of gene alteration

## External References

A list with links that lead to external databases with additional information about each gene is displayed here. These external databases include: [Entrez](https://www.ncbi.nlm.nih.gov/gquery/), [Hugo Gene Nomenclature Committee](http://www.genenames.org/), [Online Mendelian Inheritance in Man](https://www.omim.org/), and [Uniprotkb SwissProt](http://www.uniprot.org/).

## Cancer Distribution

[![Cancer Distribution](images/GDC-Gene-CancerDist.png)](images/GDC-Gene-CancerDist.png "Click to see the full image.")

A table and bar graph displayed that shows how many cases are affected by mutations within the gene as a ratio and percentage. Each row/bar represents the number of cases for each project.

## Protein Viewer

[![Protein Plot](images/GDC-Gene-ProteinGraph.png)](images/GDC-Gene-ProteinGraph.png "Click to see the full image.")

Mutations and their frequency across cases are mapped to a graphical visualization of protein-coding regions with a lollipop plot. Pfam domains are highlighted along the x-axis to assign functionality to specific protein-coding regions. The bottom track represents a view of the full gene length. Different transcripts can be selected by using the drop-down menu above the plot.

The panel to the right of the plot allows the plot to be filtered by mutation consequences or impact.  The plot will dynamically change as filters are applied.  Mutation consequence and impact is denoted in the plot by color.

The plot can be viewed at different zoom levels by clicking and dragging across the x-axis, clicking and dragging across the bottom track, or double clicking the pfam domain IDs. The `Reset` button can be used to bring the zoom level back to its original position. The plot can also be exported as a PNG image, SVG image or as JSON formatted text by choosing the `Download` button above the plot.

## Most Frequent Mutations

The 20 most frequent mutations in the gene are displayed as a bar graph that indicates the number of cases that share each mutation.

[![Gene MFM](images/GDC-Gene-MFM.png)](images/GDC-Gene-MFM.png "Click to see the full image.")

A table is displayed below that lists information about each mutation including:

* __ID:__ A UUID Code for the mutation assigned by the GDC, when clicked will bring a user to the Mutation Summary Page
* __DNA Change:__ The chromosome and starting coordinates of the mutation are displayed along with the nucleotide differences between the reference and tumor allele
* __Type:__ A general classification of the mutation
* __Consequences:__ The effects the mutation has on the gene coding for a protein (i.e. synonymous, missense, non-coding transcript)
* __# Affected Cases:__ The number of affected cases, expressed as number across all projects. Choosing the arrow next to the percentage will expand the selection with a breakdown of each affected project
* __Impact:__ A subjective classification of the severity of the variant consequence. The categories are:
  - __HIGH__: The variant is assumed to have high (disruptive) impact in the protein, probably causing protein truncation, loss of function or triggering nonsense mediated decay
  - __MODERATE__: A non-disruptive variant that might change protein effectiveness
  - __LOW__: Assumed to be mostly harmless or unlikely to change protein behavior

---
# Mutation Summary Page

The Mutation Summary Page contains information about one somatic mutation and how it affects the associated gene. Each mutation is identified by its chromosomal position and nucleotide-level change.

## Summary

[![Mutation Summary](images/GDC-Mutation-Summary.png)](images/GDC-Mutation-Summary.png "Click to see the full image.")

- __ID:__ A unique identifier (UUID) for this mutation
- __DNA Change:__ Denotes the chromosome number, position, and nucleotide change of the mutation
- __Type:__ A broad categorization of the mutation
- __Reference Genome Assembly:__ The reference genome in which the chromosomal position refers to
- __Allele in the Reference Assembly:__ The nucleotide(s) that compose the site in the reference assembly
- __Functional Impact:__ A subjective classification of the severity of the variant consequence. The categories are:
  - __HIGH__: The variant is assumed to have high (disruptive) impact in the protein, probably causing protein truncation, loss of function or triggering nonsense mediated decay
  - __MODERATE__: A non-disruptive variant that might change protein effectiveness
  - __LOW__: Assumed to be mostly harmless or unlikely to change protein behavior

### External References

A separate panel contains links to information about the gene that the mutation effects in external databases. These databases include: [Entrez](https://www.ncbi.nlm.nih.gov/gquery/), [Hugo Gene Nomenclature Committee](http://www.genenames.org/), [Online Mendelian Inheritance in Man](https://www.omim.org/), and [Uniprotkb SwissProt](http://www.uniprot.org/).

## Consequences

[![Mutation Consequences](images/GDC-Mutation-Consequences.png)](images/GDC-Mutation-Consequences.png "Click to see the full image.")

The consequences of the mutation are displayed in a table. The fields that detail each mutation are listed below:

* __Gene:__ The symbol for the affected gene
* __AA Change:__ Details on the amino acid change, including compounds and position, if applicable
* __Consequence:__ The biological consequence of each mutation
* __Coding DNA Change:__ The specific nucleotide change and position of the mutation within the gene
* __Strand:__ If the gene is located on the forward (+) or reverse (-) strand
* __Transcript(s):__ The transcript(s) affected by the mutation. Each contains a link to the [Ensembl](https://www.ensembl.org) entry for the transcript

---
## Cancer Distribution

[![Mutation Distribution](images/GDC-Mutation-CancerDist.png)](images/GDC-Mutation-CancerDist.png "Click to see the full image.")

A bar graph detailing the number of cases affected by the mutation across all projects is generated dynamically.  Below the graph is a table with information about how the mutation affects each project, which can be exported as a JSON object. The table contains the following fields:

* __Project ID__: The ID for a specific project
* __Disease Type__: The disease associated with the project
* __Site__: The anatomical site affected by the disease
* __# Affected Cases__: The number of affected cases and total number of cases displayed as a fraction and percentage

## Protein Viewer

[![Mutation Protein Graph](images/GDC-Mutation-ProteinGraph.png)](images/GDC-Mutation-ProteinGraph.png "Click to see the full image.")

The protein viewer displays a plot representing the position of mutations along the polypeptide chain associated with the mutation. The y-axis represents the number of cases that exhibit each mutation, whereas the x-axis represents the polypeptide chain sequence. [Pfam domains](http://pfam.xfam.org/) that were identified along the polypeptide chain are identified with colored rectangles labeled with pfam IDs. See the Gene Summary Page for additional details about the protein viewer.

---
# Visualization API Endpoint Additions

The GDC Visualization Suite uses the same API as the rest of the Data Portal and takes advantage of several new endpoints:

* __ssms:__ The simple somatic mutation (`ssms`) endpoint allows users to access information about each somatic point mutation. For example, a `ssm` would represent the transition of C to T at position 52000 of chromosome 1.
* __ssm_occurrences:__ A SSM entity as applied to a single instance (case). An example of a `ssm occurrence` would be that the transition of C to T at position 52000 of chromosome 1 occurred in patient TCGA-XX-XXXX.
* __genes:__ The `genes` endpoint allows users to access in-depth information about each gene.

The methods for retrieving information from these endpoint are very similar to those used for the `cases` and `files` endpoints. These methods are explored in depth in the [API Search and Retrieval](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/) documentation. The `_mapping` parameter can also be used with each of these endpoints to generate a list of potential fields.  For example:

`https://api.gdc.cancer.gov/ssms/_mapping`

While it is not an endpoint, the `observation` entity is featured in the visualization section of the API. The `observation` entity provides information from the MAF file, such as read depth and normal genotype, that supports the validity of the associated `ssm`.

## Analysis Endpoints

In addition the `ssms`, `ssm_occurrences`, and `genes` endpoints mentioned previously, several `analysis` endpoints were designed to quickly retrieve specific datasets used for visualization display.  

* __analysis/survival__
* __analysis/top_cases_counts_by_genes__
* __analysis/top_mutated_genes_by_project__
* __analysis/top_mutated_cases_by_gene__
* __analysis/top_mutated_cases_by_ssm__
* __analysis/mutated_cases_count_by_project__

For a set of API query examples, go to:  

https://github.com/NCIP/gdc-docs/blob/Visualization-Docs/docs/Visualization_Documentation/API_Endpoints.md

---
# GDC MAF Format v.1.0.0

## Introduction

Mutation Annotation Format (MAF) is a tab-delimited text file with aggregated mutation information from [VCF Files](VCF_Format.md) and are generated on a project-level. MAF files are produced through the [Somatic Aggregation Workflow](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=somatic_aggregation_workflow&_top=1)  The GDC produces MAF files at two permission levels: __protected__ and __somatic__ (or open-access). One MAF files is produced per variant calling pipeline per GDC project. MAFs are produced by aggregating the GDC annotated VCF files generated from one pipeline for one project.

Annotated VCF files often have variants reported on multiple transcripts whereas the MAF files generated from the VCFs (\*protected.maf) only report the most critically affected one. Somatic MAFs (\*somatic.maf), which are also known as [Masked Somatic Mutation](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=masked_somatic_mutation) files, are further processed to remove lower quality and potential germline variants. For tumor samples that contain variants from multiple combinations of tumor-normal aliquot pairs, only one pair is selected in the Somatic MAF based on their sample type. Somatic MAFs are publicly available and can be freely distributed within the boundaries of the [GDC Data Access Policies](https://gdc.cancer.gov/access-data/data-access-policies).

The GDC MAF file format is based on the [TCGA Mutation Annotation Format](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification) specifications, with additional columns included.

__Note:__ The criteria for allowing mutations into open-access are purposefully implemented to overcompensate and filter out germline variants. If omission of true-positive somatic mutations is a concern, the GDC recommends using protected MAFs.  

## Somatic MAF File Generation

The process for modifying a protected MAF into a somatic MAF is as follows:

*  Aliquot Selection: only one tumor-normal pair are selected for each tumor sample based on the plate number, sample type, analyte type and other features extracted from tumor TCGA aliquot barcode.
*  Low quality variant filtering and germline masking:
    1. Variants with __Mutation_Status != 'Somatic'__ or __GDC_FILTER = 'Gapfiller', 'ContEst', 'multiallelic', 'nonselectedaliquot', 'BCR_Duplicate' or 'BadSeq'__ are __removed__.
    2. Remaining variants with __GDC_Valid_Somatic = True__ are __included__ in the Somatic MAF.
    3. Remaining variants with __FILTER != 'panel_of_normals' or PASS__ are __removed__. Note that the `FILTER != panel_of_normals` value is only relevant for the variants generated from the MuTect2 pipeline.
    4. Remaining variants with __MC3_Overlap = True__ are __included__ in the Somatic MAF.
    5. Remaining variants with __GDC_FILTER = 'ndp', 'NonExonic', 'bitgt', or 'gdc_pon'__ are __removed__.
    6. Remaining variants with __SOMATIC != null__ are __included__ in the Somatic MAF.
    7. Remaining variants with __dbSNP_RS = 'novel' or null__ are __included__ in the Somatic MAF.
    8. Remaining variants are __removed__.
* Removal of the following columns:
    * vcf_info
    * vcf_format
    * vcf_tumor_gt
    * vcf_normal_gt
    * GDC_Valid_Somatic
* Set values to be blank in the following columns that may contain information about germline genotypes:
    * Match_Norm_Seq_Allele1
    * Match_Norm_Seq_Allele2
    * Match_Norm_Validation_Allele1
    * Match_Norm_Validation_Allele2
    * n_ref_count
    * n_alt_count

---

![Somatic MAF Generation](images/ProtectedMAF4.png)

---
---

## Protected MAF File Structure

The table below describes the columns in a protected MAF and their definitions. Note that the somatic (open-access) MAF structure is the same except for having the last six columns removed.

| Column | Description |
|---|---|
| 1 - Hugo_Symbol | [HUGO](http://www.genenames.org/) symbol for the gene (HUGO symbols are always in all caps). "Unknown" is used for regions that do not correspond to a gene  |
| 2 - Entrez_Gene_Id | [Entrez gene](https://www.ncbi.nlm.nih.gov/gene) ID (an integer).  "0" is used for regions that do not correspond to a gene region or Ensembl ID |
| 3 - Center| One or more genome sequencing center reporting the variant |
| 4 - NCBI_Build| The reference genome used for the alignment (GRCh38) |
| 5 - Chromosome | The affected chromosome (chr1) |
| 6 - Start_Position | Lowest numeric position of the reported variant on the genomic reference sequence. Mutation start coordinate |
| 7 - End_Position | Highest numeric genomic position of the reported variant on the genomic reference sequence. Mutation end coordinate |
| 8 - Strand | Genomic strand of the reported allele. Currently, all variants will report the positive strand: '+' |
| 9 - Variant_Classification | Translational effect of variant allele  |
| 10 - Variant_Type | Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP (di-nucleotide polymorphism) but for three consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of four or more (SNP, DNP, TNP, ONP, INS, DEL, or Consolidated)|
| 11 - Reference_Allele | The plus strand reference allele at this position. Includes the deleted sequence for a deletion or "-" for an insertion |
| 12 - Tumor_Seq_Allele1 | Primary data genotype for tumor sequencing (discovery) allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases  |
| 13 - Tumor_Seq_Allele2 | Tumor sequencing (discovery) allele 2  |
| 14 - dbSNP_RS| The rs-IDs from the [dbSNP](https://www.ncbi.nlm.nih.gov/projects/SNP/) database, "novel" if not found in any database used, or null if there is no dbSNP record, but it is found in other databases|
| 15 - dbSNP_Val_Status | The dbSNP validation status is reported as a semicolon-separated list of statuses. The union of all rs-IDs is taken when there are multiple |
| 16 - Tumor_Sample_Barcode| Aliquot barcode for the tumor sample |
| 17 - Matched_Norm_Sample_Barcode| Aliquot barcode for the matched normal sample |
| 18 - Match_Norm_Seq_Allele1 | Primary data genotype. Matched normal sequencing allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases (cleared in somatic MAF) |
| 19 - Match_Norm_Seq_Allele2| Matched normal sequencing allele 2 |
| 20 - Tumor_Validation_Allele1| Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases  |
| 21 - Tumor_Validation_Allele2 | Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 2 |
| 22 - Match_Norm_Validation_Allele1 | Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases (cleared in somatic MAF) |
| 23 - Match_Norm_Validation_Allele2 | Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 2 (cleared in somatic MAF) |
| 24 - Verification_Status| Second pass results from independent attempt using same methods as primary data source. Generally reserved for 3730 Sanger Sequencing |
| 25 - Validation_Status | Second pass results from orthogonal technology |
| 26 - Mutation_Status | An assessment of the mutation as somatic, germline, LOH, post transcriptional modification, unknown, or none. The values allowed in this field are constrained by the value in the Validation_Status field |
| 27 - Sequencing_Phase | TCGA sequencing phase (if applicable). Phase should change under any circumstance that the targets under consideration change |
| 28 - Sequence_Source | Molecular assay type used to produce the analytes used for sequencing. Allowed values are a subset of the SRA 1.5 library_strategy field values. This subset matches those used at CGHub |
| 29 - Validation_Method | The assay platforms used for the validation call |
| 30 - Score | Not in use |
| 31 - BAM_File | Not in use |
| 32 - Sequencer | Instrument used to produce primary sequence data |
| 33 - Tumor_Sample_UUID | GDC aliquot UUID for tumor sample |
| 34 - Matched_Norm_Sample_UUID | GDC aliquot UUID for matched normal sample |
| 35 - HGVSc | The coding sequence of the variant in HGVS recommended format |
| 36 - HGVSp | The protein sequence of the variant in HGVS recommended format. "p.=" signifies no change in the protein |
| 37 - HGVSp_Short | Same as the HGVSp column, but using 1-letter amino-acid codes |
| 38 - Transcript_ID| [Ensembl](http://useast.ensembl.org/index.html) ID of the transcript affected by the variant |
| 39 - Exon_Number | The exon number (out of total number) |
| 40 - t_depth | Read depth across this locus in tumor BAM |
| 41 - t_ref_count | Read depth supporting the reference allele in tumor BAM |
| 42 - t_alt_count | Read depth supporting the variant allele in tumor BAM |
| 43 - n_depth | Read depth across this locus in normal BAM |
| 44 - n_ref_count | Read depth supporting the reference allele in normal BAM (cleared in somatic MAF) |
| 45 - n_alt_count | Read depth supporting the variant allele in normal BAM (cleared in somatic MAF) |
| 46 - all_effects | A semicolon delimited list of all possible variant effects, sorted by priority ([SYMBOL,Consequence,HGVSp_Short,Transcript_ID,RefSeq]) |
|47 - Allele| The variant allele used to calculate the consequence |
|48 - Gene|Stable Ensembl ID of affected gene |
|49 - Feature|Stable Ensembl ID of feature (transcript, regulatory, motif) |
|50 - Feature_type| Type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature (or blank) |
|51 - Consequence | Consequence type of this variation; [sequence ontology](http://www.sequenceontology.org/) terms |
|52 - cDNA_position | Relative position of base pair in the cDNA sequence as a fraction. A "-" symbol is displayed as the numerator if the variant does not appear in cDNA |
|53 - CDS_position | Relative position of base pair in coding sequence. A "-" symbol is displayed as the numerator if the variant does not appear in coding sequence |
|54 - Protein_position | Relative position of affected amino acid in protein. A "-" symbol is displayed as the numerator if the variant does not appear in coding sequence |
|55 - Amino_acids | Only given if the variation affects the protein-coding sequence |
|56 - Codons | The alternative codons with the variant base in upper case |
|57 - Existing_variation | Known identifier of existing variation |
|58 - ALLELE_NUM | Allele number from input; 0 is reference, 1 is first alternate etc. |
|59 - DISTANCE | Shortest distance from the variant to transcript |
|60 - TRANSCRIPT_STRAND | The DNA strand (1 or -1) on which the transcript/feature lies |
|61 - SYMBOL | The gene symbol |
|62 - SYMBOL_SOURCE | The source of the gene symbol |
|63 - HGNC_ID | Gene identifier from the HUGO Gene Nomenclature Committee if applicable |
|64 - BIOTYPE | Biotype of transcript |
|65 - CANONICAL | A flag (YES) indicating that the VEP-based canonical transcript, the longest translation, was used for this gene. If not, the value is null |
|66 - CCDS | The [CCDS](https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi) identifier for this transcript, where applicable |
|67 - ENSP | The Ensembl protein identifier of the affected transcript |
|68 - SWISSPROT | [UniProtKB/Swiss-Prot](http://www.uniprot.org/) accession |
|69 - TREMBL | UniProtKB/TrEMBL identifier of protein product |
|70 - UNIPARC | UniParc identifier of protein product |
|71 - RefSeq | RefSeq identifier for this transcript |
|72 - SIFT | The [SIFT](http://sift.jcvi.org/) prediction and/or score, with both given as prediction (score) |
|73 - PolyPhen | The [PolyPhen](http://genetics.bwh.harvard.edu/pph2/) prediction and/or score |
|74 - EXON | The exon number (out of total number) |
|75 - INTRON | The intron number (out of total number) |
|76 - DOMAINS | The source and identifier of any overlapping protein domains |
|77 - GMAF | Non-reference allele and frequency of existing variant in [1000 Genomes](http://www.internationalgenome.org/) |
|78 - AFR_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined African population |
|79 - AMR_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined American population |
|80 - ASN_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined Asian population |
|81 - EAS_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined East Asian population |
|82 - EUR_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined European population |
|83 - SAS_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined South Asian population |
|84 - AA_MAF | Non-reference allele and frequency of existing variant in [NHLBI-ESP](http://evs.gs.washington.edu/EVS/) African American population |
|85 - EA_MAF | Non-reference allele and frequency of existing variant in NHLBI-ESP European American population|
|86 - CLIN_SIG | Clinical significance of variant from dbSNP |
|87 - SOMATIC |Somatic status of each ID reported under Existing_variation (0, 1, or null) |
|88 - PUBMED | Pubmed ID(s) of publications that cite existing variant |
|89 - MOTIF_NAME |The source and identifier of a transcription factor binding profile aligned at this position |
|90 - MOTIF_POS | The relative position of the variation in the aligned TFBP |
|91 - HIGH_INF_POS | A flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP) (Y, N, or null) |
| 92 - MOTIF_SCORE_CHANGE | The difference in motif score of the reference and variant sequences for the TFBP |
| 93 - IMPACT | The impact modifier for the consequence type |
| 94 - PICK | Indicates if this block of consequence data was picked by VEP's [pick feature](http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_pick) (1 or null) |
| 95 - VARIANT_CLASS | Sequence Ontology variant class |
| 96 - TSL | [Transcript support level](http://useast.ensembl.org/Help/Glossary?id=492), which is based on independent RNA analyses |
| 97 - HGVS_OFFSET | Indicates by how many bases the HGVS notations for this variant have been shifted |
| 98 - PHENO|Indicates if existing variant is associated with a phenotype, disease or trait (0, 1, or null) |
| 99 - MINIMISED | Alleles in this variant have been converted to minimal representation before consequence calculation (1 or null) |
| 100 - ExAC_AF | Global Allele Frequency from [ExAC](http://exac.broadinstitute.org/) |
| 101 - ExAC_AF_Adj | Adjusted Global Allele Frequency from ExAC |
| 102 - ExAC_AF_AFR | African/African American Allele Frequency from ExAC |
| 103 - ExAC_AF_AMR | American Allele Frequency from ExAC |
| 104 - ExAC_AF_EAS | East Asian Allele Frequency from ExAC |
| 105 - ExAC_AF_FIN | Finnish Allele Frequency from ExAC |
| 106 - ExAC_AF_NFE | Non-Finnish European Allele Frequency from ExAC |
| 107 - ExAC_AF_OTH | Other Allele Frequency from ExAC |
| 108 - ExAC_AF_SAS | South Asian Allele Frequency from ExAC |
| 109 - GENE_PHENO | Indicates if gene that the variant maps to is associated with a phenotype, disease or trait (0, 1, or null) |
| 110 - FILTER | Copied from input VCF |
| 111 - CONTEXT | The reference allele per VCF specs, and its five flanking base pairs |
| 112 - src_vcf_id | GDC UUID for the input VCF file |
| 113 - tumor_bam_uuid | GDC UUID for the tumor bam file |
| 114 - normal_bam_uuid | GDC UUID for the normal bam file |
| 115 - case_id | GDC UUID for the case |
| 116 - GDC_FILTER | GDC filters applied universally across all MAFs |
| 117 - COSMIC | Overlapping COSMIC variants |
| 118 - MC3_Overlap| Indicates whether this region overlaps with an MC3 variant for the same sample pair |
| 119 - GDC_Validation_Status | GDC implementation of validation checks. See notes section (#5) below for details |
| 120 - GDC_Valid_Somatic | True or False (not in somatic MAF) |
| 121 - vcf_region | Colon separated string containing the CHROM, POS, ID, REF, and ALT columns from the VCF file (e.g., chrZ:20:rs1234:A:T) (not in somatic MAF) |
| 122 - vcf_info | INFO column from VCF (not in somatic MAF) |
| 123 - vcf_format | FORMAT column from VCF (not in somatic MAF) |
| 124 - vcf_tumor_gt | Tumor sample genotype column from VCF (not in somatic MAF) |
| 125 - vcf_normal_gt | Normal sample genotype column from VCF (not in somatic MAF) |  

### Notes About GDC MAF Implementation

1. Column #4 __NCBI_Build__ is GRCh38 by default
2. Column #32 __Sequencer__ includes the sequencers used. If different sequencers were used to generate normal and tumor data, the normal sequencer is listed first.
3. Column #60 VEP name "STRAND" is changed to __TRANSCRIPT_STRAND__ to avoid confusion with Column#8 "Strand"
4. Column #122-125 __vcf_info, vcf_format, vcf_tumor_gt, and vcf_normal_gt__ are the corresponding columns from the VCF files. Including them facilitates parsing specific variant information.
5. Column #119 __GDC_Validation_Status__: GDC also collects TCGA validation sequences. It compares these with variants derived from Next-Generation Sequencing data from the same sample and populates the comparison result in "GDC_Validation_Status".
    * "Valid", if the alternative allele(s) in the tumor validation sequence is(are) the same as GDC variant call
    * "Invalid", if none of the alternative allele(s) in the tumor validation sequence is the same as GDC variant call
    * "Inconclusive" if two alternative allele exists, and one matches while the other does not
    * "Unknown" if no validation sequence exists
6. Column #120 __GDC_Valid_Somatic__ is TRUE if GDC_Validation_Status is "Valid" and the variant is "Somatic" in validation calls.  It is FALSE if these criteria are not met


---
# DAVE Release Notes

## Data Release UAT

* __GDC Product__: GDC DAVE Tools
* __Release Date__: April 18, 2017

### New updates

* None

### Bugs Fixed Since Last Release

* None

### Known Issues and Workarounds

* Features requiring login are disabled
* Some standard GDC Data Portal features are unavailable: download manifest, download clinical/biospecimen metadata, export tables, bam slicing
* Only a single consequence for a mutation is shown on the case, project, and gene pages.  Multiple consequences can be viewed on the mutation entity page
* Some queries may exhibit longer than expected load times
* Mutation data is only available for TCGA projects
* The number of cases displayed/visible on the portal is less than what is in the downloadable MAF files <!--SV-624-->
* Cancer Distribution Plot displays incorrect static title <!--SV-614-->
