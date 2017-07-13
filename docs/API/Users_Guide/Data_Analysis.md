# Data Analysis

The GDC DAVE tools use the same API as the rest of the Data Portal and takes advantage of several new endpoints. Similar to the [GDC Data Portal Exploration](http://docs.gdc.cancer.gov/Data_Portal/Users_Guide/Exploration/) feature, the GDC data analysis endpoints allow API users to programmatically explore data in the GDC using advanced filters at a gene and mutation level. Survival analysis data is also available.  

## Endpoints

The following data analysis endpoints are available from the GDC API:

| __Endpoint__ | __Description__ |
|---|---|
| __/genes__ | Allows users to access summary information about each gene using its Ensembl ID. |
| __/ssms__ | Allows users to access information about each somatic mutation. For example, a `ssm` would represent the transition of C to T at position 52000 of chromosome 1. |
| __/ssm_occurrences__ | A `ssm` entity as applied to a single instance (case). An example of a `ssm occurrence` would be that the transition of C to T at position 52000 of chromosome 1 occurred in patient TCGA-XX-XXXX. |
|__/analysis/top_cases_counts_by_genes__| Returns the number of cases with a mutation in each gene listed in the gene_ids parameter for each project. Note that this endpoint cannot be used with the `format` or `fields` parameters.|
|__/analysis/top_mutated_genes_by_project__| Returns a list of genes that have the most mutations within a given project. |
|__/analysis/top_mutated_cases_by_gene__| Generates information about the cases that are most affected by mutations in a given number of genes |
|__/analysis/mutated_cases_count_by_project__| Returns counts for the number of cases that have associated `ssm` data in each project. The number of affected cases can be found under "case_with_ssm": {"doc_count": $case_count}.|
|__/analysis/survival__| Survival plots can be generated in the Data Portal for different subsets of data, based upon many query factors such as variants, disease type and projects. This endpoint can be used to programmatically retrieve the raw data to generate these plots and apply different filters to the data. (see Survival Example)|

The methods for retrieving information from these endpoints are very similar to those used for the `cases` and `files` endpoints. These methods are explored in depth in the [API Search and Retrieval](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/) documentation. The `_mapping` parameter can also be used with each of these endpoints to generate a list of potential fields.  For example:

`https://api.gdc.cancer.gov/ssms/_mapping`

Note: While it is not an endpoint, the `observation` entity is featured in the visualization section of the API. The `observation` entity provides information from the MAF file, such as read depth and normal genotype, that supports the validity of the associated `ssm`. An example is demonstrated below:

```Shell
curl "https://api.gdc.cancer.gov/ssms/57bb3f2e-ec05-52c2-ab02-7065b7d24849?expand=occurrence.case.observation.read_depth&pretty=true"
```
```Response
{
  "data": {
    "ncbi_build": "GRCh38",
    "occurrence": [
      {
        "case": {
          "observation": [
            {
              "read_depth": {
                "t_ref_count": 321,
                "t_alt_count": 14,
                "t_depth": 335,
                "n_depth": 115
              }
            }
          ]
        }
      }
    ],
    "tumor_allele": "G",
    "mutation_type": "Simple Somatic Mutation",
    "end_position": 14304578,
    "reference_allele": "C",
    "ssm_id": "57bb3f2e-ec05-52c2-ab02-7065b7d24849",
    "start_position": 14304578,
    "mutation_subtype": "Single base substitution",
    "cosmic_id": null,
    "genomic_dna_change": "chr5:g.14304578C>G",
    "gene_aa_change": [
      "TRIO L229V",
      "TRIO L437V",
      "TRIO L447V",
      "TRIO L496V"
    ],
    "chromosome": "chr5"
  },
  "warnings": {}
}
```

## Genes Endpoint Examples

__Example 1:__ A user would like to access information about the gene `ZMPSTE24`, which has an Ensembl gene ID of `ENSG00000084073`.  This would be accomplished by appending `ENSG00000084073` (`gene_id`) to the `genes` endpoint.

```Shell
curl "https://api.gdc.cancer.gov/genes/ENSG00000084073?pretty=true"
```
```Response
{
  "data": {
    "canonical_transcript_length": 3108,
    "description": "This gene encodes a member of the peptidase M48A family. The encoded protein is a zinc metalloproteinase involved in the two step post-translational proteolytic cleavage of carboxy terminal residues of farnesylated prelamin A to form mature lamin A. Mutations in this gene have been associated with mandibuloacral dysplasia and restrictive dermopathy. [provided by RefSeq, Jul 2008]",
    "cytoband": [
      "1p34.2"
    ],
    "gene_start": 40258107,
    "canonical_transcript_length_genomic": 36078,
    "gene_id": "ENSG00000084073",
    "gene_strand": 1,
    "canonical_transcript_length_cds": 1425,
    "gene_chromosome": "1",
    "synonyms": [
      "FACE-1",
      "HGPS",
      "PRO1",
      "STE24",
      "Ste24p"
    ],
    "is_cancer_gene_census": null,
    "biotype": "protein_coding",
    "gene_end": 40294184,
    "canonical_transcript_id": "ENST00000372759",
    "symbol": "ZMPSTE24",
    "name": "zinc metallopeptidase STE24"
  },
  "warnings": {}
}
```

__Example 2:__ A user wants a subset of elements such as a list of coordinates for all genes on chromosome 7. The query can be filtered for only results from chromosome 7 using a JSON-formatted query that is URL-encoded.

```Shell
curl "https://api.gdc.cancer.gov/genes?pretty=true&fields=gene_id,symbol,gene_start,gene_end&format=tsv&size=2000&filters=%7B%0D%0A%22op%22%3A%22in%22%2C%0D%0A%22content%22%3A%7B%0D%0A%22field%22%3A%22gene_chromosome%22%2C%0D%0A%22value%22%3A%5B%0D%0A%227%22%0D%0A%5D%0D%0A%7D%0D%0A%7D"
```
```Response
gene_start      gene_end        symbol  id
28995231        29195451        CPVL    ENSG00000106066
33014114        33062797        NT5C3A  ENSG00000122643
143052320       143053347       OR6V1   ENSG00000225781
100400826       100428992       ZCWPW1  ENSG00000078487
73861159        73865893        WBSCR28 ENSG00000175877
64862999        64864370        EEF1DP4 ENSG00000213640
159231435       159233377       PIP5K1P2        ENSG00000229435
141972631       141973773       TAS2R38 ENSG00000257138
16646131        16706523        BZW2    ENSG00000136261
149239651       149255609       ZNF212  ENSG00000170260
57405025        57405090        MIR3147 ENSG00000266168
130393771       130442433       CEP41   ENSG00000106477
150800403       150805120       TMEM176A        ENSG00000002933
93591573        93911265        GNGT1   ENSG00000127928
117465784       117715971       CFTR    ENSG00000001626
5879827 5886362 OCM     ENSG00000122543
144118461       144119360       OR2A15P ENSG00000239981
30424527        30478784        NOD1    ENSG00000106100
137227341       137343865       PTN     ENSG00000105894
84876554        84876956        HMGN2P11        ENSG00000232605
107470018       107475659       GPR22   ENSG00000172209
31330711        31330896        RP11-463M14.1   ENSG00000271027
78017057        79453574        MAGI2   ENSG00000187391
55736779        55739605        CICP11  ENSG00000237799
142111749       142222324       RP11-1220K2.2   ENSG00000257743
(truncated)
```

## Simple Somatic Mutation Endpoint Examples

__Example 1__: Similar to the `/genes` endpoint, a user would like to retrieve information about the mutation based on its COSMIC ID. This would be accomplished by creating a JSON filter such as:

```Query
 {
   "op":"in",
   "content":{
      "field":"cosmic_id",
      "value":[
         "COSM4860838"
      ]
   }
}
```

```Shell
curl 'https://api.gdc.cancer.gov/ssms?pretty=true&filters=%7B%0A%22op%22%3A%22in%22%2C%0A%22content%22%3A%7B%0A%22field%22%3A%22cosmic_id%22%2C%0A%22value%22%3A%5B%0A%22COSM4860838%22%0A%5D%0A%7D%0A%7D%0A'
```

```Response
{
  "data": {
    "hits": [
      {
        "ncbi_build": "GRCh38",
        "mutation_type": "Simple Somatic Mutation",
        "mutation_subtype": "Single base substitution",
        "end_position": 62438203,
        "reference_allele": "C",
        "ssm_id": "8b3c1a7a-e4e0-5200-9d46-5767c2982145",
        "start_position": 62438203,
        "cosmic_id": [
          "COSM4860838",
          "COSM731764",
          "COSM731765"
        ],
        "id": "8b3c1a7a-e4e0-5200-9d46-5767c2982145",
        "tumor_allele": "T",
        "gene_aa_change": [
          "CADPS G1147G",
          "CADPS G1187G",
          "CADPS G1217G",
          "CADPS G1226G",
          "CADPS G127G",
          "CADPS G218G",
          "CADPS G95G"
        ],
        "chromosome": "chr3",
        "genomic_dna_change": "chr3:g.62438203C>T"
      }
    ],
    "pagination": {
      "count": 1,
      "sort": "",
      "from": 0,
      "page": 1,
      "total": 1,
      "pages": 1,
      "size": 10
    }
  },
  "warnings": {}
}
```

## Simple Somatic Mutation Occurrence Endpoint Examples

__Example 1:__ A user wants to determine the chromosome in case `TCGA-DU-6407` that contains the greatest number of `ssms`. As this relates to mutations that are observed in a case, the `ssm_occurrences` endpoint is used.

```
{  
   "op":"in",
   "content":{  
      "field":"case.submitter_id",
      "value":["TCGA-DU-6407"]
   }
}
```
```Shell
curl "https://api.gdc.cancer.gov/ssm_occurrences?format=tsv&fields=ssm.chromosome&size=5000&filters=%7B%0D%0A%22op%22%3A%22in%22%2C%0D%0A%22content%22%3A%7B%0D%0A%22field%22%3A%22case.submitter_id%22%2C%0D%0A%22value%22%3A%5B%0D%0A%22TCGA-DU-6407%22%0D%0A%5D%0D%0A%7D%0D%0A%7D"
```
```Response
ssm.chromosome	id
chr3	552c09d1-69b1-5c04-b543-524a6feae3eb
chr10	391011ff-c1fd-5e2a-a128-652bc660f64c
chr10	1378cbc4-af88-55bb-b2e5-185bb4246d7a
chr10	3a2b3870-a395-5bc3-8c8f-0d40b0f2202c
chr1	4a93d7a5-988d-5055-80da-999dc3b45d80
chrX	22a07c7c-16ba-51df-a9a9-1e41e2a45225
chr12	dbc5eafa-ea26-5f1c-946c-b6974a345b69
chr11	02ae553d-1f27-565d-96c5-2c3cfca7264a
chr2	faee73a9-4804-58ea-a91f-18c3d901774f
chr6	97c5b38b-fc96-57f5-8517-cc702b3aa70a
chr17	0010a89d-9434-5d97-8672-36ee394767d0
chr19	f08dcc53-eadc-5ceb-bf31-f6b38629e4cb
chrX	19ca262d-b354-54a0-b582-c4719e37e91d
chr19	c44a93a1-5c73-5cff-b40e-98ce7e5fe57b
chr3	b67f31b5-0341-518e-8fcc-811cd2e36af1
chr1	94abd5fd-d539-5a4a-8719-9615cf7cec5d
chr17	1476a543-2951-5ec4-b165-67551b47d810
chr2	b4822fc9-f0cc-56fd-9d97-f916234e309d
chr2	3548ecfe-5186-51e7-8f40-37f4654cd260
chr16	105e7811-4601-5ccb-ae93-e7107923599e
chr2	99b3aad4-d368-506d-99d6-047cbe5dff0f
chr13	9dc3f7cd-9efa-530a-8524-30d067e49d54
chr21	1267330b-ae6d-5e25-b19e-34e98523679e
chr16	c77f7ce5-fbe6-5da4-9a7b-b528f8e530cb
chr10	2cb06277-993e-5502-b2c5-263037c45d18
chr17	d25129ad-3ad7-584f-bdeb-fba5c3881d32
chr17	a76469cb-973c-5d4d-bf82-7cf4e8f6c129
chr10	727c9d57-7b74-556f-aa5b-e1ca1f76d119
chr15	b4a86ffd-e60c-5c9c-aaa1-9e9f02d86116
chr5	3a023e72-da92-54f7-aa18-502c1076b2b0
```

## Analysis Endpoints

In addition to the `ssms`, `ssm_occurrences`, and `genes` endpoints mentioned previously, several `/analysis` endpoints were designed to quickly retrieve specific datasets used for visualization display.  

__Example 1:__ The `/analysis/top_cases_counts_by_genes` endpoint gives the number of cases with a mutation in each gene listed in the `gene_ids` parameter for each project. Note that this endpoint cannot be used with the `format` or `fields` parameters. In this instance, the query will produce the number of cases in each projects with mutations in the gene `ENSG00000155657`.

```Shell
curl "https://api.gdc.cancer.gov/analysis/top_cases_counts_by_genes?gene_ids=ENSG00000155657&pretty=true"
```


This JSON-formatted output is broken up by project. For an example, see the following text:

```json
          "genes": {
            "my_genes": {
              "gene_id": {
                "buckets": [
                  {
                    "key": "ENSG00000155657",
                    "doc_count": 45
                  }
                ],
                "sum_other_doc_count": 0,
                "doc_count_error_upper_bound": 0
              },
              "doc_count": 45
            },
            "doc_count": 12305
          },
          "key": "TCGA-GBM",
          "doc_count": 45
        }
```

This portion of the output shows TCGA-GBM including 45 cases that have `ssms` in the gene `ENSG00000155657`.

__Example 2:__ The following demonstrates a use of the `/analysis/top_mutated_genes_by_project` endpoint.  This will output the genes that are mutated in the most cases in "TCGA-DLBC" and will count the mutations that have a `HIGH` or `MODERATE` impact on gene function. Note that the `score` field does not represent the number of mutations in a given gene, but a calculation that is used to determine which genes have the greatest number of unique mutations.  

```json
{  
   "op":"AND",
   "content":[  
      {  
         "op":"in",
         "content":{  
            "field":"case.project.project_id",
            "value":[  
               "TCGA-DLBC"
            ]
         }
      },
      {  
         "op":"in",
         "content":{  
            "field":"case.ssm.consequence.transcript.annotation.impact",
            "value":[  
               "HIGH",
               "MODERATE"
            ]
         }
      }
   ]
}
```
```Shell
curl "https://api.gdc.cancer.gov/analysis/top_mutated_genes_by_project?fields=gene_id,symbol&filters=%7B%22op%22%3A%22AND%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22case.project.project_id%22%2C%22value%22%3A%5B%22TCGA-DLBC%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22case.ssm.consequence.transcript.annotation.impact%22%2C%22value%22%3A%5B%22HIGH%22%2C%22MODERATE%22%5D%7D%7D%5D%7D&pretty=true"
```
```Response
{
  "data": {
    "hits": [
      {
        "_score": 14.0,
        "symbol": "IGHV2-70",
        "gene_id": "ENSG00000274576"
      },
      {
        "_score": 14.0,
        "symbol": "IGLV3-1",
        "gene_id": "ENSG00000211673"
      },
      {
        "_score": 14.0,
        "symbol": "IGHM",
        "gene_id": "ENSG00000211899"
      },
      {
        "_score": 11.0,
        "symbol": "KMT2D",
        "gene_id": "ENSG00000167548"
      },
      {
        "_score": 11.0,
        "symbol": "IGLL5",
        "gene_id": "ENSG00000254709"
      },
      {
        "_score": 11.0,
        "symbol": "BTG2",
        "gene_id": "ENSG00000159388"
      },
      {
        "_score": 9.0,
        "symbol": "CARD11",
        "gene_id": "ENSG00000198286"
      },
      {
        "_score": 9.0,
        "symbol": "IGHG1",
        "gene_id": "ENSG00000211896"
      },
      {
        "_score": 9.0,
        "symbol": "IGLC2",
        "gene_id": "ENSG00000211677"
      },
      {
        "_score": 9.0,
        "symbol": "LRP1B",
        "gene_id": "ENSG00000168702"
      }
    ],
    "pagination": {
      "count": 10,
      "sort": "None",
      "from": 0,
      "page": 1,
      "total": 3214,
      "pages": 322,
      "size": 10
    }
  },
  "warnings": {}
}
```

__Example 3:__ The `/analysis/top_mutated_cases_by_gene` endpoint will generate information about the cases that are most affected by mutations in a given number of genes. Below, the file count for each category is given for the cases most affected by mutations in these 50 genes.  The size of the output is limited to two cases with the `size=2` parameter, but a higher value can be set by the user.

```Shell
curl "https://api.gdc.cancer.gov/analysis/top_mutated_cases_by_gene?fields=diagnoses.days_to_death,diagnoses.age_at_diagnosis,diagnoses.vital_status,diagnoses.primary_diagnosis,demographic.gender,demographic.race,demographic.ethnicity,case_id,summary.data_categories.file_count,summary.data_categories.data_category&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%22TCGA-DLBC%22%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22genes.gene_id%22%2C%22value%22%3A%5B%22ENSG00000166710%22%2C%22ENSG00000005339%22%2C%22ENSG00000083857%22%2C%22ENSG00000168769%22%2C%22ENSG00000100906%22%2C%22ENSG00000184677%22%2C%22ENSG00000101680%22%2C%22ENSG00000101266%22%2C%22ENSG00000028277%22%2C%22ENSG00000140968%22%2C%22ENSG00000181827%22%2C%22ENSG00000116815%22%2C%22ENSG00000275221%22%2C%22ENSG00000139083%22%2C%22ENSG00000112851%22%2C%22ENSG00000112697%22%2C%22ENSG00000164134%22%2C%22ENSG00000009413%22%2C%22ENSG00000071626%22%2C%22ENSG00000135407%22%2C%22ENSG00000101825%22%2C%22ENSG00000104814%22%2C%22ENSG00000166415%22%2C%22ENSG00000142867%22%2C%22ENSG00000254585%22%2C%22ENSG00000139718%22%2C%22ENSG00000077721%22%2C%22ENSG00000130294%22%2C%22ENSG00000117245%22%2C%22ENSG00000117318%22%2C%22ENSG00000270550%22%2C%22ENSG00000163637%22%2C%22ENSG00000166575%22%2C%22ENSG00000065526%22%2C%22ENSG00000156453%22%2C%22ENSG00000128191%22%2C%22ENSG00000055609%22%2C%22ENSG00000204469%22%2C%22ENSG00000187605%22%2C%22ENSG00000185875%22%2C%22ENSG00000110888%22%2C%22ENSG00000007341%22%2C%22ENSG00000173198%22%2C%22ENSG00000115568%22%2C%22ENSG00000163714%22%2C%22ENSG00000125772%22%2C%22ENSG00000080815%22%2C%22ENSG00000189079%22%2C%22ENSG00000120837%22%2C%22ENSG00000143951%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22ssms.consequence.transcript.annotation.impact%22%2C%22value%22%3A%5B%22HIGH%22%5D%7D%7D%5D%7D&pretty=true&size=2"
```
```Response
{
  "data": {
    "hits": [
      {
        "_score": 7.0,
        "diagnoses": [
          {
            "days_to_death": null,
            "vital_status": "alive",
            "age_at_diagnosis": 18691,
            "primary_diagnosis": "c83.3"
          }
        ],
        "case_id": "eda9496e-be80-4a13-bf06-89f0cc9e937f",
        "demographic": {
          "gender": "male",
          "race": "white",
          "ethnicity": "hispanic or latino"
        },
        "summary": {
          "data_categories": [
            {
              "file_count": 1,
              "data_category": "DNA Methylation"
            },
            {
              "file_count": 5,
              "data_category": "Transcriptome Profiling"
            },
            {
              "file_count": 1,
              "data_category": "Biospecimen"
            },
            {
              "file_count": 16,
              "data_category": "Simple Nucleotide Variation"
            },
            {
              "file_count": 1,
              "data_category": "Clinical"
            },
            {
              "file_count": 4,
              "data_category": "Copy Number Variation"
            },
            {
              "file_count": 4,
              "data_category": "Raw Sequencing Data"
            }
          ]
        }
      },
      {
        "_score": 4.0,
        "diagnoses": [
          {
            "days_to_death": null,
            "vital_status": "alive",
            "age_at_diagnosis": 27468,
            "primary_diagnosis": "c83.3"
          }
        ],
        "case_id": "a43e5f0e-a21f-48d8-97e0-084d413680b7",
        "demographic": {
          "gender": "male",
          "race": "white",
          "ethnicity": "not hispanic or latino"
        },
        "summary": {
          "data_categories": [
            {
              "file_count": 1,
              "data_category": "DNA Methylation"
            },
            {
              "file_count": 5,
              "data_category": "Transcriptome Profiling"
            },
            {
              "file_count": 1,
              "data_category": "Biospecimen"
            },
            {
              "file_count": 16,
              "data_category": "Simple Nucleotide Variation"
            },
            {
              "file_count": 1,
              "data_category": "Clinical"
            },
            {
              "file_count": 4,
              "data_category": "Copy Number Variation"
            },
            {
              "file_count": 4,
              "data_category": "Raw Sequencing Data"
            }
          ]
        }
      }
    ],
    "pagination": {
      "count": 2,
      "sort": "None",
      "from": 0,
      "page": 1,
      "total": 27,
      "pages": 14,
      "size": 2
    }
  },
  "warnings": {}
}
```

__Example 4:__  The `/analysis/mutated_cases_count_by_project` endpoint produces counts for the number of cases that have associated `ssm` data in each project. The number of affected cases can be found under `"case_with_ssm": {"doc_count": $case_count}`.  

```Shell
curl "https://api.gdc.cancer.gov/analysis/mutated_cases_count_by_project?size=0&pretty=true"
```
```Response
{
  "hits": {
    "hits": [],
    "total": 14551,
    "max_score": 0.0
  },
  "_shards": {
    "successful": 9,
    "failed": 0,
    "total": 9
  },
  "took": 4,
  "aggregations": {
    "projects": {
      "buckets": [
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 216
            },
            "doc_count": 637
          },
          "key": "TARGET-NBL",
          "doc_count": 1127
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 1044
            },
            "doc_count": 7625
          },
          "key": "TCGA-BRCA",
          "doc_count": 1098
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 8
            },
            "doc_count": 579
          },
          "key": "TARGET-AML",
          "doc_count": 988
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 34
            },
            "doc_count": 290
          },
          "key": "TARGET-WT",
          "doc_count": 652
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 396
            },
            "doc_count": 3197
          },
          "key": "TCGA-GBM",
          "doc_count": 617
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 443
            },
            "doc_count": 3880
          },
          "key": "TCGA-OV",
          "doc_count": 608
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 569
            },
            "doc_count": 3874
          },
          "key": "TCGA-LUAD",
          "doc_count": 585
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 542
            },
            "doc_count": 3874
          },
          "key": "TCGA-UCEC",
          "doc_count": 560
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 339
            },
            "doc_count": 3547
          },
          "key": "TCGA-KIRC",
          "doc_count": 537
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 510
            },
            "doc_count": 3671
          },
          "key": "TCGA-HNSC",
          "doc_count": 528
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 513
            },
            "doc_count": 3606
          },
          "key": "TCGA-LGG",
          "doc_count": 516
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 496
            },
            "doc_count": 3536
          },
          "key": "TCGA-THCA",
          "doc_count": 507
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 497
            },
            "doc_count": 3520
          },
          "key": "TCGA-LUSC",
          "doc_count": 504
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 498
            },
            "doc_count": 3490
          },
          "key": "TCGA-PRAD",
          "doc_count": 500
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 470
            },
            "doc_count": 3289
          },
          "key": "TCGA-SKCM",
          "doc_count": 470
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 433
            },
            "doc_count": 3188
          },
          "key": "TCGA-COAD",
          "doc_count": 461
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 441
            },
            "doc_count": 3095
          },
          "key": "TCGA-STAD",
          "doc_count": 443
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 412
            },
            "doc_count": 2884
          },
          "key": "TCGA-BLCA",
          "doc_count": 412
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 0
            },
            "doc_count": 0
          },
          "key": "TARGET-OS",
          "doc_count": 381
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 375
            },
            "doc_count": 2635
          },
          "key": "TCGA-LIHC",
          "doc_count": 377
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 305
            },
            "doc_count": 2142
          },
          "key": "TCGA-CESC",
          "doc_count": 307
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 288
            },
            "doc_count": 2033
          },
          "key": "TCGA-KIRP",
          "doc_count": 291
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 255
            },
            "doc_count": 1821
          },
          "key": "TCGA-SARC",
          "doc_count": 261
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 149
            },
            "doc_count": 1192
          },
          "key": "TCGA-LAML",
          "doc_count": 200
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 184
            },
            "doc_count": 1293
          },
          "key": "TCGA-ESCA",
          "doc_count": 185
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 183
            },
            "doc_count": 1285
          },
          "key": "TCGA-PAAD",
          "doc_count": 185
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 179
            },
            "doc_count": 1253
          },
          "key": "TCGA-PCPG",
          "doc_count": 179
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 158
            },
            "doc_count": 1169
          },
          "key": "TCGA-READ",
          "doc_count": 172
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 150
            },
            "doc_count": 1018
          },
          "key": "TCGA-TGCT",
          "doc_count": 150
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 123
            },
            "doc_count": 867
          },
          "key": "TCGA-THYM",
          "doc_count": 124
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 66
            },
            "doc_count": 556
          },
          "key": "TCGA-KICH",
          "doc_count": 113
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 92
            },
            "doc_count": 620
          },
          "key": "TCGA-ACC",
          "doc_count": 92
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 83
            },
            "doc_count": 605
          },
          "key": "TCGA-MESO",
          "doc_count": 87
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 80
            },
            "doc_count": 560
          },
          "key": "TCGA-UVM",
          "doc_count": 80
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 0
            },
            "doc_count": 163
          },
          "key": "TARGET-RT",
          "doc_count": 75
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 48
            },
            "doc_count": 346
          },
          "key": "TCGA-DLBC",
          "doc_count": 58
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 57
            },
            "doc_count": 399
          },
          "key": "TCGA-UCS",
          "doc_count": 57
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 51
            },
            "doc_count": 306
          },
          "key": "TCGA-CHOL",
          "doc_count": 51
        },
        {
          "case_summary": {
            "case_with_ssm": {
              "doc_count": 0
            },
            "doc_count": 13
          },
          "key": "TARGET-CCSK",
          "doc_count": 13
        }
      ],
      "sum_other_doc_count": 0,
      "doc_count_error_upper_bound": 0
    }
  },
  "timed_out": false
}
```
### Survival Analysis Endpoint

[Survival plots](/Data_Portal/Projects/#Survival-Analysis) are generated for different subsets of data, based on variants or projects, in the GDC Data Portal. The `/analysis/survival` endpoint can be used to programmatically retrieve the raw data used to generate these plots and apply different filters. Note that the `fields` and `format` parameters cannot be modified.

 __Example 1:__ A user wants to download data to generate a survival plot for cases from the project TCGA-DLBC.

```Shell
curl "https://api.gdc.cancer.gov/analysis/survival?filters=%5B%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%22TCGA-DLBC%22%7D%7D%5D&pretty=true"
```
```Response
{
  "overallStats": {},
  "results": [
    {
      "donors": [
        {
          "survivalEstimate": 1,
          "id": "dc87a809-95de-4eb7-a1c2-2650475f2d7e",
          "censored": true,
          "time": 1
        },
        {
          "survivalEstimate": 1,
          "id": "4dd86ebd-ef16-4b2b-9ea0-5d1d7afef257",
          "censored": true,
          "time": 17
        },
        {
          "survivalEstimate": 1,
          "id": "0bf573ac-cd1e-42d8-90cf-b30d7b08679c",
          "censored": false,
          "time": 58
        },
        {
          "survivalEstimate": 0.9777777777777777,
          "id": "f978cb0f-d319-4c01-b4c5-23ae1403a106",
          "censored": true,
          "time": 126
        },
        {
          "survivalEstimate": 0.9777777777777777,
          "id": "a43e5f0e-a21f-48d8-97e0-084d413680b7",
          "censored": true,
          "time": 132
        },
        {
          "survivalEstimate": 0.9777777777777777,
          "id": "1843c82e-7a35-474f-9f79-c0a9af9aa09c",
          "censored": true,
          "time": 132
        },
        {
          "survivalEstimate": 0.9777777777777777,
          "id": "0030a28c-81aa-44b0-8be0-b35e1dcbf98c",
          "censored": false,
          "time": 248
        },
        {
          "survivalEstimate": 0.9539295392953929,
          "id": "f553f1a9-ecf2-4783-a609-6adca7c4c597",
          "censored": true,
          "time": 298
        },
        {
          "survivalEstimate": 0.9539295392953929,
          "id": "f784bc3a-751b-4025-aab2-0af2f6f24266",
          "censored": false,
          "time": 313
        },
        {
          "survivalEstimate": 0.929469807518588,
          "id": "29e3d122-15a1-4235-a356-b1a9f94ceb39",
          "censored": true,
          "time": 385
        },
        {
          "survivalEstimate": 0.929469807518588,
          "id": "0e251c03-bf86-4ed8-b45d-3cbc97160502",
          "censored": false,
          "time": 391
        },
        {
          "survivalEstimate": 0.9043490019099776,
          "id": "e6365b38-bc44-400c-b4aa-18ce8ff5bfce",
          "censored": true,
          "time": 427
        },
        {
          "survivalEstimate": 0.9043490019099776,
          "id": "b56bdbdb-43af-4a03-a072-54dd22d7550c",
          "censored": true,
          "time": 553
        },
        {
          "survivalEstimate": 0.9043490019099776,
          "id": "31bbad4e-3789-42ec-9faa-1cb86970f723",
          "censored": false,
          "time": 595
        },
        {
          "survivalEstimate": 0.8777505018538018,
          "id": "0e9fcccc-0630-408d-a121-2c6413824cb7",
          "censored": true,
          "time": 679
        },
        {
          "survivalEstimate": 0.8777505018538018,
          "id": "a5b188f0-a6d3-4d4a-b04f-36d47ec05338",
          "censored": false,
          "time": 708
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "ed746cb9-0f2f-48ce-923a-3a9f9f00b331",
          "censored": true,
          "time": 719
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "c85f340e-584b-4f3b-b6a5-540491fc8ad2",
          "censored": true,
          "time": 730
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "69f23725-adca-48ac-9b33-80a7aae24cfe",
          "censored": true,
          "time": 749
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "67325322-483f-443f-9ffa-2a20d108a2fb",
          "censored": true,
          "time": 751
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "eda9496e-be80-4a13-bf06-89f0cc9e937f",
          "censored": true,
          "time": 765
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "25ff86af-beb4-480c-b706-f3fe0306f7cf",
          "censored": true,
          "time": 788
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "1d0db5d7-39ca-466d-96b3-0d278c5ea768",
          "censored": true,
          "time": 791
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "c8cde9ea-89e9-4ee8-8a46-417a48f6d3ab",
          "censored": true,
          "time": 832
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "f0a326d2-1f3e-4a5d-bca8-32aaccc52338",
          "censored": true,
          "time": 946
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "a8e2df1e-4042-42af-9231-3a00e83489f0",
          "censored": true,
          "time": 965
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "e56e4d9c-052e-4ec6-a81b-dbd53e9c8ffe",
          "censored": true,
          "time": 972
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "45b0cf9f-a879-417f-8f39-7770552252c0",
          "censored": true,
          "time": 982
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "1f971af1-6772-4fe6-8d35-bbe527a037fe",
          "censored": true,
          "time": 1081
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "33365d22-cb83-4d8e-a2d1-06b675f75f6e",
          "censored": true,
          "time": 1163
        },
        {
          "survivalEstimate": 0.8503207986708705,
          "id": "6a21c948-cd85-4150-8c01-83017d7dc1ed",
          "censored": false,
          "time": 1252
        },
        {
          "survivalEstimate": 0.8003019281608192,
          "id": "f855dad1-6ffc-493e-ba6c-970874bc9210",
          "censored": true,
          "time": 1299
        },
        {
          "survivalEstimate": 0.8003019281608192,
          "id": "c1c06604-5ae2-4a53-b9c0-eb210d38e3f0",
          "censored": true,
          "time": 1334
        },
        {
          "survivalEstimate": 0.8003019281608192,
          "id": "58e66976-4507-4552-ac53-83a49a142dde",
          "censored": true,
          "time": 1373
        },
        {
          "survivalEstimate": 0.8003019281608192,
          "id": "ea54dbad-1b23-41cc-9378-d4002a8fca51",
          "censored": true,
          "time": 1581
        },
        {
          "survivalEstimate": 0.8003019281608192,
          "id": "d7df78b5-24f1-4ff4-bd9b-f0e6bec8289a",
          "censored": true,
          "time": 1581
        },
        {
          "survivalEstimate": 0.8003019281608192,
          "id": "29aff186-c321-4ff9-b81b-105e27e620ff",
          "censored": true,
          "time": 1617
        },
        {
          "survivalEstimate": 0.8003019281608192,
          "id": "5eff68ff-f6c3-40c9-9fc8-00e684a7b712",
          "censored": true,
          "time": 1739
        },
        {
          "survivalEstimate": 0.8003019281608192,
          "id": "f8cf647b-1447-4ac3-8c43-bef07765cabf",
          "censored": true,
          "time": 2131
        },
        {
          "survivalEstimate": 0.8003019281608192,
          "id": "c3d662ee-48d0-454a-bb0c-77d3338d3747",
          "censored": true,
          "time": 2983
        },
        {
          "survivalEstimate": 0.8003019281608192,
          "id": "6e9437f0-a4ed-475c-ab0e-bf1431c70a90",
          "censored": true,
          "time": 3333
        },
        {
          "survivalEstimate": 0.8003019281608192,
          "id": "fdecb74f-ac4e-46b1-b23a-5f7fde96ef9f",
          "censored": true,
          "time": 3394
        },
        {
          "survivalEstimate": 0.8003019281608192,
          "id": "a468e725-ad4b-411d-ac5c-2eacc68ec580",
          "censored": false,
          "time": 3553
        },
        {
          "survivalEstimate": 0.6402415425286554,
          "id": "1ea575f1-f731-408b-a629-f5f4abab569e",
          "censored": true,
          "time": 3897
        },
        {
          "survivalEstimate": 0.6402415425286554,
          "id": "7a589441-11ef-4158-87e7-3951d86bc2aa",
          "censored": true,
          "time": 4578
        },
        {
          "survivalEstimate": 0.6402415425286554,
          "id": "3622cf29-600f-4410-84d4-a9afeb41c475",
          "censored": true,
          "time": 5980
        },
        {
          "survivalEstimate": 0.6402415425286554,
          "id": "3f5a897d-1eaa-4d4c-8324-27ac07c90927",
          "censored": false,
          "time": 6425
        }
      ],
      "meta": {
        "id": 140429063094496
      }
    }
  ]
}
```

__Example 2:__ Here the survival endpoint is used to compare two survival plots for TCGA-BRCA cases.  One plot will display survival information about cases with a particular mutation (in this instance: `chr3:g.179234297A>G`) and the other plot will display information about cases without that mutation. This type of query will also print the results of a chi-squared analysis between the two subsets of cases.  

```json
[  
   {  
      "op":"and",
      "content":[  
         {  
            "op":"=",
            "content":{  
               "field":"cases.project.project_id",
               "value":"TCGA-BRCA"
            }
         },
         {  
            "op":"=",
            "content":{  
               "field":"gene.ssm.ssm_id",
               "value":"937a26c2-089c-51de-a0f9-70567d965c38"
            }
         }
      ]
   },
   {  
      "op":"and",
      "content":[  
         {  
            "op":"=",
            "content":{  
               "field":"cases.project.project_id",
               "value":"TCGA-BRCA"
            }
         },
         {  
            "op":"excludeifany",
            "content":{  
               "field":"gene.ssm.ssm_id",
               "value":"937a26c2-089c-51de-a0f9-70567d965c38"
            }
         }
      ]
   }
]
```
```Shell
curl "https://api.gdc.cancer.gov/analysis/survival?filters=%5B%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%22TCGA-BRCA%22%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22gene.ssm.ssm_id%22%2C%22value%22%3A%22937a26c2-089c-51de-a0f9-70567d965c38%22%7D%7D%5D%7D%2C%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%22TCGA-BRCA%22%7D%7D%2C%7B%22op%22%3A%22excludeifany%22%2C%22content%22%3A%7B%22field%22%3A%22gene.ssm.ssm_id%22%2C%22value%22%3A%22937a26c2-089c-51de-a0f9-70567d965c38%22%7D%7D%5D%7D%5D&pretty=true"
```
```json2
{
  "overallStats": {
    "degreesFreedom": 1,
    "chiSquared": 0.8577589072612264,
    "pValue": 0.35436660628146011
  },
  "results": [
    {
      "donors": [
        {
          "survivalEstimate": 1,
          "id": "a991644b-3ee6-4cda-acf0-e37de48a49fc",
          "censored": true,
          "time": 10
        },
        {
          "survivalEstimate": 1,
          "id": "2e1e3bf0-1708-4b65-936c-48b89eb8966a",
          "censored": true,
          "time": 19
        },
(truncated)
],
"meta": {
  "id": 140055251282040
}
},
{
"donors": [
  {
    "survivalEstimate": 1,
    "id": "5e4187c9-98f8-4bdb-a8da-6a914e96f47a",
    "censored": true,
    "time": -31
  },
(truncated)
```

The output represents two sets of coordinates delimited as objects with the `donors` tag. One set of coordinates will generate a survival plot representing TCGA-BRCA cases that have the mutation of interest and the other will generate a survival plot for the remaining cases in TCGA-BRCA.

__Example 3:__ Custom survival plots can be generated using the GDC API.  For example, a user could generate survival plot data comparing patients with a mutation in genes associated with a biological pathway with patients without mutations in that pathway.  

``` Query
[  
   {  
      "op":"and",
      "content":[  
         {  
            "op":"=",
            "content":{  
               "field":"cases.project.project_id",
               "value":"TCGA-BRCA"
            }
         },
         {  
            "op":"=",
            "content":{  
               "field":"gene.gene_id",
               "value":["ENSG00000141510","ENSG00000155657"]
            }
         }
      ]
   },
   {  
      "op":"and",
      "content":[  
         {  
            "op":"=",
            "content":{  
               "field":"cases.project.project_id",
               "value":"TCGA-BRCA"
            }
         },
         {  
            "op":"excludeifany",
            "content":{  
               "field":"gene.gene_id",
               "value":["ENSG00000141510","ENSG00000155657"]
            }
         }
      ]
   }
]
```
```Shell
curl "https://api.gdc.cancer.gov/analysis/survival?filters=%5B%0D%0A%7B%0D%0A%22op%22%3A%22and%22%2C%0D%0A%22content%22%3A%5B%0D%0A%7B%0D%0A%22op%22%3A%22%3D%22%2C%0D%0A%22content%22%3A%7B%0D%0A%22field%22%3A%22cases.project.project_id%22%2C%0D%0A%22value%22%3A%22TCGA-BRCA%22%0D%0A%7D%0D%0A%7D%2C%0D%0A%7B%0D%0A%22op%22%3A%22%3D%22%2C%0D%0A%22content%22%3A%7B%0D%0A%22field%22%3A%22gene.gene_id%22%2C%0D%0A%22value%22%3A%5B%22ENSG00000141510%22%2C%22ENSG00000155657%22%5D%0D%0A%7D%0D%0A%7D%0D%0A%5D%0D%0A%7D%2C%0D%0A%7B%0D%0A%22op%22%3A%22and%22%2C%0D%0A%22content%22%3A%5B%0D%0A%7B%0D%0A%22op%22%3A%22%3D%22%2C%0D%0A%22content%22%3A%7B%0D%0A%22field%22%3A%22cases.project.project_id%22%2C%0D%0A%22value%22%3A%22TCGA-BRCA%22%0D%0A%7D%0D%0A%7D%2C%0D%0A%7B%0D%0A%22op%22%3A%22excludeifany%22%2C%0D%0A%22content%22%3A%7B%0D%0A%22field%22%3A%22gene.gene_id%22%2C%0D%0A%22value%22%3A%5B%22ENSG00000141510%22%2C%22ENSG00000155657%22%5D%0D%0A%7D%0D%0A%7D%0D%0A%5D%0D%0A%7D%0D%0A%5D&pretty=true"
```
