# Data Analysis

The GDC DAVE tools use the same API as the rest of the Data Portal and takes advantage of several new endpoints. Similar to the [GDC Data Portal Exploration](http://docs.gdc.cancer.gov/Data_Portal/Users_Guide/Exploration/) feature, the GDC data analysis endpoints allow API users to programmatically explore data in the GDC using advanced filters at a gene and mutation level. Survival analysis data is also available.  

## Endpoints

The following data analysis endpoints are available from the GDC API:

|__Node__| __Endpoint__ | __Description__ |
|---|---|---|
|__Genes__| __/genes__ | Allows users to access summary information about each gene using its Ensembl ID. |
|__SSMS__| __/ssms__ | Allows users to access information about each somatic mutation. For example, a `ssm` would represent the transition of C to T at position 52000 of chromosome 1. |
||__/ssms/`<ssm_id>`__|Get information about a specific ssm using a `<ssm_id>`, often supplemented with the `expand` option to show fields of interest. |
|| __/ssm_occurrences__ | A `ssm` entity as applied to a single instance (case). An example of a `ssm occurrence` would be that the transition of C to T at position 52000 of chromosome 1 occurred in patient TCGA-XX-XXXX. |
||__/ssm_occurrences/`<ssm_occurrences_id>`__|Get information about a specific ssm occurrence using a `<ssm_occurrences_id>`, often supplemented with the `expand` option to show fields of interest. |
|__CNVS__|__/cnvs__|Allows users to access data about copy number variations (cnvs). This data will be specifc to cnvs and not a specific case. |
||__/cnvs/`<cnv_id>`__|Get information about a specific copy number variation using a `<cnv_id>`, often supplemented with the `expand` option to show fields of interest. |
||__/cnvs/ids__|This endpoint will retrieve nodes that contain the queried cnv_id. This is accomplished by adding the query parameter: /cnvs/ids?query=`<cnv_id>`.|
||__/cnv_occurrences__|A `cnv` entity as applied to a single case.|
||__/cnv_occurrences/`<cnv_occurrence_id>`__|Get information about a specific copy number variation occurrence using a `<cnv_occurrence_id>`, often supplemented with the `expand` option to show fields of interest. |
||__/cnv_occurrences/ids__|This endpoint will retrieve nodes that contain the queried cnv_occurrence_id. This is accomplished by adding the query parameter: /cnv_occurrences/ids?query=`<cnv_occurrences_id>`|
|__Analysis__|__/analysis/top_cases_counts_by_genes__| Returns the number of cases with a mutation in each gene listed in the gene_ids parameter for each project. Note that this endpoint cannot be used with the `format` or `fields` parameters.|
||__/analysis/top_mutated_genes_by_project__| Returns a list of genes that have the most mutations within a given project. |
||__/analysis/top_mutated_cases_by_gene__| Generates information about the cases that are most affected by mutations in a given number of genes |
||__/analysis/mutated_cases_count_by_project__| Returns counts for the number of cases that have associated `ssm` data in each project. The number of affected cases can be found under "case_with_ssm": {"doc_count": $case_count}.|
||__/analysis/survival__| Survival plots can be generated in the Data Portal for different subsets of data, based upon many query factors such as variants, disease type and projects. This endpoint can be used to programmatically retrieve the raw data to generate these plots and apply different filters to the data. (see Survival Example)|


The methods for retrieving information from these endpoints are very similar to those used for the `cases` and `files` endpoints. These methods are explored in depth in the [API Search and Retrieval](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/) documentation. The `_mapping` parameter can also be used with each of these endpoints to generate a list of potential fields.  For example:

`https://api.gdc.cancer.gov/ssms/_mapping`

While it is not an endpoint, the `observation` entity is featured in the visualization section of the API. The `observation` entity provides information from the MAF file, such as read depth and normal genotype, that supports the validity of the associated `ssm`. An example is demonstrated below:

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

__Example 1__: Similar to the `/genes` endpoint, a user would like to retrieve information about the mutation based on its COSMIC ID. This would be accomplished by creating a JSON filter, which will then be encoded to URL for the `curl` command.

```Filter
 {
   "op":"in",
   "content":{
      "field":"cosmic_id",
      "value":[
         "COSM1135366"
      ]
   }
}
```

```Shell
curl 'https://api.gdc.cancer.gov/ssms?pretty=true&filters=%7B%0A%22op%22%3A%22in%22%2C%0A%22content%22%3A%7B%0A%22field%22%3A%22cosmic_id%22%2C%0A%2value%22%3A%5B%0A%22COSM1135366%22%0A%5D%0A%7D%0A%7D%0A'
```

```Response
{
  "data": {
    "hits": [
      {
        "id": "edd1ae2c-3ca9-52bd-a124-b09ed304fcc2",
        "start_position": 25245350,
        "gene_aa_change": [
          "KRAS G12D"
        ],
        "reference_allele": "C",
        "ncbi_build": "GRCh38",
        "cosmic_id": [
          "COSM1135366",
          "COSM521"
        ],
        "mutation_subtype": "Single base substitution",
        "mutation_type": "Simple Somatic Mutation",
        "chromosome": "chr12",
        "ssm_id": "edd1ae2c-3ca9-52bd-a124-b09ed304fcc2",
        "genomic_dna_change": "chr12:g.25245350C>T",
        "tumor_allele": "T",
        "end_position": 25245350
      }
    ],
    "pagination": {
      "count": 1,
      "total": 1,
      "size": 10,
      "from": 0,
      "sort": "",
      "page": 1,
      "pages": 1
    }
  },
  "warnings": {}
}
```

__Example 2:__ Based on the previous example's `ssm_id` (`8b3c1a7a-e4e0-5200-9d46-5767c2982145`), a user would like to look at the consequences and the VEP impact due to this ssm.

```Shell
curl 'https://api.gdc.cancer.gov/ssms/edd1ae2c-3ca9-52bd-a124-b09ed304fcc2?pretty=true&expand=consequence.transcript&fields=consequence.transcript.annotation.vep_impact'
```

```JSON
{
  "data": {
    "consequence": [
      {
        "transcript": {
          "annotation": {
            "vep_impact": "MODERATE"
          },
          "transcript_id": "ENST00000557334",
          "aa_end": 12,
          "consequence_type": "missense_variant",
          "aa_start": 12,
          "is_canonical": false,
          "aa_change": "G12D",
          "ref_seq_accession": ""
        }
      },
      {
        "transcript": {
          "annotation": {
            "vep_impact": "MODERATE"
          },
          "transcript_id": "ENST00000256078",
          "aa_end": 12,
          "consequence_type": "missense_variant",
          "aa_start": 12,
          "is_canonical": true,
          "aa_change": "G12D",
          "ref_seq_accession": "NM_001369786.1&NM_033360.4"
        }
      },
      {
        "transcript": {
          "annotation": {
            "vep_impact": "MODERATE"
          },
          "transcript_id": "ENST00000311936",
          "aa_end": 12,
          "consequence_type": "missense_variant",
          "aa_start": 12,
          "is_canonical": false,
          "aa_change": "G12D",
          "ref_seq_accession": "NM_001369787.1&NM_004985.5"
        }
      },
      {
        "transcript": {
          "annotation": {
            "vep_impact": "MODERATE"
          },
          "transcript_id": "ENST00000556131",
          "aa_end": 12,
          "consequence_type": "missense_variant",
          "aa_start": 12,
          "is_canonical": false,
          "aa_change": "G12D",
          "ref_seq_accession": ""
        }
      }
    ]
  },
  "warnings": {}
}
```

## Simple Somatic Mutation Occurrence Endpoint Examples

__Example 1:__ A user wants to determine the chromosome in case `TCGA-DU-6407` that contains the greatest number of `ssms`. As this relates to mutations that are observed in a case, the `ssm_occurrences` endpoint is used.

```Filter
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
```tsv
id	ssm.chromosome
105e7811-4601-5ccb-ae93-e7107923599e	chr16
faee73a9-4804-58ea-a91f-18c3d901774f	chr2
99b3aad4-d368-506d-99d6-047cbe5dff0f	chr2
2cb06277-993e-5502-b2c5-263037c45d18	chr10
f08dcc53-eadc-5ceb-bf31-f6b38629e4cb	chr19
97c5b38b-fc96-57f5-8517-cc702b3aa70a	chr6
19ca262d-b354-54a0-b582-c4719e37e91d	chrX
b4822fc9-f0cc-56fd-9d97-f916234e309d	chr2
22a07c7c-16ba-51df-a9a9-1e41e2a45225	chrX
0010a89d-9434-5d97-8672-36ee394767d0	chr17
3a023e72-da92-54f7-aa18-502c1076b2b0	chr5
391011ff-c1fd-5e2a-a128-652bc660f64c	chr10
3548ecfe-5186-51e7-8f40-37f4654cd260	chr2
b67f31b5-0341-518e-8fcc-811cd2e36af1	chr3
4a93d7a5-988d-5055-80da-999dc3b45d80	chr1
9dc3f7cd-9efa-530a-8524-30d067e49d54	chr13
552c09d1-69b1-5c04-b543-524a6feae3eb	chr3
dbc5eafa-ea26-5f1c-946c-b6974a345b69	chr12
d25129ad-3ad7-584f-bdeb-fba5c3881d32	chr17
1378cbc4-af88-55bb-b2e5-185bb4246d7a	chr10
c44a93a1-5c73-5cff-b40e-98ce7e5fe57b	chr19
1267330b-ae6d-5e25-b19e-34e98523679e	chr21
1476a543-2951-5ec4-b165-67551b47d810	chr17
727c9d57-7b74-556f-aa5b-e1ca1f76d119	chr10
94abd5fd-d539-5a4a-8719-9615cf7cec5d	chr1
a76469cb-973c-5d4d-bf82-7cf4e8f6c129	chr17
```
__Example 2:__ A user has retrieved a `ssm_occurrence`, and would like to determine if that case also has diagnostic information.

```Shell
curl 'https://api.gdc.cancer.gov/ssm_occurrences/6fd8527d-5c40-5604-8fa9-0ce798eec231?pretty=true&expand=case.diagnoses'
```

```Json
{
  "data": {
    "ssm_occurrence_id": "6fd8527d-5c40-5604-8fa9-0ce798eec231",
    "case": {
      "diagnoses": [
        {
          "ajcc_pathologic_t": "T3b",
          "synchronous_malignancy": "No",
          "morphology": "8720/3",
          "ajcc_pathologic_stage": "Stage IIB",
          "ajcc_pathologic_n": "N0",
          "ajcc_pathologic_m": "M0",
          "submitter_id": "TCGA-Z2-A8RT_diagnosis",
          "days_to_diagnosis": 0,
          "last_known_disease_status": "not reported",
          "tissue_or_organ_of_origin": "Skin, NOS",
          "days_to_last_follow_up": 839.0,
          "age_at_diagnosis": 15342,
          "primary_diagnosis": "Malignant melanoma, NOS",
          "classification_of_tumor": "not reported",
          "prior_malignancy": "no",
          "year_of_diagnosis": 2012,
          "diagnosis_id": "1d06a202-c51a-52e2-805f-eeb5f7fac14e",
          "icd_10_code": "C44.6",
          "site_of_resection_or_biopsy": "Skin of upper limb and shoulder",
          "prior_treatment": "No",
          "state": "released",
          "tumor_grade": "not reported",
          "progression_or_recurrence": "not reported",
          "ajcc_staging_system_edition": "7th"
        }
      ]
    }
  },
  "warnings": {}
}
```

## Copy Number Variation Endpoint Examples

__Example 1:__ A user is interested in finding the first 30 cnvs found on chromosome 4 that have a cnv loss.

```Filter
{
    "op": "and",
    "content": [
        {
            "op": "in",
            "content": {
                "field": "chromosome",
                "value": [
                    "4"
                ]
            }
        },
        {
            "op": "in",
            "content": {
                "field": "cnv_change",
                "value": [
                    "Loss"
                ]
            }
        }
    ]
}
```

```Shell
curl 'https://api.gdc.cancer.gov/cnvs?filters=%7B%0D%0A+++%22op%22%3A+%22and%22%2C%0D%0A++++%22content%22%3A+%5B%0D%0A++++++++%7B%0D%0A++++++++++++%22op%22%3A+%22in%22%2C%0D%0A++++++++++++%22content%22%3A+%7B%0D%0A++++++++++++++++%22field%22%3A+%22chromosome%22%2C%0D%0A++++++++++++++++%22value%22%3A+%5B%0D%0A++++++++++++++++++++%224%22%0D%0A++++++++++++++++%5D%0D%0A++++++++++++%7D%0D%0A++++++++%7D%2C%0D%0A++++++++%7B%0D%0A++++++++++++%22op%22%3A+%22in%22%2C%0D%0A++++++++++++%22content%22%3A+%7B%0D%0A++++++++++++++++%22field%22%3A+%22cnv_change%22%2C%0D%0A++++++++++++++++%22value%22%3A+%5B%0D%0A++++++++++++++++++++%22Loss%22%0D%0A++++++++++++++++%5D%0D%0A++++++++++++%7D%0D%0A++++++++%7D%0D%0A++++%5D%0D%0A%7D&size=30&sort=start_position&format=tsv'
```

```tsv
chromosome	cnv_change	cnv_id	end_position	gene_level_cn	id	ncbi_build	start_position
4	Loss	11381600-f064-5c42-90d2-a5c79c8b23e1	88208	True	11381600-f064-5c42-90d2-a5c79c8b23e1	GRCh38	53286
4	Loss	edef0f2f-c1a7-507c-842f-e1f8a568df9d	202303	True	edef0f2f-c1a7-507c-842f-e1f8a568df9d	GRCh38	124501
4	Loss	eba92f9a-b045-54a8-948a-451e439ed418	305474	True	eba92f9a-b045-54a8-948a-451e439ed418	GRCh38	270675
4	Loss	89319453-2a3f-5ebe-be30-8af0426e0343	384868	True	89319453-2a3f-5ebe-be30-8af0426e0343	GRCh38	337814
4	Loss	6567929c-4b6f-582b-aedf-acde2c0ec736	499156	True	6567929c-4b6f-582b-aedf-acde2c0ec736	GRCh38	425815
4	Loss	2daff58b-5065-50cd-8239-253180eaee81	540200	True	2daff58b-5065-50cd-8239-253180eaee81	GRCh38	499210
4	Loss	2b42c8d4-6d85-5352-96e1-9e52e722c248	576295	True	2b42c8d4-6d85-5352-96e1-9e52e722c248	GRCh38	573880
4	Loss	2646cdc7-7602-59a4-ae4f-d171352bae88	670782	True	2646cdc7-7602-59a4-ae4f-d171352bae88	GRCh38	625573
4	Loss	c11ad392-949f-593f-a3ab-d834b2f82809	674330	True	c11ad392-949f-593f-a3ab-d834b2f82809	GRCh38	672436
4	Loss	f31be658-4de0-549e-81be-e79759879acf	682033	True	f31be658-4de0-549e-81be-e79759879acf	GRCh38	673580
4	Loss	d72c62f2-fc29-5b83-9839-7f6b03970aff	689271	True	d72c62f2-fc29-5b83-9839-7f6b03970aff	GRCh38	681829
4	Loss	45448d47-6e13-5d30-824d-96150a7f55c6	770640	True	45448d47-6e13-5d30-824d-96150a7f55c6	GRCh38	705748
4	Loss	517e65ea-9084-54c2-abe0-b1b47e9f872c	826129	True	517e65ea-9084-54c2-abe0-b1b47e9f872c	GRCh38	784957
4	Loss	b5a09c9b-d842-5b76-a500-56f18252c29d	932373	True	b5a09c9b-d842-5b76-a500-56f18252c29d	GRCh38	849276
4	Loss	e3a3b61d-2881-5ad4-90bf-58ef29ae9ecb	958656	True	e3a3b61d-2881-5ad4-90bf-58ef29ae9ecb	GRCh38	932387
4	Loss	8630a1b6-3215-5b71-903a-ad9845505afc	986895	True	8630a1b6-3215-5b71-903a-ad9845505afc	GRCh38	958887
4	Loss	f748b06f-1fb7-53a9-a7d6-2c22a3ae6de5	993440	True	f748b06f-1fb7-53a9-a7d6-2c22a3ae6de5	GRCh38	979073
4	Loss	a5e4a63f-c5f6-5f0f-a6b6-f51bfb643533	1004564	True	a5e4a63f-c5f6-5f0f-a6b6-f51bfb643533	GRCh38	986997
4	Loss	73f6fbbe-6fd9-524c-a7c8-a7cf3f08ada4	1026898	True	73f6fbbe-6fd9-524c-a7c8-a7cf3f08ada4	GRCh38	1009936
4	Loss	adad579a-b002-5022-823a-570c59549065	1113564	True	adad579a-b002-5022-823a-570c59549065	GRCh38	1056250
4	Loss	d5a5c45e-594b-5cbc-97d5-75fc5155d021	1208962	True	d5a5c45e-594b-5cbc-97d5-75fc5155d021	GRCh38	1166932
4	Loss	6c910993-faa8-5abc-b433-b3afcc5e9e11	1249953	True	6c910993-faa8-5abc-b433-b3afcc5e9e11	GRCh38	1211445
4	Loss	4453b4cb-7d8a-5e26-a856-eac62eec287a	1340147	True	4453b4cb-7d8a-5e26-a856-eac62eec287a	GRCh38	1289887
4	Loss	6db1001a-a41b-518d-9491-2bf41544d90f	1395989	True	6db1001a-a41b-518d-9491-2bf41544d90f	GRCh38	1345691
4	Loss	6bef981a-ead1-5aa7-8a69-8d38e576e5c0	1406442	True	6bef981a-ead1-5aa7-8a69-8d38e576e5c0	GRCh38	1402932
4	Loss	af6e0b49-922a-587e-b353-4b9414605cf1	1684261	True	af6e0b49-922a-587e-b353-4b9414605cf1	GRCh38	1617915
4	Loss	400352ad-8526-562a-bbf4-29b90a48f46f	1712344	True	400352ad-8526-562a-bbf4-29b90a48f46f	GRCh38	1692731
4	Loss	8811414d-2434-56c6-afe5-a998c9b18d47	1745171	True	8811414d-2434-56c6-afe5-a998c9b18d47	GRCh38	1712891
4	Loss	169c4409-0256-5841-9314-f1a4dd2bcc38	1721358	True	169c4409-0256-5841-9314-f1a4dd2bcc38	GRCh38	1715952
4	Loss	1712ccac-6e70-5fb3-b71e-1a029eaf047c	1808872	True	1712ccac-6e70-5fb3-b71e-1a029eaf047c	GRCh38	1793293
```

__Example 2:__ A user wants to determine the location and identity of the gene affected by the cnv `544c4896-0152-5787-8d77-894a16f0ded0`, and determine whether the gene is found within the Cancer Gene Census.

```Shell
curl 'https://api.gdc.cancer.gov/cnvs/544c4896-0152-5787-8d77-894a16f0ded0?pretty=true&expand=consequence.gene'
```

```Json
{
  "data": {
    "start_position": 27100354,
    "consequence": [
      {
        "gene": {
          "biotype": "protein_coding",
          "symbol": "HOXA2",
          "gene_id": "ENSG00000105996"
        }
      }
    ],
    "gene_level_cn": true,
    "cnv_change": "Gain",
    "ncbi_build": "GRCh38",
    "chromosome": "7",
    "cnv_id": "544c4896-0152-5787-8d77-894a16f0ded0",
    "end_position": 27102686
  },
  "warnings": {}
}
```

## Copy Number Variation Occurrence Enpoint Examples

__Example 1:__ A user is interested in finding cases that have both cnv and ssm data for females diagnosed with Squamous Cell Neoplasms and have a cnv gain change on chromosome 9. It is important to note that for a case like this, where multiple arguments are need for one filtered field, it is easier for the API to have multiple filters for the same field, `case.available_variation_data` in this example, than having one filter with multiple arguments.

```Filter
{
    "op": "and",
    "content": [
        {
            "op": "in",
            "content": {
                "field": "cnv.cnv_change",
                "value": [
                    "Gain"
                ]
            }
        },
        {
            "op": "in",
            "content": {
                "field": "case.demographic.gender",
                "value": [
                    "female"
                ]
            }
        },
        {
            "op": "in",
            "content": {
                "field": "case.available_variation_data",
                "value": [
                    "cnv"
                ]
            }
        },
        {
            "op": "in",
            "content": {
                "field": "case.available_variation_data",
                "value": [
                    "ssm"
                ]
            }
        },
        {
            "op": "in",
            "content": {
                "field": "cnv.chromosome",
                "value": [
                    "9"
                ]
            }
        },
        {
            "op": "in",
            "content": {
                "field": "case.disease_type",
                "value": [
                    "Squamous Cell Neoplasms"
                ]
            }
        }
    ]
}

```

```Shell
curl 'https://api.gdc.cancer.gov/cnv_occurrences?filters=%7B%0D%0A++++%22op%22%3A+%22and%22%2C%0D%0A++++%22content%22%3A+%5B%0D%0A++++++++%7B%0D%0A++++++++++++%22op%22%3A+%22in%22%2C%0D%0A++++++++++++%22content%22%3A+%7B%0D%0A++++++++++++++++%22field%22%3A+%22cnv.cnv_change%22%2C%0D%0A++++++++++++++++%22value%22%3A+%5B%0D%0A++++++++++++++++++++%22Gain%22%0D%0A++++++++++++++++%5D%0D%0A++++++++++++%7D%0D%0A++++++++%7D%2C%0D%0A++++++++%7B%0D%0A++++++++++++%22op%22%3A+%22in%22%2C%0D%0A++++++++++++%22content%22%3A+%7B%0D%0A++++++++++++++++%22field%22%3A+%22case.demographic.gender%22%2C%0D%0A++++++++++++++++%22value%22%3A+%5B%0D%0A++++++++++++++++++++%22female%22%0D%0A++++++++++++++++%5D%0D%0A++++++++++++%7D%0D%0A++++++++%7D%2C%0D%0A++++++++%7B%0D%0A++++++++++++%22op%22%3A+%22in%22%2C%0D%0A++++++++++++%22content%22%3A+%7B%0D%0A++++++++++++++++%22field%22%3A+%22case.available_variation_data%22%2C%0D%0A++++++++++++++++%22value%22%3A+%5B%0D%0A++++++++++++++++++++%22cnv%22%0D%0A++++++++++++++++%5D%0D%0A++++++++++++%7D%0D%0A++++++++%7D%2C%0D%0A++++++++%7B%0D%0A++++++++++++%22op%22%3A+%22in%22%2C%0D%0A++++++++++++%22content%22%3A+%7B%0D%0A++++++++++++++++%22field%22%3A+%22case.available_variation_data%22%2C%0D%0A++++++++++++++++%22value%22%3A+%5B%0D%0A++++++++++++++++++++%22ssm%22%0D%0A++++++++++++++++%5D%0D%0A++++++++++++%7D%0D%0A++++++++%7D%2C%0D%0A++++++++%7B%0D%0A++++++++++++%22op%22%3A+%22in%22%2C%0D%0A++++++++++++%22content%22%3A+%7B%0D%0A++++++++++++++++%22field%22%3A+%22cnv.chromosome%22%2C%0D%0A++++++++++++++++%22value%22%3A+%5B%0D%0A++++++++++++++++++++%229%22%0D%0A++++++++++++++++%5D%0D%0A++++++++++++%7D%0D%0A++++++++%7D%2C%0D%0A++++++++%7B%0D%0A++++++++++++%22op%22%3A+%22in%22%2C%0D%0A++++++++++++%22content%22%3A+%7B%0D%0A++++++++++++++++%22field%22%3A+%22case.disease_type%22%2C%0D%0A++++++++++++++++%22value%22%3A+%5B%0D%0A++++++++++++++++++++%22Squamous+Cell+Neoplasms%22%0D%0A++++++++++++++++%5D%0D%0A++++++++++++%7D%0D%0A++++++++%7D%0D%0A++++%5D%0D%0A%7D&fields=case.available_variation_data,case.case_id&format=tsv'
```

```tsv
case.case_id	case.available_variation_data.1	case.available_variation_data.0	id
638035f6-2909-4a44-980f-468ac5d74e18	ssm	cnv	e76d2aaf-f951-5a51-a949-a241dba61f73
ad98977b-e159-410a-b8c2-f4e8a07f9784	ssm	cnv	ff3506b8-ee80-570f-ad2d-4ab4a7363b82
c83c52f4-3815-4f49-8218-cf80aaa62e2f	ssm	cnv	e73696c5-386f-5cae-aa10-f8628f32ee0e
dac27c24-cdbf-4527-9214-178fde3d098a	ssm	cnv	77885824-fae1-5116-9851-694255249cc8
0e91d7b5-ce35-4671-ab9f-cfd5369b557c	ssm	cnv	526529ae-8e59-597e-aea1-cc0b06a82e76
ea34663c-f40e-4a3e-9ac0-65d5e9eef12b	ssm	cnv	e4a0c034-44d4-5dea-912a-ce331d9a9512
05026179-b1da-411e-a286-89727b1ae380	ssm	cnv	30bdc04c-54a5-53ca-bdd0-b808f23da266
f1a1bbf9-4751-4fb4-8a2b-19f8d4ba57bd	ssm	cnv	02e3fbb3-da8f-5983-8d10-189e641ddf11
a6ec75d4-1c90-4527-bfae-aa91d2dae082	ssm	cnv	94b0e8be-1130-5b88-9103-6756bdabf67b
107f6b9a-2883-4499-a40a-ec25bc834a06	ssm	cnv	ad831f27-e6f5-5b78-8a15-0b652621ea4c
```

__Example 2:__ A user is interested in the first cnv occurrence (`e76d2aaf-f951-5a51-a949-a241dba61f73`) from the previous example, and would like to know more about the case exposures and demographics.

```Shell
curl 'https://api.gdc.cancer.gov/cnv_occurrences/e76d2aaf-f951-5a51-a949-a241dba61f73?pretty=true&expand=cnv,case,case.exposures,case.demographic'
```

```Json
{
  "data": {
    "cnv": {
      "ncbi_build": "GRCh38",
      "cnv_id": "0d475712-c11e-51fb-b6e6-407d12978057",
      "gene_level_cn": true,
      "cnv_change": "Gain",
      "end_position": 133348131,
      "variant_status": "Tumor only",
      "start_position": 133338323,
      "chromosome": "9"
    },
    "case": {
      "disease_type": "Squamous Cell Neoplasms",
      "updated_datetime": "2018-09-06T11:07:45.510627-05:00",
      "created_datetime": null,
      "demographic": {
        "updated_datetime": "2018-09-06T11:07:45.510627-05:00",
        "created_datetime": null,
        "gender": "female",
        "year_of_birth": 1954,
        "submitter_id": "TCGA-EA-A3HR_demographic",
        "state": "released",
        "race": "white",
        "demographic_id": "dd8576a8-bd62-55e7-b0df-7233ceded2fb",
        "ethnicity": "not hispanic or latino",
        "year_of_death": null
      },
      "submitter_id": "TCGA-EA-A3HR",
      "state": "released",
      "case_id": "638035f6-2909-4a44-980f-468ac5d74e18",
      "primary_site": "Cervix uteri",
      "available_variation_data": [
        "cnv",
        "ssm"
      ],
      "exposures": [
        {
          "cigarettes_per_day": null,
          "weight": 86,
          "updated_datetime": "2018-09-06T11:07:45.510627-05:00",
          "created_datetime": null,
          "alcohol_intensity": null,
          "bmi": 40,
          "years_smoked": null,
          "submitter_id": "TCGA-EA-A3HR_exposure",
          "alcohol_history": null,
          "state": "released",
          "tobacco_smoking_status": null,
          "tobacco_smoking_onset_year": null,
          "tobacco_smoking_quit_year": null,
          "exposure_id": "0e7265ab-bf65-50c7-bf33-96a7ac452d7c",
          "height": 146,
          "pack_years_smoked": null
        }
      ]
    },
    "cnv_occurrence_id": "e76d2aaf-f951-5a51-a949-a241dba61f73"
  }
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

[Survival plots](/Data_Portal/Users_Guide/Exploration/#survival-analysis) are generated for different subsets of data, based on variants or projects, in the GDC Data Portal. The `/analysis/survival` endpoint can be used to programmatically retrieve the raw data used to generate these plots and apply different filters. Note that the `fields` and `format` parameters cannot be modified.

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
          "value":"edd1ae2c-3ca9-52bd-a124-b09ed304fcc2"
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
          "value":"edd1ae2c-3ca9-52bd-a124-b09ed304fcc2"
        }
      }
    ]
  }
]
```
```Shell
curl "https://api.gdc.cancer.gov/analysis/survival?filters=%5B%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%22TCGA-BRCA%22%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22gene.ssm.ssm_id%22%2C%22value%22%3A%22edd1ae2c-3ca9-52bd-a124-b09ed304fcc2%22%7D%7D%5D%7D%2C%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%22TCGA-BRCA%22%7D%7D%2C%7B%22op%22%3A%22excludeifany%22%2C%22content%22%3A%7B%22field%22%3A%22gene.ssm.ssm_id%22%2C%22value%22%3A%22edd1ae2c-3ca9-52bd-a124-b09ed304fcc2%22%7D%7D%5D%7D%5D&pretty=true"
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

__Example 3:__ Custom survival plots can be generated using the GDC API.  For example, a user could generate survival plot data comparing patients with a mutation in genes associated with a biological pathway with patients without mutations in that pathway. The following example compares a patient with at least one mutation in either gene `ENSG00000141510` or `ENSG00000155657` with patients that do not have mutations in these genes.

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

__Example 4:__ Survival plots can be even more customizable when sets of case IDs are used. Two sets of case IDs (or barcodes) can be retrieved in a separate step based on custom criteria and compared in a survival plot. See below for an example query.

```Json
[{  
   "op":"=",
   "content":{  
      "field":"cases.submitter_id",
      "value":["TCGA-HT-A74J","TCGA-43-A56U","TCGA-GM-A3XL","TCGA-A1-A0SQ","TCGA-K1-A6RV","TCGA-J2-A4AD","TCGA-XR-A8TE"]
   }
},
{  
   "op":"=",
   "content":{  
      "field":"cases.submitter_id",
      "value":["TCGA-55-5899","TCGA-55-6642","TCGA-55-7907","TCGA-67-6216","TCGA-75-5146","TCGA-49-4510","TCGA-78-7159"]
   }
}]
```
```Shell
curl "https://api.gdc.cancer.gov/analysis/survival?filters=%5B%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.submitter_id%22%2C%22value%22%3A%5B%22TCGA-HT-A74J%22%2C%22TCGA-43-A56U%22%2C%22TCGA-GM-A3XL%22%2C%22TCGA-A1-A0SQ%22%2C%22TCGA-K1-A6RV%22%2C%22TCGA-J2-A4AD%22%2C%22TCGA-XR-A8TE%22%5D%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.submitter_id%22%2C%22value%22%3A%5B%22TCGA-55-5899%22%2C%22TCGA-55-6642%22%2C%22TCGA-55-7907%22%2C%22TCGA-67-6216%22%2C%22TCGA-75-5146%22%2C%22TCGA-49-4510%22%2C%22TCGA-78-7159%22%5D%7D%7D%5D"
```
