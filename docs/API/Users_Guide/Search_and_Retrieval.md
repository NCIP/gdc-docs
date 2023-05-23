# Search and Retrieval

## Introducing Search and Retrieval Requests

The GDC API provides endpoints that search and retrieve information stored in the GDC according to the [GDC Data Model](../../Data/Data_Model/GDC_Data_Model.md). The general format of requests to search & retrieval endpoints is described below.

>**Note:** Queries described in this section work for datasets that have been released to the GDC Data Portal. Unreleased data that is in the process of being submitted to GDC cannot be queried using these methods. See [Submission](Submission.md) to learn how to query unreleased data using GraphQL.

### Components of a Request

A typical search and retrieval API request specifies the following parameters:

- a `filters` parameter, that specifies the search terms for the query
- several parameters that specify the API response, such as:
	- `format` &mdash; specifies response format (JSON, TSV, XML)
	- `fields` &mdash; specifies the which data elements should be returned in the response, if available
	- `size` &mdash; specifies the the maximum number of results to include in the response
	- other parameters are described below.

Requests can be executed using HTTP GET or HTTP POST. GET requests are limited by maximum URL length, so the POST method is recommended for large queries.

### POST Example

The following is an example of an HTTP POST request to the `files` endpoint of the GDC API. It looks for Gene Expression Quantification files associated with specific TCGA cases (represented by TCGA barcodes) and retrieves the associated biospecimen metadata in TSV format.

#### Request

	curl --request POST --header "Content-Type: application/json" --data @Payload 'https://api.gdc.cancer.gov/files' > response.tsv

#### Payload

	{
	    "filters":{
	        "op":"and",
	        "content":[
	            {
	                "op":"in",
	                "content":{
	                    "field":"cases.submitter_id",
	                    "value":[
	                        "TCGA-CK-4948",
	                        "TCGA-D1-A17N",
	                        "TCGA-4V-A9QX",
	                        "TCGA-4V-A9QM"
	                    ]
	                }
	            },
	            {
	                "op":"=",
	                "content":{
	                    "field":"files.data_type",
	                    "value":"Gene Expression Quantification"
	                }
	            }
	        ]
	    },
	    "format":"tsv",
	    "fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,analysis.workflow_type,cases.project.project_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id",
	    "size":"1000"
	}

Each component of the request is explained below.

### GET Example

The above request can be executed as an HTTP GET:

	https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.submitter_id%22%2C%22value%22%3A%5B%22TCGA-CK-4948%22%2C%22TCGA-D1-A17N%22%2C%22TCGA-4V-A9QX%22%2C%22TCGA-4V-A9QM%22%5D%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%22Gene%20Expression%20Quantification%22%7D%7D%5D%7D&format=tsv&fields=file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,analysis.workflow_type,cases.project.project_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id&size=1000

Each component of the request is explained below.


## Endpoints

The following search and retrieval endpoints are available in the GDC API:

| Endpoints | Description |
| --- | --- |
| [files](/API/Users_Guide/Search_and_Retrieval/#files-endpoint) | Information about files stored in the GDC |
| [cases](/API/Users_Guide/Search_and_Retrieval/#cases-endpoint) | Information related to cases, or sample donors |
| [history](/API/Users_Guide/Search_and_Retrieval/#history-endpoint) | Information related to file version history |
| [projects](/API/Users_Guide/Search_and_Retrieval/#project-endpoint) | Information about projects |
| [annotations](/API/Users_Guide/Search_and_Retrieval/#annotations-endpoint) | Information about annotations to GDC data |
| \_mapping | Information about elements that can be used to query other endpoints |

The choice of endpoint determines what is listed in the search results. The `files` endpoint will generate a list of files, whereas the `cases` endpoint will generate a list of cases. Each of the above endpoints, other than `_mapping`, can query and return any of the related fields in the [GDC Data Model](../../Data/Data_Model/GDC_Data_Model.md). So the `cases` endpoint can be queried for file fields (e.g. to look for cases that have certain types of experimental data), and the `files` endpoint can be queried for clinical metadata associated with a case (e.g. to look for files from cases diagnosed with a specific cancer type).

### `Project` Endpoint
The `projects` endpoint provides access to project records, the highest level of data organization in the GDC.

#### Example
This example is a query for projects contained in the GDC. It uses the [from](#from), [size](#size), [sort](#sort), and [pretty](#pretty) parameters, and returns the first two projects sorted by project id.

```shell
curl 'https://api.gdc.cancer.gov/projects?from=0&size=2&sort=project_id:asc&pretty=true'
```
``` Output
{
  "data": {
    "hits": [
      {
        "id": "APOLLO-LUAD",
        "primary_site": [
          "Bronchus and lung"
        ],
        "dbgap_accession_number": "phs003011",
        "project_id": "APOLLO-LUAD",
        "disease_type": [
          "Adenomas and Adenocarcinomas"
        ],
        "name": "APOLLO1: Proteogenomic characterization of lung adenocarcinoma",
        "releasable": false,
        "state": "open",
        "released": true
      },
      {
        "id": "BEATAML1.0-COHORT",
        "primary_site": [
          "Hematopoietic and reticuloendothelial systems"
        ],
        "dbgap_accession_number": "phs001657",
        "project_id": "BEATAML1.0-COHORT",
        "disease_type": [
          "Myelodysplastic Syndromes",
          "Leukemias, NOS",
          "Unknown",
          "Myeloid Leukemias",
          "Plasma Cell Tumors",
          "Chronic Myeloproliferative Disorders"
        ],
        "name": "Functional Genomic Landscape of Acute Myeloid Leukemia",
        "releasable": true,
        "state": "open",
        "released": true
      }
    ],
    "pagination": {
      "count": 2,
      "total": 78,
      "size": 2,
      "from": 0,
      "sort": "None",
      "page": 1,
      "pages": 39
    }
  },
  "warnings": {}
}
```

#### Retrieval of project metadata using project_id

The `project` endpoint supports a simple query format that retrieves the metadata of a single project using its `project_id`:

```shell
curl 'https://api.gdc.cancer.gov/projects/TARGET-NBL?expand=summary,summary.experimental_strategies,summary.data_categories&pretty=true'
```
```Response
{
  "data": {
    "summary": {
      "file_count": 5705,
      "data_categories": [
        {
          "file_count": 943,
          "case_count": 278,
          "data_category": "Sequencing Reads"
        },
        {
          "file_count": 3080,
          "case_count": 220,
          "data_category": "Simple Nucleotide Variation"
        },
        {
          "file_count": 3,
          "case_count": 1119,
          "data_category": "Clinical"
        },
        {
          "file_count": 705,
          "case_count": 225,
          "data_category": "DNA Methylation"
        },
        {
          "file_count": 2,
          "case_count": 1132,
          "data_category": "Biospecimen"
        },
        {
          "file_count": 324,
          "case_count": 155,
          "data_category": "Transcriptome Profiling"
        },
        {
          "file_count": 648,
          "case_count": 155,
          "data_category": "Structural Variation"
        }
      ],
      "experimental_strategies": [
        {
          "file_count": 1458,
          "case_count": 155,
          "experimental_strategy": "RNA-Seq"
        },
        {
          "file_count": 15,
          "case_count": 8,
          "experimental_strategy": "WGS"
        },
        {
          "file_count": 3522,
          "case_count": 222,
          "experimental_strategy": "WXS"
        },
        {
          "file_count": 705,
          "case_count": 225,
          "experimental_strategy": "Methylation Array"
        }
      ],
      "case_count": 1132,
      "file_size": 16968781125824
    },
    "primary_site": [
      "Stomach",
      "Bones, joints and articular cartilage of limbs",
      "Heart, mediastinum, and pleura",
      "Peripheral nerves and autonomic nervous system",
      "Uterus, NOS",
      "Bones, joints and articular cartilage of other and unspecified sites",
      "Other endocrine glands and related structures",
      "Renal pelvis",
      "Retroperitoneum and peritoneum",
      "Liver and intrahepatic bile ducts",
      "Meninges",
      "Connective, subcutaneous and other soft tissues",
      "Adrenal gland",
      "Unknown",
      "Spinal cord, cranial nerves, and other parts of central nervous system",
      "Skin",
      "Other and ill-defined sites",
      "Kidney",
      "Lymph nodes",
      "Hematopoietic and reticuloendothelial systems"
    ],
    "dbgap_accession_number": "phs000467",
    "project_id": "TARGET-NBL",
    "disease_type": [
      "Neuroepitheliomatous Neoplasms",
      "Not Applicable"
    ],
    "name": "Neuroblastoma",
    "releasable": true,
    "state": "open",
    "released": true
  },
  "warnings": {}
}

```

### `Files` Endpoint

The GDC Files Endpoint `https://api.gdc.cancer.gov/files` enables search and retrieval of information relating to files stored in the GDC, including file properties such as `file_name`, `md5sum`, `data_format`, and others.

#### Example

This example is a query for files contained in the GDC. It uses the [from](#from), [size](#size), [sort](#sort), and [pretty](#pretty) parameters, and returns only the first two files, sorted by file size, from smallest to largest.

```shell
curl 'https://api.gdc.cancer.gov/files?from=0&size=2&sort=file_size:asc&pretty=true'
```
``` Output
{
  "data": {
    "hits": [
      {
        "id": "0ab5e358-b1ff-4433-8959-c37c5890d9aa",
        "data_format": "BEDPE",
        "access": "controlled",
        "file_name": "090e2828-079c-48e6-97cb-735c763da8d3.wgs.BRASS.rerun_structural_variation.bedpe.gz",
        "submitter_id": "247c3c9a-58b9-4b70-bda8-cb197acb5609",
        "data_category": "Somatic Structural Variation",
        "acl": [
          "phs001287"
        ],
        "type": "structural_variation",
        "file_size": 20,
        "created_datetime": "2022-04-08T20:27:04.633842-05:00",
        "md5sum": "7029066c27ac6f5ef18d660d5741979a",
        "updated_datetime": "2022-07-07T11:02:27.204310-05:00",
        "file_id": "0ab5e358-b1ff-4433-8959-c37c5890d9aa",
        "data_type": "Structural Rearrangement",
        "state": "released",
        "experimental_strategy": "WGS",
        "version": "1",
        "data_release": "34.0 - 37.0"
      },
      {
        "id": "a8bc2405-b57d-48bb-b241-18b3e28caa56",
        "data_format": "BEDPE",
        "access": "controlled",
        "file_name": "eae76f14-8aa7-427f-a90c-4e0ed095e0c2.wgs.BRASS.rerun_structural_variation.bedpe.gz",
        "submitter_id": "618cd251-ddcb-4a7e-9a6d-efb132b0bd7a",
        "data_category": "Somatic Structural Variation",
        "acl": [
          "phs001287"
        ],
        "type": "structural_variation",
        "file_size": 20,
        "created_datetime": "2022-04-08T20:43:16.505747-05:00",
        "md5sum": "7029066c27ac6f5ef18d660d5741979a",
        "updated_datetime": "2022-07-07T11:00:43.345766-05:00",
        "file_id": "a8bc2405-b57d-48bb-b241-18b3e28caa56",
        "data_type": "Structural Rearrangement",
        "state": "released",
        "experimental_strategy": "WGS",
        "version": "1",
        "data_release": "34.0 - 37.0"
      }
    ],
    "pagination": {
      "count": 2,
      "total": 931947,
      "size": 2,
      "from": 0,
      "sort": "None",
      "page": 1,
      "pages": 465974
    }
  },
  "warnings": {}
}

```

#### Retrieval of file metadata using individual UUIDs:

The `/files` endpoint supports a simple query format that retrieves the metadata of a single file using its UUID.  Note that the `/files` endpoint is inactive when querying for earlier file versions.  In that case, the `/history` or `/files/versions` endpoints should be used instead.

```Shell
curl 'https://api.gdc.cancer.gov/files/20f45e04-3c10-4f11-b57b-719880eab69e?pretty=true'
```
``` Output
{
  "data": {
    "data_format": "VCF",
    "access": "controlled",
    "file_name": "TCGA_BRCA.8d9cb5ae-e568-41fc-8b53-14467c2623dc.wxs.MuTect2.somatic_annotation.vcf.gz",
    "submitter_id": "675f31dd-70e5-4a72-8139-423b14b31564",
    "data_category": "Simple Nucleotide Variation",
    "acl": [
      "phs000178"
    ],
    "type": "annotated_somatic_mutation",
    "file_size": 6894331,
    "created_datetime": "2022-02-07T08:48:39.178606-06:00",
    "md5sum": "f4c0f52ba5bf24f0ca6ba3a923fecc5e",
    "updated_datetime": "2022-02-09T12:11:12.781445-06:00",
    "file_id": "20f45e04-3c10-4f11-b57b-719880eab69e",
    "data_type": "Annotated Somatic Mutation",
    "state": "released",
    "experimental_strategy": "WXS",
    "version": "2",
    "data_release": "32.0 - 37.0"
  },
  "warnings": {}
}
```

>__Note:__ The `file_size` field associated with each file is reported in bytes.  


#### Example of retrieving file version information:

The `https://api.gdc.cancer.gov/files/versions` endpoint enables search and retrieval of version information about a file.  A file may be versioned if a file is updated by the GDC (e.g. using a new alignment algorithm or fixing a file that contained an error). `Version` refers to the instance of a particular file.  Inputs can either be a list of UUIDs as shown in example 1 or a download manifest as shown in example 2.  Output includes information about the current and latest version for any given file.  While `/files` also returns information about a file version this endpoint will only work for the most recent version of a file whereas `/files/versions` will work for all previous and current versions of a file.  In both examples below the output format can be modified by adding the `format=tsv` parameter.

```Shell1
curl 'https://api.gdc.cancer.gov/files/versions/1dd28069-5777-4ff9-bd2b-d1ba68e88b06,2a03abac-f1a2-49a9-a57c-7543739dd862?pretty=true'
```

``` Output1
[
  {
    "id": "1dd28069-5777-4ff9-bd2b-d1ba68e88b06",
    "filename": "1dd28069-5777-4ff9-bd2b-d1ba68e88b06.vcf.gz",
    "version": "1",
    "md5": "c2f9b196e154906a70c7ec46492a859d",
    "size": 332092,
    "state": "validated",
    "release": "12.0",
    "latest_id": "76b3f4d8-c6b7-4662-ac42-1d27d4684281",
    "latest_filename": "def1cc5b-55f0-4372-a3ff-df3ea93cf3e7.wxs.somaticsniper.raw_somatic_mutation.vcf.gz",
    "latest_version": "2",
    "latest_md5": "bbc66201eeb12e8f63fc6dcc156dbac9",
    "latest_size": 357706,
    "latest_state": "validated",
    "latest_release": [
      "32.0",
      "33.0",
      "33.1",
      "34.0",
      "35.0",
      "36.0",
      "37.0"
    ]
  },
  {
    "id": "2a03abac-f1a2-49a9-a57c-7543739dd862",
    "filename": "a5d86cde-32ca-4ed6-b1a5-5a47575f2ac6_gdc_realn_rehead.bam",
    "version": "1",
    "md5": "48686fcd84ac713d44261ca9e26b89fb",
    "size": 6653119038,
    "state": "validated",
    "release": "12.0",
    "latest_id": "de0ce84d-c286-405c-a556-39dac14c7c74",
    "latest_filename": "d45c33cc-88e2-4de5-a578-f7e31a6c0738.rna_seq.genomic.gdc_realn.bam",
    "latest_version": "2",
    "latest_md5": "5d4b0c13f1c1a235ed064d94278bc196",
    "latest_size": 6223445806,
    "latest_state": "validated",
    "latest_release": [
      "32.0",
      "33.0",
      "33.1",
      "34.0",
      "35.0",
      "36.0",
      "37.0"
    ]
  }
]
```
```Shell2
curl --request POST --header "Content-Type: text/tsv"  https://api.gdc.cancer.gov/files/versions/manifest?pretty=true --data-binary @gdc_manifest_20180809_154816.txt
```

``` Output2
[
  {
    "id": "0b20e27c-9a09-4f15-923f-d5b4f185dc22",
    "filename": "nationwidechildrens.org_clinical.TCGA-13-1500.xml",
    "version": "1",
    "md5": "597aa4df24c4d544b6c25cbd8b25a33e",
    "size": 44857,
    "state": "validated",
    "release": "12.0",
    "latest_id": "0b20e27c-9a09-4f15-923f-d5b4f185dc22",
    "latest_filename": "nationwidechildrens.org_clinical.TCGA-13-1500.xml",
    "latest_version": "1",
    "latest_md5": "597aa4df24c4d544b6c25cbd8b25a33e",
    "latest_size": 44857,
    "latest_state": "validated",
    "latest_release": [
      "12.0",
      "13.0",
      "14.0",
      "15.0",
      "16.0",
      "17.0",
      "18.0",
      "19.0",
      "20.0",
      "21.0",
      "22.0",
      "23.0",
      "24.0",
      "25.0",
      "26.0",
      "27.0",
      "28.0",
      "29.0",
      "30.0",
      "31.0",
      "32.0",
      "33.0",
      "33.1",
      "34.0",
      "35.0",
      "36.0",
      "37.0"
    ]
  },
  {
    "id": "3edc7084-013c-4493-8507-c00b0e9962d8",
    "filename": "BUCKS_p_TCGA_272_273_N_GenomeWideSNP_6_G05_1320676.grch38.seg.v2.txt",
    "version": "1",
    "md5": "35a18d990a05eedfaf96e753bee0b96d",
    "size": 27620,
    "state": "validated",
    "release": "12.0",
    "latest_id": "3edc7084-013c-4493-8507-c00b0e9962d8",
    "latest_filename": "BUCKS_p_TCGA_272_273_N_GenomeWideSNP_6_G05_1320676.grch38.seg.v2.txt",
    "latest_version": "1",
    "latest_md5": "35a18d990a05eedfaf96e753bee0b96d",
    "latest_size": 27620,
    "latest_state": "validated",
    "latest_release": [
      "12.0",
      "13.0",
      "14.0",
      "15.0",
      "16.0",
      "17.0",
      "18.0",
      "19.0",
      "20.0",
      "21.0",
      "22.0",
      "23.0",
      "24.0",
      "25.0",
      "26.0",
      "27.0",
      "28.0",
      "29.0",
      "30.0",
      "31.0",
      "32.0",
      "33.0",
      "33.1",
      "34.0",
      "35.0",
      "36.0",
      "37.0"
    ]
  },
  {
    "id": "a22f5e32-b16e-458f-a412-7e438056ece6",
    "filename": "a22f5e32-b16e-458f-a412-7e438056ece6.vep.vcf.gz",
    "version": "1",
    "md5": "68b2433b31679bbbc6681919a1b81762",
    "size": 2346,
    "state": "validated",
    "release": "12.0",
    "latest_id": "55491171-6170-45cb-af9d-d99345b289e5",
    "latest_filename": "4b89bb97-41f6-43c4-a481-287556f7bb4a.targeted_sequencing.annotated_somatic_mutation.vcf.gz",
    "latest_version": "2",
    "latest_md5": "5b2faa8780317744331b8b852f7d9ef1",
    "latest_size": 2618,
    "latest_state": "validated",
    "latest_release": [
      "32.0",
      "33.0",
      "33.1",
      "34.0",
      "35.0",
      "36.0",
      "37.0"
    ]
  }
]
```

### `Cases` Endpoint

The GDC Cases Endpoint `https://api.gdc.cancer.gov/cases` enables search and retrieval of information related to a specific case.

The `cases` endpoint is designed to retrieve the metadata associated with one or more cases, including all nested biospecimen entities. Filters can be applied to retrieve information for entire cases, but not for lower-level biospecimen entities. For example, a sample within a case cannot be used to query for aliquots that are associated only with that sample. All aliquots associated with the case would be retrieved.


#### Example

This example is a query for files contained in GDC. It returns case where submitter id is `TCGA-BH-A0EA`, using the [pretty](#pretty) and [filters](#filters) parameters and the following [filtering operators](#filtering-operators):

	{"op":"and","content":[{"op":"in","content":{"field":"submitter_id","value":["TCGA-BH-A0EA"]}}]}

Command:

```shell
curl 'https://api.gdc.cancer.gov/cases?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22submitter_id%22%2C%22value%22%3A%5B%22TCGA-BH-A0EA%22%5D%7D%7D%5D%7D%0A%0A&pretty=true'
```
``` Output
{
  "data": {
    "hits": [
      {
        "id": "1f601832-eee3-48fb-acf5-80c4a454f26e",
        "slide_ids": [
          "a0826f0d-986a-491b-8c6f-b34f8929f3ee",
          "90154ea1-6b76-4445-870e-d531d6fa1239",
          "1dd1cab5-5a81-428a-8153-91e8c4cf9905"
        ],
        "submitter_slide_ids": [
          "TCGA-BH-A0EA-01Z-00-DX1",
          "TCGA-BH-A0EA-01A-01-MSA",
          "TCGA-BH-A0EA-01A-01-TSA"
        ],
        "disease_type": "Ductal and Lobular Neoplasms",
        "analyte_ids": [
          "f19f408a-815f-43d9-8032-e9482b796371",
          "fe678556-acf4-4bde-a95e-860bb0150a95",
          "69ddc092-88a0-4839-a2bb-9f1c9e760409",
          "66ed0f86-5ca5-4dec-ba76-7ee4dcf31831",
          "30cb470f-66d4-4085-8c30-83a42e8453d4"
        ],
        "submitter_id": "TCGA-BH-A0EA",
        "submitter_analyte_ids": [
          "TCGA-BH-A0EA-10A-01D",
          "TCGA-BH-A0EA-01A-11D",
          "TCGA-BH-A0EA-01A-11R",
          "TCGA-BH-A0EA-10A-01W",
          "TCGA-BH-A0EA-01A-11W"
        ],
        "aliquot_ids": [
          "cde982b7-3b0a-49eb-8710-a599cb0e44c1",
          "b1a3739d-d554-4202-b96f-f25a444e2042",
          "eef9dce1-6ba6-432b-bbe2-53c7dbe64fe7",
          "97c64d6a-7dce-4d0f-9cb3-b3e4eb4719c5",
          "262715e1-835c-4f16-8ee7-6900e26f7cf5",
          "cfbd5476-e83a-401d-9f9a-639c73a0e35b",
          "bcb7fc6d-60a0-48b7-aa81-14c0dda72d76",
          "561b8777-801a-49ed-a306-e7dafeb044b6",
          "edad5bd3-efe0-4c5f-b05c-2c0c2951c45a",
          "42d050e4-e8ee-4442-b9c0-0ee14706b138",
          "2beb34c4-d493-4a73-b21e-de77d43251ff",
          "ca71ca96-cbb7-4eab-9487-251dda34e107"
        ],
        "submitter_aliquot_ids": [
          "TCGA-BH-A0EA-10A-01W-A12U-09",
          "TCGA-BH-A0EA-01A-11D-A111-01",
          "TCGA-BH-A0EA-01A-11D-A314-09",
          "TCGA-BH-A0EA-01A-11D-A10X-02",
          "TCGA-BH-A0EA-10A-01D-A10Z-02",
          "TCGA-BH-A0EA-10A-01D-A110-09",
          "TCGA-BH-A0EA-01A-11D-A10Y-09",
          "TCGA-BH-A0EA-10A-01D-A113-01",
          "TCGA-BH-A0EA-01A-11D-A112-05",
          "TCGA-BH-A0EA-01A-11R-A115-07",
          "TCGA-BH-A0EA-01A-11W-A12T-09",
          "TCGA-BH-A0EA-01A-11R-A114-13"
        ],
        "created_datetime": null,
        "diagnosis_ids": [
          "84654ad5-2a2c-5c3b-8340-ecac6a5550fe"
        ],
        "sample_ids": [
          "55864d86-dab8-47bb-a3e3-8cfb198b06c1",
          "9a6c71a6-82cd-42b1-a93f-f569370848d6",
          "7f791228-dd77-4ab0-8227-d784a4c7fea1"
        ],
        "submitter_sample_ids": [
          "TCGA-BH-A0EA-01A",
          "TCGA-BH-A0EA-01Z",
          "TCGA-BH-A0EA-10A"
        ],
        "primary_site": "Breast",
        "submitter_diagnosis_ids": [
          "TCGA-BH-A0EA_diagnosis"
        ],
        "updated_datetime": "2019-08-06T14:15:54.128069-05:00",
        "case_id": "1f601832-eee3-48fb-acf5-80c4a454f26e",
        "state": "released",
        "portion_ids": [
          "cb6086d1-3416-4310-b109-e8fa6e8b72d4",
          "8629bf5a-cdaf-4f6a-90bb-27dd4a7565c5",
          "ae4f5816-f97a-4605-9b05-9ab820467dee"
        ],
        "submitter_portion_ids": [
          "TCGA-BH-A0EA-10A-01",
          "TCGA-BH-A0EA-01A-21-A13C-20",
          "TCGA-BH-A0EA-01A-11"
        ]
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

#### Retrieval of case metadata using individual UUIDs:

The `cases` endpoint supports a simple query format that retrieves the metadata of a single case using its UUID:

```shell
curl 'https://api.gdc.cancer.gov/cases/1f601832-eee3-48fb-acf5-80c4a454f26e?pretty=true&expand=diagnoses'
```
```Response

  "data": {
    "slide_ids": [
      "90154ea1-6b76-4445-870e-d531d6fa1239",
      "1dd1cab5-5a81-428a-8153-91e8c4cf9905",
      "a0826f0d-986a-491b-8c6f-b34f8929f3ee"
    ],
    "submitter_slide_ids": [
      "TCGA-BH-A0EA-01A-01-MSA",
      "TCGA-BH-A0EA-01A-01-TSA",
      "TCGA-BH-A0EA-01Z-00-DX1"
    ],
    "disease_type": "Ductal and Lobular Neoplasms",
    "analyte_ids": [
      "fe678556-acf4-4bde-a95e-860bb0150a95",
      "66ed0f86-5ca5-4dec-ba76-7ee4dcf31831",
      "30cb470f-66d4-4085-8c30-83a42e8453d4",
      "69ddc092-88a0-4839-a2bb-9f1c9e760409",
      "f19f408a-815f-43d9-8032-e9482b796371"
    ],
    "submitter_id": "TCGA-BH-A0EA",
    "submitter_analyte_ids": [
      "TCGA-BH-A0EA-01A-11D",
      "TCGA-BH-A0EA-01A-11R",
      "TCGA-BH-A0EA-10A-01W",
      "TCGA-BH-A0EA-01A-11W",
      "TCGA-BH-A0EA-10A-01D"
    ],
    "aliquot_ids": [
      "eef9dce1-6ba6-432b-bbe2-53c7dbe64fe7",
      "2beb34c4-d493-4a73-b21e-de77d43251ff",
      "b1a3739d-d554-4202-b96f-f25a444e2042",
      "262715e1-835c-4f16-8ee7-6900e26f7cf5",
      "cfbd5476-e83a-401d-9f9a-639c73a0e35b",
      "edad5bd3-efe0-4c5f-b05c-2c0c2951c45a",
      "bcb7fc6d-60a0-48b7-aa81-14c0dda72d76",
      "42d050e4-e8ee-4442-b9c0-0ee14706b138",
      "97c64d6a-7dce-4d0f-9cb3-b3e4eb4719c5",
      "561b8777-801a-49ed-a306-e7dafeb044b6",
      "ca71ca96-cbb7-4eab-9487-251dda34e107",
      "cde982b7-3b0a-49eb-8710-a599cb0e44c1"
    ],
    "submitter_aliquot_ids": [
      "TCGA-BH-A0EA-01A-11R-A115-07",
      "TCGA-BH-A0EA-01A-11D-A112-05",
      "TCGA-BH-A0EA-10A-01W-A12U-09",
      "TCGA-BH-A0EA-01A-11D-A10X-02",
      "TCGA-BH-A0EA-10A-01D-A113-01",
      "TCGA-BH-A0EA-10A-01D-A110-09",
      "TCGA-BH-A0EA-01A-11D-A314-09",
      "TCGA-BH-A0EA-01A-11D-A10Y-09",
      "TCGA-BH-A0EA-01A-11D-A111-01",
      "TCGA-BH-A0EA-10A-01D-A10Z-02",
      "TCGA-BH-A0EA-01A-11R-A114-13",
      "TCGA-BH-A0EA-01A-11W-A12T-09"
    ],
    "diagnoses": [
      {
        "synchronous_malignancy": "Not Reported",
        "ajcc_pathologic_stage": "Stage IIA",
        "days_to_diagnosis": 0,
        "created_datetime": null,
        "last_known_disease_status": "not reported",
        "tissue_or_organ_of_origin": "Breast, NOS",
        "days_to_last_follow_up": null,
        "age_at_diagnosis": 26548,
        "primary_diagnosis": "Infiltrating duct carcinoma, NOS",
        "updated_datetime": "2019-08-08T16:25:42.215495-05:00",
        "prior_malignancy": "yes",
        "year_of_diagnosis": 2008,
        "state": "released",
        "prior_treatment": "No",
        "days_to_last_known_disease_status": null,
        "ajcc_staging_system_edition": "6th",
        "ajcc_pathologic_t": "T1c",
        "days_to_recurrence": null,
        "morphology": "8500/3",
        "ajcc_pathologic_n": "N1a",
        "ajcc_pathologic_m": "M0",
        "submitter_id": "TCGA-BH-A0EA_diagnosis",
        "classification_of_tumor": "not reported",
        "diagnosis_id": "84654ad5-2a2c-5c3b-8340-ecac6a5550fe",
        "icd_10_code": "C50.9",
        "site_of_resection_or_biopsy": "Breast, NOS",
        "tumor_grade": "not reported",
        "progression_or_recurrence": "not reported"
      }
    ],
    "created_datetime": null,
    "diagnosis_ids": [
      "84654ad5-2a2c-5c3b-8340-ecac6a5550fe"
    ],
    "sample_ids": [
      "55864d86-dab8-47bb-a3e3-8cfb198b06c1",
      "7f791228-dd77-4ab0-8227-d784a4c7fea1",
      "9a6c71a6-82cd-42b1-a93f-f569370848d6"
    ],
    "submitter_sample_ids": [
      "TCGA-BH-A0EA-01A",
      "TCGA-BH-A0EA-01Z",
      "TCGA-BH-A0EA-10A"
    ],
    "primary_site": "Breast",
    "submitter_diagnosis_ids": [
      "TCGA-BH-A0EA_diagnosis"
    ],
    "updated_datetime": "2019-08-06T14:15:54.128069-05:00",
    "case_id": "1f601832-eee3-48fb-acf5-80c4a454f26e",
    "state": "released",
    "portion_ids": [
      "1ef8b20e-43e5-49d7-ac9a-03ce14f58daa",
      "cb6086d1-3416-4310-b109-e8fa6e8b72d4",
      "8629bf5a-cdaf-4f6a-90bb-27dd4a7565c5",
      "ae4f5816-f97a-4605-9b05-9ab820467dee"
    ],
    "submitter_portion_ids": [
      "TCGA-BH-A0EA-01A-21-A13C-20",
      "TCGA-BH-A0EA-10A-01",
      "TCGA-BH-A0EA-01A-21",
      "TCGA-BH-A0EA-01A-11"
    ]
  },
  "warnings": {}
}
```

### `Annotations` Endpoint

The GDC Annotation Endpoint `https://api.gdc.cancer.gov/annotations` enables search and retrieval of annotations stored in the GDC.


#### Example

This example is a query for any annotations **directly** associated with the following GDC entities:

* the case with UUID e0d36cc0-652c-4224-bb10-09d15c7bd8f1
* the sample with UUID 25ebc29a-7598-4ae4-ba7f-618d448882cc
* the aliquot with UUID fe660d7c-2746-4b50-ab93-b2ed99960553

The query uses the [filters](#filters) parameter to specify entity UUIDs. Code samples below include the bare and percent-encoded filter JSON.

```Filter-JSON
{
   "op":"in",
   "content":{
      "field":"entity_id",
      "value":[
         "e0d36cc0-652c-4224-bb10-09d15c7bd8f1",
         "25ebc29a-7598-4ae4-ba7f-618d448882cc",
         "fe660d7c-2746-4b50-ab93-b2ed99960553"
      ]
   }
}
```
```Filter-JSON-percent-encoded
%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22entity_id%22%2C%22value%22%3A%5B%22e0d36cc0-652c-4224-bb10-09d15c7bd8f1%22%2C%2225ebc29a-7598-4ae4-ba7f-618d448882cc%22%2C%22fe660d7c-2746-4b50-ab93-b2ed99960553%22%5D%7D%7D
```
```shell
curl 'https://api.gdc.cancer.gov/annotations?filters=%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22entity_id%22%2C%22value%22%3A%5B%22e0d36cc0-652c-4224-bb10-09d15c7bd8f1%22%2C%2225ebc29a-7598-4ae4-ba7f-618d448882cc%22%2C%22fe660d7c-2746-4b50-ab93-b2ed99960553%22%5D%7D%7D&pretty=true'
```
``` Output
{
  "data": {
    "hits": [
      {
        "category": "Item flagged DNU",
        "status": "Approved",
        "entity_id": "fe660d7c-2746-4b50-ab93-b2ed99960553",
        "classification": "CenterNotification",
        "entity_type": "aliquot",
        "created_datetime": "2015-09-28T00:00:00",
        "annotation_id": "5ddadefe-8b57-5ce2-b8b2-918d63d99a59",
        "notes": "The aliquot failed Broad pipeline QC and not all files are suitable for use. Consult the SDRF file to determine which files are usable.",
        "updated_datetime": "2017-03-09T13:20:38.962182-06:00",
        "submitter_id": "29087",
        "state": "submitted",
        "case_id": "41b59716-116f-4942-8b63-409870a87e26",
        "case_submitter_id": "TCGA-DK-A3IM",
        "entity_submitter_id": "TCGA-DK-A3IM-10A-01D-A20B-01",
        "id": "5ddadefe-8b57-5ce2-b8b2-918d63d99a59"
      },
      {
        "category": "Item is noncanonical",
        "status": "Approved",
        "entity_id": "25ebc29a-7598-4ae4-ba7f-618d448882cc",
        "classification": "Notification",
        "entity_type": "sample",
        "created_datetime": "2012-07-12T00:00:00",
        "annotation_id": "d6500f94-618f-5334-a810-ade76b887ec9",
        "notes": "No Matching Normal",
        "updated_datetime": "2017-03-09T13:47:18.182075-06:00",
        "submitter_id": "8009",
        "state": "submitted",
        "case_id": "bd114e05-5a97-41e2-a0d5-5d39a1e9d461",
        "case_submitter_id": "TCGA-08-0514",
        "entity_submitter_id": "TCGA-08-0514-01A",
        "id": "d6500f94-618f-5334-a810-ade76b887ec9"
      },
      {
        "category": "Prior malignancy",
        "status": "Approved",
        "entity_id": "e0d36cc0-652c-4224-bb10-09d15c7bd8f1",
        "classification": "Notification",
        "entity_type": "case",
        "created_datetime": "2013-03-12T00:00:00",
        "annotation_id": "33336cdf-2cf0-5af2-bb52-fecd3427f180",
        "notes": "Patient had a prior lymphoma. Unknown radiation or systemic chemotherapy.",
        "updated_datetime": "2017-03-09T12:11:31.786013-06:00",
        "submitter_id": "15630",
        "state": "submitted",
        "case_id": "e0d36cc0-652c-4224-bb10-09d15c7bd8f1",
        "case_submitter_id": "TCGA-FS-A1ZF",
        "entity_submitter_id": "TCGA-FS-A1ZF",
        "id": "33336cdf-2cf0-5af2-bb52-fecd3427f180"
      }
    ],
    "pagination": {
      "count": 3,
      "sort": "",
      "from": 0,
      "page": 1,
      "total": 3,
      "pages": 1,
      "size": 10
    }
  },
  "warnings": {}
}
```
### `History` Endpoint

The GDC History Endpoint `https://api.gdc.cancer.gov/history` enables search and retrieval of version and release information about a file.  This endpoint will return the entire provenance of all versions of a file.  A file may be versioned if a file is updated by the GDC (e.g. using a new alignment algorithm or fixing a file that contained an error). `Version` refers to the instance of a particular file. `Release` refers to which data release a file was part of.  A file may be a part of many different data releases with no change in version number or content.  

#### Example

This example is a query for versioning information associated with the follow with file `1dd28069-5777-4ff9-bd2b-d1ba68e88b06`.


```shell
curl 'https://api.gdc.cancer.gov/history/1dd28069-5777-4ff9-bd2b-d1ba68e88b06'
```
``` Output
[{"uuid": "1dd28069-5777-4ff9-bd2b-d1ba68e88b06", "version": "1", "file_change": "released", "release_date": "2018-08-23", "data_release": "12.0"}]
```


### `_mapping` Endpoint

Each search and retrieval endpoint is equipped with a ```_mapping``` endpoint that provides information about available fields. For example, `files/_mapping` endpoint provides information about fields and field groups available at the `files` endpoint: `https://api.gdc.cancer.gov/files/_mapping`.

The high-level structure of a response to a `_mapping` query is as follows:

	"_mapping": {}
	, "defaults": []
	, "expand": []
	, "fields": []
	, "multi": []
	, "nested": []

[//]: # (_)

Each part of the response is described below:

| Part | Description |
|------|-------------|
| `_mapping` | All available fields and their descriptions. The endpoint-agnostic field names provided here are compatible with the `filters` parameter but are not always compatible with the `fields` parameter |
| `defaults` | The default set of fields included in the API response when the `fields` parameter is not used in the request |
| `expand` | Field group names for use with the `expand` parameter |
| `fields` | All available fields in an endpoint-specific format that is compatible with both the `filters` and `fields` parameters |
| `multi` | GDC internal use |
| `nested` | Nested fields |


#### Example

```shell
curl 'https://api.gdc.cancer.gov/projects/_mapping'
```
```output
This output was put thought a json format application for easier viewability.
{
    ...

      "_mapping": {
        "projects.disease_type": {
          "doc_type": "projects",
          "field": "disease_type",
          "type": "id"
        },
        "projects.name": {
          "doc_type": "projects",
          "field": "name",
          "type": "id"
        }
      }

    ...

}
```

Similar information can be obtained using the `fields` parameter; `fields` queries provide additional information in the response, such as the name of the Elastic Search document (`doc_type`), the field name and the type of value. A list of supported types (such as `string`, `long`, `float`, ...) can be obtained from [Elastic Search Documentation](https://www.elastic.co/guide/en/elasticsearch/reference/current/mapping-types.html).


## Request Parameters

The GDC API supports the following search & retrieval request parameters:

Parameter | Default | Description
--------- | ------- | -----------
filters| null | Specifies search parameters
format | JSON | Specifies the API response format: JSON, XML, or TSV
pretty | false | Returns response with indentations and line breaks in a human-readable format
fields | null | Specifies which fields to include in the response
expand | null | Returns multiple related fields
size | 10 | Specifies the number of results to return
from   | 0 | Specifies the first record to return from a set of search results
sort | null | Specifies sorting for the search results
facets | null | Provides all existing values for a given field and the number of records having this value.


### Filters: Specifying the Query

The `filters` parameter enables passing of complex search queries to the GDC API. The parameter carries a query in the form of a JSON object.

#### Query Format

A `filters` query consists of an operator (or a nested set of operators) with a set of `field` and `value` operands.

The following `filters` query operators are supported by the GDC API:

| Operator | Description                                      | Number of Operands | Logic example                                                |
|----------|--------------------------------------------------|--------------------|--------------------------------------------------------------|
| =        | equals (string or number)                        | one                | gender = "female"                                            |
| !=       | does not equal (string or number)                | one                | project_id != "TARGET-AML"                                   |
| <        | less than (number)                               | one                | age at diagnosis < 90y                                       |
| <=       | less than or equal (number)                      | one                | age at diagnosis <= 17                                       |
| >        | greater than (number)                            | one                | age at diagnosis > 50                                        |
| >=       | greater than or equal (number)                   | one                | age at diagnosis >= 18                                       |
| is       | is (missing)                                     | one                | gender is missing                                            |
| not      | not (missing)                                    | one                | race not missing                                             |
| in       | matches a string or number in (a list)           | multiple           | primary_site in [Brain, Lung]                                |
| exclude  | does not match any strings or values in (a list) | multiple           | experimental_strategy exclude [WXS, WGS, "Genotyping array"] |
| and      | (operation1) and (operation2)                    | multiple           | {primary_site in [Brain, Lung]} and {gender = "female"}      |
| or       | (operation1) or (operation2)                     | multiple           | {project_id != "TARGET-AML"} or {age at diagnosis < 90y}     |

The `field` operand specifies a field that corresponds to a property defined in the [GDC Data Dictionary](../../Data_Dictionary/viewer.md). A list of supported fields is provided in [Appendix A](Appendix_A_Available_Fields.md); the list can also be accessed programmatically at the [_mapping endpoint](#95mapping-endpoint).

The `value` operand specifies the search terms. Users can get a list of available values for a specific property by making a call to the appropriate API endpoint using the `facets` parameter, e.g. `https://api.gdc.cancer.gov/v0/cases?facets=demographic.gender&size=0&pretty=true`. See [Facets](#facets) for details.

A simple query with a single operator looks like this:

	{
	    "op":"=",
	    "content":{
	        "field":"cases.demographic.gender",
	        "value":[
	            "male"
	        ]
	    }
	}

A more complex query with multiple operators looks like this:

	{
	    "op":"and",
	    "content":[
	        {
	            "op":"in",
	            "content":{
	                "field":"cases.submitter_id",
	                "value":[
	                    "TCGA-CK-4948",
	                    "TCGA-D1-A17N",
	                    "TCGA-4V-A9QX",
	                    "TCGA-4V-A9QM"
	                ]
	            }
	        },
	        {
	            "op":"=",
	            "content":{
	                "field":"files.data_type",
	                "value":"Gene Expression Quantification"
	            }
	        }
	    ]
	}


#### Example: HTTP GET Request

This example requests `male` cases using HTTP GET.

The JSON object to be passed to the GDC API looks like:

    {
        "op": "=",
        "content": {
            "field": "cases.demographic.gender",
            "value": [
                "male"
           ]
        }
    }

URL-encoding the above JSON object using [Percent-(URL)-encoding tool](https://www.freeformatter.com/url-encoder.html) results in the following string:

    %7B%0D%0A++++%22op%22%3A+%22%3D%22%2C%0D%0A++++%22content%22%3A+%7B%0D%0A++++++++%22field%22%3A+%22cases.demographic.gender%22%2C%0D%0A++++++++%22value%22%3A+%5B%0D%0A++++++++++++%22male%22%0D%0A++++++++%5D%0D%0A++++%7D%0D%0A%7D

The above string can now be passed to the GDC API using the `filters` parameter:

```shell
 curl  'https://api.gdc.cancer.gov/cases?filters=%7b%22op%22%3a+%22%3d%22%2c%0d%0a++++++%22content%22%3a+%7b%0d%0a++++++++++%22field%22%3a+%22cases.demographic.gender%22%2c%0d%0a++++++++++%22value%22%3a+%5b%22male%22%5d%0d%0a++++++%7d%0d%0a%7d&pretty=true'
```
```python
import requests
import json
cases_endpt = 'https://api.gdc.cancer.gov/cases'
filt = {"op":"=",
        "content":{
            "field": "cases.demographic.gender",
            "value": ["male"]
        }
}
params = {'filters':json.dumps(filt), 'sort':'demographic.gender:asc'}
# requests URL-encodes automatically
response = requests.get(cases_endpt, params = params)
print json.dumps(response.json(), indent=2)
```
``` Output
{
  "data": {
    "hits": [
      {
        "id": "f8970455-bfb2-4b1d-ab71-3c5d619898ad",
        "slide_ids": [
          "324684b5-8f18-4aa3-9b32-78382b96760b",
          "542a84f2-35e5-4843-9e98-c3d4bf0efe34"
        ],
        "submitter_slide_ids": [
          "TCGA-ZN-A9VQ-01Z-00-DX1",
          "TCGA-ZN-A9VQ-01A-01-TS1"
        ],
        "disease_type": "Mesothelial Neoplasms",
        "analyte_ids": [
          "73491451-b44e-41b3-be8c-1d7e44e54d08",
          "eda45603-7e65-4160-9963-e8907e7248b2",
          "f8f4c5d9-b09d-46d4-9fc1-afbebae1a81d"
        ],
        "submitter_id": "TCGA-ZN-A9VQ",
        "submitter_analyte_ids": [
          "TCGA-ZN-A9VQ-01A-11D",
          "TCGA-ZN-A9VQ-01A-11R",
          "TCGA-ZN-A9VQ-10A-01D"
        ],
        "submitter_aliquot_ids": [
          "TCGA-ZN-A9VQ-01A-11R-A40A-07",
          "TCGA-ZN-A9VQ-10A-01D-A39T-01",
          "TCGA-ZN-A9VQ-01A-11D-A39S-05",
          "TCGA-ZN-A9VQ-01A-11D-A40F-26",
          "TCGA-ZN-A9VQ-01A-11D-A39Q-01",
          "TCGA-ZN-A9VQ-01A-11D-A761-36",
          "TCGA-ZN-A9VQ-01A-11R-A404-13",
          "TCGA-ZN-A9VQ-01A-11D-A39R-32",
          "TCGA-ZN-A9VQ-10A-01D-A39U-32",
          "TCGA-ZN-A9VQ-10A-01D-A761-36",
          "TCGA-ZN-A9VQ-10A-01D-A40G-26"
        ],
        "aliquot_ids": [
          "26e89986-bfd5-4d4e-a3ff-5ad612c48358",
          "3d7835c0-388f-4cdc-9fe3-2dae11b71daa",
          "3f076394-cb86-46ea-8cdc-88f385c6b54e",
          "0b2cbe7c-8392-4072-95af-9ded20aa3888",
          "6762b668-a952-4c66-8a63-af909dbdc3ec",
          "dd85d29b-5883-4292-849a-8706698ff32b",
          "0e3a5bc1-9fe1-49e6-9ae8-fe7e86198acd",
          "4c1a209f-b005-4137-a344-c0befa66047c",
          "418d6b09-cffc-4ca4-84eb-7c9b2b5aacaf",
          "a09f740a-9529-4688-be12-978f13054e1e",
          "c32fb3d0-2894-4859-be14-30f03e0d8997"
        ],
        "created_datetime": null,
        "diagnosis_ids": [
          "68b5b7ee-bbbe-502c-b087-f325c4ccde09"
        ],
        "sample_ids": [
          "95830203-ccab-4a0d-8daf-e2b67ab95b86",
          "089c6901-5fe6-48b0-97ab-39f00609255c",
          "546ba1f1-7e16-4701-875a-8e9dd426fb76"
        ],
        "submitter_sample_ids": [
          "TCGA-ZN-A9VQ-10A",
          "TCGA-ZN-A9VQ-01Z",
          "TCGA-ZN-A9VQ-01A"
        ],
        "submitter_diagnosis_ids": [
          "TCGA-ZN-A9VQ_diagnosis"
        ],
        "primary_site": "Heart, mediastinum, and pleura",
        "updated_datetime": "2019-08-06T14:39:56.656272-05:00",
        "case_id": "f8970455-bfb2-4b1d-ab71-3c5d619898ad",
        "portion_ids": [
          "5b4b99dc-a44e-4679-8b19-d1d78020aa9f",
          "3e8db2c5-d5b3-4a9b-9ef5-948f16fe5cac",
          "ae23fda4-25ac-44ad-bfba-75e393b4bec5"
        ],
        "state": "released",
        "submitter_portion_ids": [
          "TCGA-ZN-A9VQ-01A-21-A45O-20",
          "TCGA-ZN-A9VQ-01A-11",
          "TCGA-ZN-A9VQ-10A-01"
        ]
      },
      {
        "id": "c739fd61-22b2-412d-bcf3-89bda45a2c0f",
        "slide_ids": [
          "6aa9b64a-9624-424e-87be-2ca9902794a2",
          "1c5ef953-b382-4b73-8bda-f2cfe8c86874"
        ],
        "submitter_slide_ids": [
          "TCGA-3H-AB3X-01A-01-TS1",
          "TCGA-3H-AB3X-01Z-00-DX1"
        ],
        "disease_type": "Mesothelial Neoplasms",
        "analyte_ids": [
          "0be4a69c-50e5-4634-b3ad-cf75f5cea8c5",
          "e208e95c-f999-41c2-b8df-6ee04af51f83",
          "7f0ed3c8-f3b2-47bf-a911-fe6463502315"
        ],
        "submitter_id": "TCGA-3H-AB3X",
        "submitter_analyte_ids": [
          "TCGA-3H-AB3X-01A-11R",
          "TCGA-3H-AB3X-10A-01D",
          "TCGA-3H-AB3X-01A-11D"
        ],
        "aliquot_ids": [
          "ba8caacf-7a47-48db-b375-58b2a417d073",
          "4026d79f-155e-48e9-91de-dedbf201f55a",
          "c895dbf4-5753-489a-b83a-dc7d80456388",
          "d0253b2c-99e0-49ff-b7f6-314bf1729cbf",
          "4a21ef3b-2d39-4588-8763-d2c26254932a",
          "537d8c34-d69a-4982-92d1-4d2d48e8b9c5",
          "41b7c449-47f2-470f-9f92-f3a4a02c8549",
          "3127573c-6700-4b82-9262-2bcc9c72b56c",
          "57dbf875-55a2-4a18-bccb-3a16df3ddbfa"
        ],
        "submitter_aliquot_ids": [
          "TCGA-3H-AB3X-10A-01D-A39U-32",
          "TCGA-3H-AB3X-01A-11D-A39S-05",
          "TCGA-3H-AB3X-01A-11R-A40A-07",
          "TCGA-3H-AB3X-01A-11D-A39Q-01",
          "TCGA-3H-AB3X-01A-11D-A40F-26",
          "TCGA-3H-AB3X-01A-11R-A404-13",
          "TCGA-3H-AB3X-01A-11D-A39R-32",
          "TCGA-3H-AB3X-10A-01D-A39T-01",
          "TCGA-3H-AB3X-10A-01D-A40G-26"
        ],
        "created_datetime": null,
        "diagnosis_ids": [
          "9252ba53-0bd5-5f49-b154-aa1dd6473fd5"
        ],
        "sample_ids": [
          "276c6d7b-712a-465c-913f-320de285cad4",
          "b06b0a8a-e992-42aa-a6aa-e85685bbe3f9",
          "ae48bc7d-fc5a-4d63-9399-783e1d7f50d6"
        ],
        "submitter_sample_ids": [
          "TCGA-3H-AB3X-01Z",
          "TCGA-3H-AB3X-10A",
          "TCGA-3H-AB3X-01A"
        ],
        "submitter_diagnosis_ids": [
          "TCGA-3H-AB3X_diagnosis"
        ],
        "primary_site": "Bronchus and lung",
        "updated_datetime": "2019-08-06T14:39:45.057305-05:00",
        "case_id": "c739fd61-22b2-412d-bcf3-89bda45a2c0f",
        "portion_ids": [
          "7ec81582-1cb3-4f34-a68f-23a9d9804658",
          "72bd54de-7e9f-46de-b62e-7994c9d4ad4d"
        ],
        "state": "released",
        "submitter_portion_ids": [
          "TCGA-3H-AB3X-10A-01",
          "TCGA-3H-AB3X-01A-11"
        ]
      },
      {
        "id": "ae90972d-5bc2-4b53-b5ff-1b8c31f39342",
        "slide_ids": [
          "e94c1249-1f03-4a7f-9161-d1714ef546b2",
          "e8bdc1d8-01c0-40a7-b308-cb5681f28b2f"
        ],
        "submitter_slide_ids": [
          "TCGA-CQ-A4CH-01A-01-TSA",
          "TCGA-CQ-A4CH-01Z-00-DX1"
        ],
        "disease_type": "Squamous Cell Neoplasms",
        "analyte_ids": [
          "be4729cb-4f7d-4c9b-8df7-2de2328c40d4",
          "53e0774f-2a36-4100-a73a-5fd8a45c4169",
          "51e4cd09-e408-4077-9ec3-7c805296a016",
          "d99212f6-8c23-44c4-8c30-87fc4166c9d0",
          "f3adfcba-e014-43e0-8c9b-09039b84613e"
        ],
        "submitter_id": "TCGA-CQ-A4CH",
        "submitter_analyte_ids": [
          "TCGA-CQ-A4CH-10A-01W",
          "TCGA-CQ-A4CH-01A-11D",
          "TCGA-CQ-A4CH-01A-11R",
          "TCGA-CQ-A4CH-10A-01D",
          "TCGA-CQ-A4CH-01A-11W"
        ],
        "submitter_aliquot_ids": [
          "TCGA-CQ-A4CH-01A-11R-A266-07",
          "TCGA-CQ-A4CH-10A-01D-A25X-01",
          "TCGA-CQ-A4CH-01A-11R-A25Z-13",
          "TCGA-CQ-A4CH-01A-11W-A296-08",
          "TCGA-CQ-A4CH-10A-01W-A296-08",
          "TCGA-CQ-A4CH-10A-01D-A25Y-08",
          "TCGA-CQ-A4CH-01A-11D-A25X-01",
          "TCGA-CQ-A4CH-01A-11D-A265-05",
          "TCGA-CQ-A4CH-01A-11D-A25Y-08"
        ],
        "aliquot_ids": [
          "9eb34130-b3b9-4c92-8238-4053f8c6d06b",
          "0c99ef22-b1ef-42c9-a184-4c98f674c7be",
          "b77aa7fa-403b-4ee1-b537-695a799c80f5",
          "3aeb1426-58f0-47bd-825e-8d1578dda18b",
          "afd797fb-97d2-4482-ae15-753f5f66b828",
          "73074dac-057b-456d-b0b8-761b149d10dc",
          "48bbbc5e-34c0-46be-a3bb-2d28c2e7d357",
          "818db1ea-5ae7-4d50-87cf-a85c62830566",
          "92e8b340-322e-4411-95d8-7db81767f660"
        ],
        "created_datetime": null,
        "diagnosis_ids": [
          "7e574325-c82b-5054-bf4a-038f528c4110"
        ],
        "sample_ids": [
          "273490d1-9862-4480-86aa-12522b35fe24",
          "c7dc33be-768c-42a2-a1ab-2eb5b67be87e",
          "465d6949-f5a4-4a9d-a082-1f1403876d85"
        ],
        "submitter_sample_ids": [
          "TCGA-CQ-A4CH-01Z",
          "TCGA-CQ-A4CH-01A",
          "TCGA-CQ-A4CH-10A"
        ],
        "submitter_diagnosis_ids": [
          "TCGA-CQ-A4CH_diagnosis"
        ],
        "primary_site": "Other and unspecified parts of tongue",
        "updated_datetime": "2019-08-06T14:25:53.026261-05:00",
        "case_id": "ae90972d-5bc2-4b53-b5ff-1b8c31f39342",
        "portion_ids": [
          "3bf8dee3-ac07-4ab2-bb1c-c8b564e6b3e1",
          "8c802162-a804-437a-9f54-93de7a4c21b3",
          "5131b158-ca97-4afe-b236-5d6fc41f70fd"
        ],
        "state": "released",
        "submitter_portion_ids": [
          "TCGA-CQ-A4CH-01A-11",
          "TCGA-CQ-A4CH-10A-01",
          "TCGA-CQ-A4CH-01A-21-A45L-20"
        ]
      },
      {
        "id": "24c1cf70-bd67-431e-a623-d20f8d3f52b2",
        "slide_ids": [
          "0fb2f319-414d-4bd4-bd88-e61303208dfb",
          "847e5a3f-c88a-40a3-926e-563de0a26ca0",
          "15737384-02ba-47f3-8655-c11ba975edc5",
          "50044cf9-fee9-4745-8904-8f5240d4d18d"
        ],
        "submitter_slide_ids": [
          "TCGA-CJ-4881-01Z-00-DX1",
          "TCGA-CJ-4881-01A-01-BS1",
          "TCGA-CJ-4881-11A-01-TS1",
          "TCGA-CJ-4881-01A-01-TS1"
        ],
        "disease_type": "Adenomas and Adenocarcinomas",
        "analyte_ids": [
          "2f46942c-7f9c-4dab-ac3b-fb5abf0ca11c",
          "986f008d-d05c-4dda-88fc-858bf4f1f379",
          "e90c3295-7c25-4e34-9b89-eeb6838081e4",
          "2fe03fe3-1036-40d2-81ea-9933be17e4d9",
          "020b312a-7b12-49b0-9135-189eb5e1e42a",
          "21fddcbd-2f48-4733-845e-089a4e9076a7",
          "4a25c9f2-f1ff-4cb2-ba9a-77e55961c488"
        ],
        "submitter_id": "TCGA-CJ-4881",
        "submitter_analyte_ids": [
          "TCGA-CJ-4881-11A-01X",
          "TCGA-CJ-4881-01A-01R",
          "TCGA-CJ-4881-11A-01D",
          "TCGA-CJ-4881-01A-01X",
          "TCGA-CJ-4881-01A-01W",
          "TCGA-CJ-4881-01A-01D",
          "TCGA-CJ-4881-11A-01W"
        ],
        "submitter_aliquot_ids": [
          "TCGA-CJ-4881-01A-01D-1303-05",
          "TCGA-CJ-4881-01A-01R-1305-07",
          "TCGA-CJ-4881-01A-01X-1371-10",
          "TCGA-CJ-4881-11A-01X-1371-10",
          "TCGA-CJ-4881-01A-01R-1762-13",
          "TCGA-CJ-4881-01A-01D-2098-10",
          "TCGA-CJ-4881-01A-01W-1369-10",
          "TCGA-CJ-4881-01A-01R-1304-13",
          "TCGA-CJ-4881-11A-01D-1302-01",
          "TCGA-CJ-4881-11A-01D-2098-10",
          "TCGA-CJ-4881-11A-01D-1303-05",
          "TCGA-CJ-4881-11A-01W-1369-10",
          "TCGA-CJ-4881-01A-01D-1373-10",
          "TCGA-CJ-4881-01A-01D-1301-02",
          "TCGA-CJ-4881-11A-01D-1301-02",
          "TCGA-CJ-4881-11A-01D-1373-10",
          "TCGA-CJ-4881-01A-01D-1302-01"
        ],
        "aliquot_ids": [
          "581cae5b-a55e-4bd4-a9da-48b5480615c0",
          "00d54c43-aba8-4503-827a-30444e38c704",
          "290dd57c-0f01-431d-8b72-5f25f1a00ca7",
          "63f47ca8-5a98-4c54-8d83-9f0f9c9f4559",
          "ece4609c-ffe3-4330-a92b-c4847b618d77",
          "30ec4e95-d07b-4a2d-9869-7b9008eb8d8b",
          "e373d58b-dccf-49e3-ab23-c177883ea2bc",
          "521d43d9-882f-42dd-ad21-ef5b6df5dae9",
          "67d59b58-b34c-4616-8661-b58cbb32e726",
          "9efca5db-e210-468a-a38e-9fcc83d3f113",
          "26281187-0846-4578-8e8a-ad9886493af7",
          "495dbab8-1e78-4cb4-b3e4-0ffda17c823a",
          "ae20e2d0-5d39-4a94-a9ff-dee71503cbfe",
          "d883e20c-5237-420c-a58b-98ca359f6b2a",
          "ec4d0eff-cbe4-4dbb-8319-1f4b0b4a5d35",
          "5d2fada9-5a0f-41b9-b602-9675623191ca",
          "ffbf81f1-b5a5-4739-9f35-716169650023"
        ],
        "created_datetime": null,
        "diagnosis_ids": [
          "f9c1962f-13de-5ffe-b201-9d125adc590f"
        ],
        "sample_ids": [
          "f00de7e8-d54d-4f07-85bc-8fa3f1d56b0c",
          "1f8ce8ea-74b6-4d91-81e7-b689d11c26cd",
          "534a3864-e8f4-462f-81a5-c0a3895cb68f"
        ],
        "submitter_sample_ids": [
          "TCGA-CJ-4881-01Z",
          "TCGA-CJ-4881-01A",
          "TCGA-CJ-4881-11A"
        ],
        "submitter_diagnosis_ids": [
          "TCGA-CJ-4881_diagnosis"
        ],
        "primary_site": "Kidney",
        "updated_datetime": "2019-08-06T14:29:28.932622-05:00",
        "case_id": "24c1cf70-bd67-431e-a623-d20f8d3f52b2",
        "portion_ids": [
          "4f1a0956-ee24-4ba1-b8a3-b23fb8a26601",
          "1a7f6ad8-97f3-416d-8d1e-1752bb6638f7",
          "367e98f4-181c-4656-b28e-855fa6f265af"
        ],
        "state": "released",
        "submitter_portion_ids": [
          "TCGA-CJ-4881-11A-01",
          "TCGA-CJ-4881-01A-01",
          "TCGA-CJ-4881-01A-21-1739-20"
        ]
      },
      {
        "id": "12adefc4-c9fd-46d3-904b-8fc52d5f1913",
        "slide_ids": [
          "802aa0c2-d722-46b1-9a66-29333abba34c",
          "6a956f2d-0040-498f-b0e6-68839fe8f668",
          "3d0c023e-7023-4e11-931f-61c72105d828",
          "cbd62e61-cc59-4533-bca2-518ce48c1ab3",
          "d56c89b6-cda6-45f7-965a-1cf6dc037ad1"
        ],
        "submitter_slide_ids": [
          "TCGA-B0-4696-01A-01-TS1",
          "TCGA-B0-4696-11A-01-TS1",
          "TCGA-B0-4696-11A-01-BS1",
          "TCGA-B0-4696-01A-01-BS1",
          "TCGA-B0-4696-01Z-00-DX1"
        ],
        "disease_type": "Adenomas and Adenocarcinomas",
        "analyte_ids": [
          "a5399e3a-54ba-4971-b72a-a7bea9312511",
          "41488c12-a3cd-4118-bc9d-abf8d28eac27",
          "3da18be3-3bf7-475f-947e-5b2c22e7c046",
          "6bb85779-d0a1-47ca-aacc-540c907604f4",
          "06fd554d-32bb-4c13-aec8-33619b77fea6",
          "474d8e2f-76b8-4cb3-b166-62266be78c42",
          "26bf50f5-c51a-4c06-b68a-c23f271d0077"
        ],
        "submitter_id": "TCGA-B0-4696",
        "submitter_analyte_ids": [
          "TCGA-B0-4696-01A-01W",
          "TCGA-B0-4696-01A-01X",
          "TCGA-B0-4696-11A-01D",
          "TCGA-B0-4696-11A-01W",
          "TCGA-B0-4696-01A-01R",
          "TCGA-B0-4696-11A-01X",
          "TCGA-B0-4696-01A-01D"
        ],
        "aliquot_ids": [
          "3c57f0c5-340c-478f-af17-55c0bda50cce",
          "1701e9f3-22ba-4fd9-9a88-d310a3c9dc1e",
          "c9c0689c-da29-476a-84bb-f52ae5698b50",
          "f04d25de-d349-4ed4-94c7-aaf24528e77b",
          "f00ce2ee-45c4-4232-8758-d4de6d85ea46",
          "5ddfdefe-9b9d-438b-b55b-2953a2541378",
          "504da228-bbe6-490f-9f59-94f6f00adb31",
          "e9aedb52-652a-48b2-a9ef-e3a1010e5874",
          "e3847aa2-047f-46ba-bb5c-e18d725defc2",
          "92590459-3a3e-4160-9c57-30e078258103",
          "57dac8d5-58bf-4beb-ad24-d34437d1b3e7",
          "9e5b6046-1d2c-44af-8430-421ec0036bb6",
          "a744ca75-7994-453c-a283-3db5778844f9",
          "3128f70f-0a05-41dc-8657-2312f8af4104",
          "70544bed-6e2a-45f0-9b7c-e4f44705935f",
          "65197a42-c28f-47ee-909f-64555fb8478e"
        ],
        "submitter_aliquot_ids": [
          "TCGA-B0-4696-01A-01R-1277-07",
          "TCGA-B0-4696-01A-01D-1275-05",
          "TCGA-B0-4696-11A-01D-1275-05",
          "TCGA-B0-4696-11A-01D-1274-01",
          "TCGA-B0-4696-01A-01D-2096-10",
          "TCGA-B0-4696-11A-01D-1273-02",
          "TCGA-B0-4696-01A-01X-1360-10",
          "TCGA-B0-4696-01A-01D-1361-10",
          "TCGA-B0-4696-01A-01D-1274-01",
          "TCGA-B0-4696-11A-01X-1360-10",
          "TCGA-B0-4696-11A-01W-1359-10",
          "TCGA-B0-4696-01A-01W-1359-10",
          "TCGA-B0-4696-01A-01D-1273-02",
          "TCGA-B0-4696-11A-01D-1361-10",
          "TCGA-B0-4696-11A-01D-2096-10",
          "TCGA-B0-4696-01A-01R-1276-13"
        ],
        "created_datetime": null,
        "diagnosis_ids": [
          "34b6cdff-d3b0-55e6-9e92-4ff8e9f7803f"
        ],
        "sample_ids": [
          "1ee83ac8-c90b-4364-9e4f-2ed76978ccda",
          "1784f393-cf44-4bba-a0fd-98d0f6c50b83",
          "9e30a07a-f227-4b6a-9c80-619abed6f7bd"
        ],
        "submitter_sample_ids": [
          "TCGA-B0-4696-01A",
          "TCGA-B0-4696-11A",
          "TCGA-B0-4696-01Z"
        ],
        "submitter_diagnosis_ids": [
          "TCGA-B0-4696_diagnosis"
        ],
        "primary_site": "Kidney",
        "updated_datetime": "2019-08-06T14:28:04.080572-05:00",
        "case_id": "12adefc4-c9fd-46d3-904b-8fc52d5f1913",
        "portion_ids": [
          "1f4528b1-44c0-459d-a71e-b3c11575b6f4",
          "e8148187-d209-4ee6-a261-c2860402f02f",
          "2f3b0ffd-223c-4573-86e4-afdef6116698"
        ],
        "state": "released",
        "submitter_portion_ids": [
          "TCGA-B0-4696-01A-02-1738-20",
          "TCGA-B0-4696-01A-01",
          "TCGA-B0-4696-11A-01"
        ]
      },
      {
        "id": "c0edde5e-d229-4061-8820-14afc712c5b6",
        "slide_ids": [
          "b57e7804-588d-4e71-909b-257fa877e2ad",
          "fc91a512-85a5-4479-ab07-9b8bd96fdf91"
        ],
        "submitter_slide_ids": [
          "TCGA-HD-A633-01A-01-TS1",
          "TCGA-HD-A633-01Z-00-DX1"
        ],
        "disease_type": "Squamous Cell Neoplasms",
        "analyte_ids": [
          "c8cd9c1f-6688-4c3a-9414-ebede9a06622",
          "cac22161-1b5d-4436-8ba3-12fba6aae76f",
          "6b410bfe-a420-4268-a2a3-f3aac24e8761",
          "22a7918c-6a39-444e-b995-2ba9e60eae90",
          "a07ff300-0c9c-4b70-88a8-b4ac22936e0d"
        ],
        "submitter_id": "TCGA-HD-A633",
        "submitter_analyte_ids": [
          "TCGA-HD-A633-10A-01D",
          "TCGA-HD-A633-10A-01W",
          "TCGA-HD-A633-01A-11W",
          "TCGA-HD-A633-01A-11R",
          "TCGA-HD-A633-01A-11D"
        ],
        "submitter_aliquot_ids": [
          "TCGA-HD-A633-01A-11W-A316-08",
          "TCGA-HD-A633-01A-11D-A28R-08",
          "TCGA-HD-A633-10A-01D-A28U-08",
          "TCGA-HD-A633-01A-11R-A28V-07",
          "TCGA-HD-A633-01A-11D-A28Q-01",
          "TCGA-HD-A633-01A-11R-A28Z-13",
          "TCGA-HD-A633-10A-01D-A28T-01",
          "TCGA-HD-A633-10A-01W-A317-08",
          "TCGA-HD-A633-01A-11D-A28S-05"
        ],
        "aliquot_ids": [
          "19657bce-8994-4dc6-8923-6d3959874dec",
          "10a6dce3-7446-4dc5-8df0-b36558aa8449",
          "0eb29e0f-dff7-47fe-a184-f4d3247298dc",
          "b492257f-17a5-48de-aee3-a21a9bedf5ee",
          "db856f4b-d854-4bab-b92e-755dfdb8acea",
          "d67c2939-353f-4299-88ea-e3d640fc4ac8",
          "31641d8f-733a-469b-ae84-907430e8bf0a",
          "36aead4f-8d43-426b-ad2c-7dcb249f579e",
          "895828e4-db59-4f6a-9df3-c423ec3ec6b7"
        ],
        "created_datetime": null,
        "diagnosis_ids": [
          "3c27a0ae-0496-5d58-b9c6-1b6362b5da05"
        ],
        "sample_ids": [
          "ec6326e4-1a3e-4afd-ae55-6a9fdb872d93",
          "6ad70a8a-2d7e-406a-a419-b805b72eba56",
          "1b65e4f9-45cb-47ed-bbae-e3f189b1c8f8"
        ],
        "submitter_sample_ids": [
          "TCGA-HD-A633-10A",
          "TCGA-HD-A633-01A",
          "TCGA-HD-A633-01Z"
        ],
        "submitter_diagnosis_ids": [
          "TCGA-HD-A633_diagnosis"
        ],
        "primary_site": "Other and unspecified parts of mouth",
        "updated_datetime": "2019-08-06T14:26:51.527876-05:00",
        "case_id": "c0edde5e-d229-4061-8820-14afc712c5b6",
        "portion_ids": [
          "dd7bf6a5-9e94-470a-8fb4-9f2a2a56c173",
          "d4475b19-78d8-4b71-ac0e-86a4230aa0cd",
          "9c1c0a78-2d57-4157-bff1-00d5bb818c26"
        ],
        "state": "released",
        "submitter_portion_ids": [
          "TCGA-HD-A633-01A-11",
          "TCGA-HD-A633-10A-01",
          "TCGA-HD-A633-01A-21-A45M-20"
        ]
      },
      {
        "id": "ae2051a3-f851-4cb3-af18-029f4574b179",
        "slide_ids": [
          "90f06679-0397-4688-8c75-b1855efa4bd1",
          "f8c610bf-7789-4b57-8a74-cf29106be142",
          "30e6a59b-fdc9-4976-9040-e1a851232b1b",
          "e6b71383-bd0f-4f2f-9483-f0048522c234"
        ],
        "submitter_slide_ids": [
          "TCGA-CV-7245-01A-01-TS1",
          "TCGA-CV-7245-01A-01-BS1",
          "TCGA-CV-7245-01Z-00-DX1",
          "TCGA-CV-7245-11A-01-TS1"
        ],
        "disease_type": "Squamous Cell Neoplasms",
        "analyte_ids": [
          "4215b417-de72-49da-87f5-efee18f6166f",
          "2acff637-98f6-454d-8217-70d42e5612c1",
          "a55e555d-4858-4153-9047-4ce50794efe8",
          "b30b89ec-e0f2-49b8-b47d-489730000460",
          "b8c80ab1-d3f7-43da-8968-6c846f03724d",
          "17ae48ec-047f-4cca-8182-a2635ce772be",
          "f8115c36-56c0-4e68-bc96-a26331972298",
          "d4c8a05a-a194-48fb-8c5d-2f9391613dcb"
        ],
        "submitter_id": "TCGA-CV-7245",
        "submitter_analyte_ids": [
          "TCGA-CV-7245-11A-01R",
          "TCGA-CV-7245-01A-11R",
          "TCGA-CV-7245-01A-11D",
          "TCGA-CV-7245-10A-01D",
          "TCGA-CV-7245-11A-01D",
          "TCGA-CV-7245-11A-01W",
          "TCGA-CV-7245-01A-11W",
          "TCGA-CV-7245-10A-01W"
        ],
        "submitter_aliquot_ids": [
          "TCGA-CV-7245-01A-11R-2015-13",
          "TCGA-CV-7245-11A-01R-2015-13",
          "TCGA-CV-7245-10A-01W-2033-08",
          "TCGA-CV-7245-11A-01D-2010-01",
          "TCGA-CV-7245-11A-01R-2016-07",
          "TCGA-CV-7245-01A-11D-2014-05",
          "TCGA-CV-7245-10A-01D-2011-01",
          "TCGA-CV-7245-11A-01D-2014-05",
          "TCGA-CV-7245-01A-11R-2016-07",
          "TCGA-CV-7245-11A-01D-2012-08",
          "TCGA-CV-7245-10A-01D-2013-08",
          "TCGA-CV-7245-01A-11D-2012-08",
          "TCGA-CV-7245-11A-01W-2033-08",
          "TCGA-CV-7245-01A-11W-2032-08",
          "TCGA-CV-7245-01A-11D-2010-01"
        ],
        "aliquot_ids": [
          "c1a21df1-8265-4040-843a-2785b98b3463",
          "95b2722d-f9c7-4d7f-8a70-7603080ecac5",
          "265cfe54-10da-45c9-b0d1-2d5c65308be1",
          "d73e8c20-b392-475e-936b-10d6bb173e6c",
          "ce0c0c2e-902a-4ecb-b753-e696750ca407",
          "0ead3f08-6ce9-4eea-8ac2-c8021f427f2b",
          "0f19fa35-7bf4-4cec-b463-8b1e01dfefa0",
          "56291b3c-595c-4388-a264-9037a48401d8",
          "f84f8f63-06bf-40c7-8ba8-439aaf360c89",
          "bdc58b10-1756-45bf-a34f-fa4095429476",
          "081cec95-360f-4127-9fbc-e6bfd212c23f",
          "61897747-36dc-4aad-b17a-fe935c605ca3",
          "f28ead58-aa1c-490b-8b4c-e7dc2d76a563",
          "b87cdd87-3663-4e8a-99f2-7d475f939cb6",
          "fe4ca5c3-0767-4264-be28-2bb9891c0c02"
        ],
        "created_datetime": null,
        "diagnosis_ids": [
          "cceca036-354c-59e4-ad7a-886865aada46"
        ],
        "sample_ids": [
          "aaa6607c-e333-4f99-b1dc-390c3cad4f1b",
          "c8ac9749-82ff-4348-8d91-5ce677697a15",
          "f0684e03-2d91-4e9b-9d11-a3063ccc2f34",
          "d455327a-a947-4a67-8fb5-44438abe5976"
        ],
        "submitter_sample_ids": [
          "TCGA-CV-7245-10A",
          "TCGA-CV-7245-01A",
          "TCGA-CV-7245-01Z",
          "TCGA-CV-7245-11A"
        ],
        "submitter_diagnosis_ids": [
          "TCGA-CV-7245_diagnosis"
        ],
        "primary_site": "Larynx",
        "updated_datetime": "2019-08-06T14:26:16.536997-05:00",
        "case_id": "ae2051a3-f851-4cb3-af18-029f4574b179",
        "portion_ids": [
          "9d94c9a4-9d71-4843-86f7-be076d939415",
          "90abd8dd-a40a-4433-9ae1-856021300eb4",
          "3cfc5361-e915-4237-bc4d-a2b06f48cc3e",
          "8bb93715-41cd-436a-8ee5-a35ea050b6c0"
        ],
        "state": "released",
        "submitter_portion_ids": [
          "TCGA-CV-7245-01A-13-2074-20",
          "TCGA-CV-7245-11A-01",
          "TCGA-CV-7245-10A-01",
          "TCGA-CV-7245-01A-11"
        ]
      },
      {
        "id": "57959b73-534d-450c-b09d-f70ef1ffee25",
        "slide_ids": [
          "f7e2faec-0fd2-41f2-82dc-8a888a1c01a5",
          "c9a0789f-cb9b-48f5-a8ce-489071fd618c",
          "bb49cb73-66de-4c20-a6ef-8483da3fc0e2",
          "f09ecda2-cab1-4f2b-a616-26f9ce4fd479",
          "09b3f4a8-1366-41a0-b2ea-2dc010f431f8"
        ],
        "submitter_slide_ids": [
          "TCGA-BP-4346-01A-01-TS1",
          "TCGA-BP-4346-11A-01-BS1",
          "TCGA-BP-4346-01A-01-BS1",
          "TCGA-BP-4346-01Z-00-DX1",
          "TCGA-BP-4346-11A-01-TS1"
        ],
        "disease_type": "Adenomas and Adenocarcinomas",
        "analyte_ids": [
          "ce37810b-cf9e-4c0c-962f-432047dae3e7",
          "ac583369-7b4e-4b96-ae8e-ae5769c075a5",
          "5d9e9e98-e84d-4e33-8c73-b2ea951bbdbd",
          "4cf5f5d3-436d-4f40-9beb-c3d2299c008b",
          "a2debb38-75ec-456c-944a-d90aaa9fa032",
          "f75a4f5c-5d8b-423e-b2b8-718b818a025b",
          "8ab159b0-57f8-44b5-916e-3d0d709ffc38"
        ],
        "submitter_id": "TCGA-BP-4346",
        "submitter_analyte_ids": [
          "TCGA-BP-4346-01A-01W",
          "TCGA-BP-4346-01A-01D",
          "TCGA-BP-4346-11A-01D",
          "TCGA-BP-4346-11A-01X",
          "TCGA-BP-4346-01A-01X",
          "TCGA-BP-4346-01A-01R",
          "TCGA-BP-4346-11A-01W"
        ],
        "submitter_aliquot_ids": [
          "TCGA-BP-4346-11A-01X-1364-10",
          "TCGA-BP-4346-11A-01D-1283-01",
          "TCGA-BP-4346-01A-01D-1283-01",
          "TCGA-BP-4346-01A-01D-1282-02",
          "TCGA-BP-4346-01A-01D-1284-05",
          "TCGA-BP-4346-01A-01X-1364-10",
          "TCGA-BP-4346-11A-01D-1366-10",
          "TCGA-BP-4346-01A-01D-1366-10",
          "TCGA-BP-4346-01A-01W-1362-10",
          "TCGA-BP-4346-11A-01D-1284-05",
          "TCGA-BP-4346-11A-01D-2097-10",
          "TCGA-BP-4346-01A-01R-1289-07",
          "TCGA-BP-4346-01A-01R-1288-13",
          "TCGA-BP-4346-01A-01D-2097-10",
          "TCGA-BP-4346-11A-01D-1282-02",
          "TCGA-BP-4346-11A-01W-1362-10"
        ],
        "aliquot_ids": [
          "81ad330a-2e60-4a6e-b301-6564875f41c4",
          "b3ceb69a-6e16-41cf-9aa1-539b48aacca2",
          "956d13ed-7b60-4e15-a6d2-adfc3ecff4f8",
          "f1b05fbc-cda4-4015-a496-6d30c592fa3d",
          "e680ec7b-d88c-4da8-b339-da043fda3dc6",
          "a8b5426b-ff3a-435f-bffc-90bb4a06d77d",
          "20735253-1ba8-41c9-b4ab-682c1af79b9d",
          "391da188-dd73-40a8-afae-2f523354b95a",
          "60bdedba-bdea-4734-b7ff-a171f274e7d4",
          "109292f2-4bc0-4eff-8e5f-c829f51b835b",
          "753e5d3b-0334-4202-8af4-45bc9aad100e",
          "f58cda65-eaf7-4f4b-af5a-5c05618814d3",
          "d758addf-a866-45b6-819a-159f814206cb",
          "4f350288-6283-4f87-b1b4-6afae2edfa85",
          "69e027cf-071a-4b3e-97de-c450555f38c1",
          "0ca5fa21-f692-40cc-b501-cb3b670530b8"
        ],
        "created_datetime": null,
        "diagnosis_ids": [
          "c1d7ae9c-c61e-5509-aa35-7555f1d35493"
        ],
        "sample_ids": [
          "a6180591-a96d-4d38-9937-5963fac6b3df",
          "18c6c530-51e7-4e43-a1b9-f152d17c7c47",
          "5306e996-869d-4cee-a990-dea934af0cb8"
        ],
        "submitter_sample_ids": [
          "TCGA-BP-4346-01Z",
          "TCGA-BP-4346-11A",
          "TCGA-BP-4346-01A"
        ],
        "submitter_diagnosis_ids": [
          "TCGA-BP-4346_diagnosis"
        ],
        "primary_site": "Kidney",
        "updated_datetime": "2019-08-06T14:28:51.268056-05:00",
        "case_id": "57959b73-534d-450c-b09d-f70ef1ffee25",
        "portion_ids": [
          "8092a672-1506-4e71-a75d-ee15c8f4cbc4",
          "386408e6-2197-4edc-a459-21bb444fbc7f",
          "28206fd9-253b-4f2b-94f1-1e558edd32ed"
        ],
        "state": "released",
        "submitter_portion_ids": [
          "TCGA-BP-4346-01A-03-1737-20",
          "TCGA-BP-4346-11A-01",
          "TCGA-BP-4346-01A-01"
        ]
      },
      {
        "id": "714496c5-d221-4397-9c5a-cd2d22603e6f",
        "slide_ids": [
          "db5c34b9-d2c2-4bc7-a4e4-eb47eb806a1b",
          "80e8f54c-2976-4dcb-a9c3-0954e841aec3",
          "9c045280-8ea3-4ca5-a20b-4483dacc4c5c"
        ],
        "submitter_slide_ids": [
          "TCGA-CN-6013-01A-01-TS1",
          "TCGA-CN-6013-01Z-00-DX1",
          "TCGA-CN-6013-01A-01-BS1"
        ],
        "disease_type": "Squamous Cell Neoplasms",
        "analyte_ids": [
          "8a90c5f9-c866-4109-8058-b7da2ae79d4a",
          "3eeca312-02c3-4e20-898d-04e3b8aa0aad",
          "1b1ac5e8-63e9-4f9f-9231-00e0f91ae576",
          "016a6a0b-6d11-4074-a8b4-98c74599a586",
          "17b12617-8570-4d5f-85d6-39dc5a47dca6"
        ],
        "submitter_id": "TCGA-CN-6013",
        "submitter_analyte_ids": [
          "TCGA-CN-6013-01A-11R",
          "TCGA-CN-6013-01A-11W",
          "TCGA-CN-6013-10A-01W",
          "TCGA-CN-6013-10A-01D",
          "TCGA-CN-6013-01A-11D"
        ],
        "submitter_aliquot_ids": [
          "TCGA-CN-6013-10A-01D-1683-08",
          "TCGA-CN-6013-01A-11W-1767-08",
          "TCGA-CN-6013-01A-11D-1681-02",
          "TCGA-CN-6013-01A-11D-1683-08",
          "TCGA-CN-6013-01A-11D-1684-05",
          "TCGA-CN-6013-10A-01D-1681-02",
          "TCGA-CN-6013-01A-11R-1685-13",
          "TCGA-CN-6013-10A-01D-1682-01",
          "TCGA-CN-6013-01A-11R-1686-07",
          "TCGA-CN-6013-01A-11D-1682-01",
          "TCGA-CN-6013-10A-01W-1767-08"
        ],
        "aliquot_ids": [
          "e4a0cd96-81aa-41c8-b32c-1b6c57a08679",
          "992de9b5-c394-48e7-b4e3-4c4aeacb4a23",
          "0f8d7f78-2640-444d-a068-a348f50cde8f",
          "2de7dbc6-2a3e-4e67-bb96-552e27137618",
          "e996c35e-3614-4258-b99e-0eb444398931",
          "78c2da9e-3837-4098-b662-3c10d8b5d14e",
          "157aba21-4c81-45b9-914c-925813c537f6",
          "762ca7da-6d6e-4709-b1aa-13621d962d35",
          "1906f4b3-d226-4570-bbce-f2cbbb1d5ec2",
          "9496355f-6aa6-4f6a-ac06-b17966b4208c",
          "ffff755e-285f-466c-85c6-b372f4a7ae14"
        ],
        "created_datetime": null,
        "diagnosis_ids": [
          "7747a88f-1a9b-5833-9442-af0b1350bb4f"
        ],
        "sample_ids": [
          "c0c84016-117e-498c-9ed9-52e279e89e33",
          "d018cd44-6b01-44f2-9181-81153225aab8",
          "8b510bc8-7b3c-4861-9a40-c2c30b7099e0"
        ],
        "submitter_sample_ids": [
          "TCGA-CN-6013-10A",
          "TCGA-CN-6013-01A",
          "TCGA-CN-6013-01Z"
        ],
        "submitter_diagnosis_ids": [
          "TCGA-CN-6013_diagnosis"
        ],
        "primary_site": "Gum",
        "updated_datetime": "2019-08-06T14:25:25.511101-05:00",
        "case_id": "714496c5-d221-4397-9c5a-cd2d22603e6f",
        "portion_ids": [
          "4add0374-495b-4d4c-b96b-8494600fea26",
          "436822b9-08d5-493e-b213-b9b5dfc161e9",
          "acdd7eae-c26c-46c9-9dca-036ed3317718"
        ],
        "state": "released",
        "submitter_portion_ids": [
          "TCGA-CN-6013-10A-01",
          "TCGA-CN-6013-01A-11",
          "TCGA-CN-6013-01A-13-2072-20"
        ]
      },
      {
        "id": "74f89bda-bc2a-4cd3-9aa3-f89aa7282dd5",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Myeloid Leukemias",
        "submitter_id": "TARGET-20-PAUXKI",
        "submitter_aliquot_ids": [
          "TARGET-20-PAUXKI-09A-01R"
        ],
        "aliquot_ids": [
          "192d5f48-ba11-44cf-91d5-5eec5daea281"
        ],
        "sample_ids": [
          "22b9ad16-7791-4c41-a1c5-0392b4d4721c"
        ],
        "created_datetime": "2019-02-25T10:13:06.478422-06:00",
        "diagnosis_ids": [
          "2f7f2b8c-35fa-412f-894a-6c4b09cf8077"
        ],
        "submitter_sample_ids": [
          "TARGET-20-PAUXKI-09A"
        ],
        "primary_site": "Hematopoietic and reticuloendothelial systems",
        "submitter_diagnosis_ids": [
          "TARGET-20-PAUXKI_diagnosis"
        ],
        "updated_datetime": "2019-10-24T08:22:10.208559-05:00",
        "case_id": "74f89bda-bc2a-4cd3-9aa3-f89aa7282dd5",
        "index_date": null,
        "state": "released"
      }
    ],
    "pagination": {
      "count": 10,
      "total": 39092,
      "size": 10,
      "from": 0,
      "sort": "",
      "page": 1,
      "pages": 3910
    }
  },
  "warnings": {}
}
```



#### Example: HTTP POST Request

This example demonstrates how to obtain metadata in TSV format for a set of files using their UUIDs (e.g. UUIDs obtained from a [download manifest file generated by the GDC Data Portal](/Data_Portal/Users_Guide/Cart/#gdc-data-transfer-tool)).

The first step is to construct a JSON query object, including `filters`, `fields`, `format`, and `size` parameters. The object is then submitted as HTTP POST payload to the GDC API using curl, in order to retrieve a TSV file with the requested metadata.

```Payload_txt
{
    "filters":{
        "op":"in",
        "content":{
            "field":"files.file_id",
            "value":[
                "0001801b-54b0-4551-8d7a-d66fb59429bf",
                "002c67f2-ff52-4246-9d65-a3f69df6789e",
                "003143c8-bbbf-46b9-a96f-f58530f4bb82",
                "0043d981-3c6b-463f-b512-ab1d076d3e62",
                "004e2a2c-1acc-4873-9379-ef1aa12283b6",
                "005239a8-2e63-4ff1-9cd4-714f81837a61",
                "006b8839-31e5-4697-b912-8e3f4124dd15",
                "006ce9a8-cf38-462e-bb99-7f08499244ab",
                "007ce9b5-3268-441e-9ffd-b40d1127a319",
                "0084a614-780b-42ec-b85f-7a1b83128cd3",
                "00a5e471-a79f-4d56-8a4c-4847ac037400",
                "00ab2b5a-b59e-4ec9-b297-76f74ff1d3fb",
                "00c5f14e-a398-4076-95d1-25f320ee3a37",
                "00c74a8b-10aa-40cc-991e-3365ea1f3fce",
                "00df5a50-bce3-4edf-a078-641e54800dcb"
            ]
        }
    },
    "format":"TSV",
    "fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id",
    "size":"100"
}
```
```Shell
curl --request POST --header "Content-Type: application/json" --data @Payload.txt 'https://api.gdc.cancer.gov/files' > File_metadata.txt
```
```File_metadata_txt
cases_0_submitter_id	cases_0_case_id	data_type	cases_0_samples_0_sample_type	cases_0_samples_0_tissue_type	file_name	cases_0_samples_0_submitter_id	cases_0_samples_0_portions_0_analytes_0_aliquots_0_aliquot_id	cases_0_samples_0_sample_id	file_id	data_category	cases_0_samples_0_tumor_descriptor	cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id
TCGA-B0-5094	8aaa4e25-5c12-4ace-96dc-91aaa0c4457c	Aligned Reads	Solid Tissue Normal		C345.TCGA-B0-5094-11A-01D-1421-08.5_gdc_realn.bam	TCGA-B0-5094-11A	b4e4630a-b38c-4b62-b0e8-d73f0e3b4e47	7519d7a8-c3ee-417b-9cfc-111bc5ad0637	0001801b-54b0-4551-8d7a-d66fb59429bf	Raw Sequencing Data		TCGA-B0-5094-11A-01D-1421-08
TCGA-B0-5117	ae55b2d3-62a1-419e-9f9a-5ddfac356db4	Aligned Reads	Solid Tissue Normal		C345.TCGA-B0-5117-11A-01D-1421-08.5_gdc_realn.bam	TCGA-B0-5117-11A	45c68b6b-0bed-424d-9a77-4f87bbaa3649	b1116541-bece-4df3-b3dd-cec50aeb277b	003143c8-bbbf-46b9-a96f-f58530f4bb82	Raw Sequencing Data		TCGA-B0-5117-11A-01D-1421-08
TCGA-G7-6790	e7a1cbe2-793c-4747-8412-8be794f2382b	Aligned Reads	Blood Derived Normal		C489.TCGA-G7-6790-10A-01D-1962-08.2_gdc_realn.bam	TCGA-G7-6790-10A	66cbb40f-14b3-40c0-a332-e8a8e21bca11	4be83d0f-8b09-4e9e-8318-358371d34332	004e2a2c-1acc-4873-9379-ef1aa12283b6	Raw Sequencing Data		TCGA-G7-6790-10A-01D-1962-08
TCGA-B9-A69E	a4225cb2-7b4b-4122-b6b9-629c26e3ea56	Aligned Reads	Blood Derived Normal		TCGA-B9-A69E-10A-01D-A31X-10_Illumina_gdc_realn.bam	TCGA-B9-A69E-10A	f4799bdc-b207-4053-9a4b-5a26ebf8ab91	5d6d6cd4-6a7b-499d-936a-1be9bf74b07f	0084a614-780b-42ec-b85f-7a1b83128cd3	Raw Sequencing Data		TCGA-B9-A69E-10A-01D-A31X-10
TCGA-EE-A2GU	24faa36a-268d-4a13-b3ae-eacd431a2bcc	Aligned Reads	Blood Derived Normal		C828.TCGA-EE-A2GU-10A-01D-A198-08.2_gdc_realn.bam	TCGA-EE-A2GU-10A	c3feacc2-5a26-4bb2-a312-8b2ee53ccad1	cc4a5ed8-376a-4842-a25d-ffb07d8e1ca0	00c74a8b-10aa-40cc-991e-3365ea1f3fce	Raw Sequencing Data		TCGA-EE-A2GU-10A-01D-A198-08
TCGA-CE-A484	e62a728d-390f-428a-bea1-fc8c9814fb11	Aligned Reads	Blood Derived Normal		C499.TCGA-CE-A484-10A-01D-A23U-08.3_gdc_realn.bam	TCGA-CE-A484-10A	641a0220-6eec-434a-b606-e256113b65da	27a8008e-044a-4966-b518-cc6905e292ca	00df5a50-bce3-4edf-a078-641e54800dcb	Raw Sequencing Data		TCGA-CE-A484-10A-01D-A23U-08
TCGA-DA-A1IB	8fc9cc74-f388-49f0-b957-debb62638634	Aligned Reads	Blood Derived Normal		C828.TCGA-DA-A1IB-10A-01D-A198-08.2_gdc_realn.bam	TCGA-DA-A1IB-10A	30919a1a-df9f-4604-835e-f66ac7bcacdf	432952c5-6505-4220-a581-f65270a45281	00ab2b5a-b59e-4ec9-b297-76f74ff1d3fb	Raw Sequencing Data		TCGA-DA-A1IB-10A-01D-A198-08
TCGA-AX-A2HG	7a2cf5ce-8317-4fff-946e-b9937afab815	Aligned Reads	Blood Derived Normal		6c2a8ea343da8d6cc0fd2043492f16df_gdc_realn.bam	TCGA-AX-A2HG-10A	8c34ffe2-9012-4b4a-b610-a42a9c6a9780	ef4b80ec-b453-48ec-8ad8-ccac83e1e4db	00c5f14e-a398-4076-95d1-25f320ee3a37	Raw Sequencing Data		TCGA-AX-A2HG-10A-01D-A17D-09
TCGA-EC-A24G	b5c1e511-baf2-45b3-9919-110e8941e3c2	Aligned Reads	Blood Derived Normal		671333b193812fc2bd2744053b383459_gdc_realn.bam	TCGA-EC-A24G-10A	2a8cb8fe-b64f-453e-8139-7ede12f3fc51	61cf2e54-1b8d-40a0-9c73-a7449cbd570a	00a5e471-a79f-4d56-8a4c-4847ac037400	Raw Sequencing Data		TCGA-EC-A24G-10A-01D-A16D-09
TCGA-B5-A0K0	29c8f468-5ac1-4d6c-8376-e36e6d246926	Aligned Reads	Blood Derived Normal		TCGA-B5-A0K0-10A-01W-A062-09_IlluminaGA-DNASeq_exome_gdc_realn.bam	TCGA-B5-A0K0-10A	02e65074-ffda-4795-b8f5-1bfd20bd1019	1df69e2e-f392-465f-8e61-4671ba2fcd35	007ce9b5-3268-441e-9ffd-b40d1127a319	Raw Sequencing Data		TCGA-B5-A0K0-10A-01W-A062-09
TCGA-C8-A27B	f0d8a1fe-e313-44f1-99cc-b965cbeeff0e	Aligned Reads	Blood Derived Normal		3c99d98ea8eb6acbf819e67fc77623d9_gdc_realn.bam	TCGA-C8-A27B-10A	922226ba-6244-4953-ad42-f4daa474c288	31139082-7978-45aa-9d8f-ac4789ac5cec	006b8839-31e5-4697-b912-8e3f4124dd15	Raw Sequencing Data		TCGA-C8-A27B-10A-01D-A167-09
TCGA-E9-A295	fec0da58-1047-44d2-b6d1-c18cceed43dc	Aligned Reads	Blood Derived Normal		fd4421a6bbf3efd4e3d5c17fdd610314_gdc_realn.bam	TCGA-E9-A295-10A	cd761feb-9a20-4495-8943-c6243532a5cf	e74183e1-f0b4-412a-8dac-a62d404add78	002c67f2-ff52-4246-9d65-a3f69df6789e	Raw Sequencing Data		TCGA-E9-A295-10A-01D-A16D-09
TCGA-EB-A44O	c787c4da-c564-44f1-89eb-dd9da107acb1	Aligned Reads	Blood Derived Normal		C828.TCGA-EB-A44O-10A-01D-A25O-08.3_gdc_realn.bam	TCGA-EB-A44O-10A	c723584a-c404-4c88-bfea-e40f5dbba542	5b738547-1825-4684-81bd-864bf2eb43ef	006ce9a8-cf38-462e-bb99-7f08499244ab	Raw Sequencing Data		TCGA-EB-A44O-10A-01D-A25O-08
TCGA-A2-A3XX	53886143-c1c6-40e9-88e6-e4e5e0271fc8	Aligned Reads	Blood Derived Normal		b40998d4778f18ed80d6dd8bff0eb761_gdc_realn.bam	TCGA-A2-A3XX-10A	e96d5811-4736-40dd-966d-e0e172aeb0af	c6eb6218-ad71-40a6-88b7-a4f1a015b816	0043d981-3c6b-463f-b512-ab1d076d3e62	Raw Sequencing Data		TCGA-A2-A3XX-10A-01D-A23C-09
TCGA-EB-A3XB	a9255dcb-b236-4777-ac43-555e3a5386c3	Aligned Reads	Blood Derived Normal		C828.TCGA-EB-A3XB-10B-01D-A23B-08.1_gdc_realn.bam	TCGA-EB-A3XB-10B	9f4ffc2f-d006-4d86-b3b1-b25020481893	0e1d4c7c-204d-4765-b090-68ed4cd83835	005239a8-2e63-4ff1-9cd4-714f81837a61	Raw Sequencing Data		TCGA-EB-A3XB-10B-01D-A23B-08
```


### Format

Specifies the format of the API response: JSON (default), `TSV` or `XML`.

#### Examples

```shell1
curl  'https://api.gdc.cancer.gov/cases?fields=submitter_id&size=5&format=TSV'
```
```python1
import requests

cases_endpt = 'https://api.gdc.cancer.gov/cases'
params = {'fields':'submitter_id',
          'format':'TSV'}
response = requests.get(cases_endpt, params = params)
print response.content
```
```response1
id	submitter_id
375436b3-66ac-4d5e-b495-18a96d812a69	TCGA-F5-6810
74543fa4-ce73-46e4-9c59-224e8242b4a2	TCGA-AG-A01W
f8970455-bfb2-4b1d-ab71-3c5d619898ad	TCGA-ZN-A9VQ
c739fd61-22b2-412d-bcf3-89bda45a2c0f	TCGA-3H-AB3X
340fef21-55d8-433f-b00a-51276b849356	TCGA-MQ-A4LI

```
```shell2
curl  'https://api.gdc.cancer.gov/cases?fields=submitter_id&size=5&format=XML&pretty=true'
```
```python2
import requests

cases_endpt = 'https://api.gdc.cancer.gov/cases'
params = {'fields':'submitter_id',
          'format':'XML',
          'pretty':'true'}
response = requests.get(cases_endpt, params = params)
print response.content
```
```Output2
<?xml version="1.0" ?>
<response>
	<data>
		<hits>
			<item>
				<id>375436b3-66ac-4d5e-b495-18a96d812a69</id>
				<submitter_id>TCGA-F5-6810</submitter_id>
			</item>
			<item>
				<id>74543fa4-ce73-46e4-9c59-224e8242b4a2</id>
				<submitter_id>TCGA-AG-A01W</submitter_id>
			</item>
			<item>
				<id>f8970455-bfb2-4b1d-ab71-3c5d619898ad</id>
				<submitter_id>TCGA-ZN-A9VQ</submitter_id>
			</item>
			<item>
				<id>c739fd61-22b2-412d-bcf3-89bda45a2c0f</id>
				<submitter_id>TCGA-3H-AB3X</submitter_id>
			</item>
			<item>
				<id>340fef21-55d8-433f-b00a-51276b849356</id>
				<submitter_id>TCGA-MQ-A4LI</submitter_id>
			</item>
		</hits>
		<pagination>
			<count>5</count>
			<total>84392</total>
			<size>5</size>
			<from>0</from>
			<sort/>
			<page>1</page>
			<pages>16879</pages>
		</pagination>
	</data>
	<warnings/>
</response>
```

### Pretty

Returns when the `pretty` parameter is set to `true`, the API response is formatted with additional whitespace to improve legibility.

#### Example

```Request1
curl  'https://api.gdc.cancer.gov/cases?fields=submitter_id&sort=submitter_id:asc&size=5'
```
```Response1
{"data": {"hits": [{"id": "be37f1f7-2f98-4f74-bc04-6dd2ae2afcad", "submitter_id": "01BR001"}, {"id": "e6915db0-7c89-484d-8f9f-15cca68b82fc", "submitter_id": "01BR008"}, {"id": "16614d46-172b-479c-992b-e80a8e9a2c59", "submitter_id": "01BR009"}, {"id": "567fc9e3-17a6-42b1-a896-5e9a9507d1d8", "submitter_id": "01BR010"}, {"id": "54e89878-a1bc-4f5a-9d68-4842a469586e", "submitter_id": "01BR015"}], "pagination": {"count": 5, "total": 84392, "size": 5, "from": 0, "sort": "submitter_id:asc", "page": 1, "pages": 16879}}, "warnings": {}}
```
```Request2
curl  'https://api.gdc.cancer.gov/cases?fields=submitter_id&sort=submitter_id:asc&size=5&pretty=true'
```
```Response2
{
  "data": {
    "hits": [
      {
        "id": "be37f1f7-2f98-4f74-bc04-6dd2ae2afcad",
        "submitter_id": "01BR001"
      },
      {
        "id": "e6915db0-7c89-484d-8f9f-15cca68b82fc",
        "submitter_id": "01BR008"
      },
      {
        "id": "16614d46-172b-479c-992b-e80a8e9a2c59",
        "submitter_id": "01BR009"
      },
      {
        "id": "567fc9e3-17a6-42b1-a896-5e9a9507d1d8",
        "submitter_id": "01BR010"
      },
      {
        "id": "54e89878-a1bc-4f5a-9d68-4842a469586e",
        "submitter_id": "01BR015"
      }
    ],
    "pagination": {
      "count": 5,
      "total": 84392,
      "size": 5,
      "from": 0,
      "sort": "submitter_id:asc",
      "page": 1,
      "pages": 16879
    }
  },
  "warnings": {}
}
```

### Fields

This query parameter specifies which fields are to be included in the API response. The fields in the API response will be unordered. A listing of available fields for each endpoint is provided in [Appendix A](Appendix_A_Available_Fields.md).

#### Example

The following example requests case submitter ID, file UUID, file name and file size from the `files` endpoint.

```shell
curl 'https://api.gdc.cancer.gov/files?fields=cases.submitter_id,file_id,file_name,file_size&pretty=true'
```
```python
import requests
import json

files_endpt = 'https://api.gdc.cancer.gov/files'
params = {'fields':'cases.submitter_id,file_id,file_name,file_size'}
response = requests.get(files_endpt, params = params)
print json.dumps(response.json(), indent=2)
```
```Response
{
  "data": {
    "hits": [
      {
        "id": "c2cefea6-74d4-4859-8fe2-822767d6f68d",
        "cases": [
          {
            "submitter_id": "HCM-BROD-0003-C71"
          }
        ],
        "file_name": "30f53128-5def-4d1c-b203-9717e9cf4401_wxs_gdc_realn.bam",
        "file_id": "c2cefea6-74d4-4859-8fe2-822767d6f68d",
        "file_size": 35753708766
      },
      {
        "id": "070d2103-4350-477b-8bd7-ee529d9d24fb",
        "cases": [
          {
            "submitter_id": "HCM-CSHL-0142-C18"
          }
        ],
        "file_name": "ca5a3304-2af0-4a5f-9479-881918520921.wxs.varscan2.raw_somatic_mutation.vcf.gz",
        "file_id": "070d2103-4350-477b-8bd7-ee529d9d24fb",
        "file_size": 53177
      },
      {
        "id": "9133f158-4bea-4036-b94a-60c25385ed36",
        "cases": [
          {
            "submitter_id": "HCM-BROD-0028-C71"
          }
        ],
        "file_name": "03203e90-6180-4994-bfc8-1521669a6a49.wxs.MuSE.somatic_annotation.vcf.gz",
        "file_id": "9133f158-4bea-4036-b94a-60c25385ed36",
        "file_size": 158597
      },
      {
        "id": "22a04866-c605-4b2d-a48e-816058028c6f",
        "cases": [
          {
            "submitter_id": "HCM-BROD-0002-C71"
          }
        ],
        "file_name": "07d8e937-79e3-4fac-83c1-2e67e7a6ae14.wxs.MuSE.aliquot.maf.gz",
        "file_id": "22a04866-c605-4b2d-a48e-816058028c6f",
        "file_size": 125050
      },
      {
        "id": "2843fe16-b371-44a9-b9ab-1a93c26d24db",
        "cases": [
          {
            "submitter_id": "HCM-CSHL-0366-C50"
          }
        ],
        "file_name": "43f4ba37-92ed-4d30-86f5-e1eeb0109d9a.wgs.sanger_raw_pindel.raw_somatic_mutation.vcf.gz",
        "file_id": "2843fe16-b371-44a9-b9ab-1a93c26d24db",
        "file_size": 62812030
      },
      {
        "id": "aea93c80-0551-459e-8408-6d16148a7210",
        "cases": [
          {
            "submitter_id": "HCM-CSHL-0461-D12"
          }
        ],
        "file_name": "c45b5b42-acfb-457c-9244-2e70368c29c5.FPKM-UQ.txt.gz",
        "file_id": "aea93c80-0551-459e-8408-6d16148a7210",
        "file_size": 406727
      },
      {
        "id": "d46795b7-2166-44af-97a9-c825585878d3",
        "cases": [
          {
            "submitter_id": "HCM-CSHL-0248-C19"
          }
        ],
        "file_name": "e837ac90-f482-454f-a1e0-f16cea1d9f95.wgs.CaVEMan.raw_somatic_mutation.vcf.gz",
        "file_id": "d46795b7-2166-44af-97a9-c825585878d3",
        "file_size": 2658954
      },
      {
        "id": "7ba38c35-6491-48f7-811f-336b8487021f",
        "cases": [
          {
            "submitter_id": "HCM-CSHL-0057-C18"
          }
        ],
        "file_name": "9720b990-3ff4-4cb9-a900-ad57da72cff4.FPKM.txt.gz",
        "file_id": "7ba38c35-6491-48f7-811f-336b8487021f",
        "file_size": 358324
      },
      {
        "id": "3565c301-a03c-4334-9e06-3bb01f92c3f0",
        "cases": [
          {
            "submitter_id": "HCM-BROD-0002-C71"
          }
        ],
        "file_name": "52f32e5d-d56d-4e54-b962-9a0ed377afd3.wgs.BRASS.raw_structural_variation.bedpe.gz",
        "file_id": "3565c301-a03c-4334-9e06-3bb01f92c3f0",
        "file_size": 9715
      },
      {
        "id": "2fe7b061-48d7-45af-b435-84919ce68e47",
        "cases": [
          {
            "submitter_id": "HCM-BROD-0012-C71"
          }
        ],
        "file_name": "53984f05-821c-492a-8a0b-6e2c0b340e92.rna_seq.star_splice_junctions.tsv.gz",
        "file_id": "2fe7b061-48d7-45af-b435-84919ce68e47",
        "file_size": 2865433
      }
    ],
    "pagination": {
      "count": 10,
      "total": 596758,
      "size": 10,
      "from": 0,
      "sort": "",
      "page": 1,
      "pages": 59676
    }
  },
  "warnings": {}
}
```

### Expand

The `expand` parameter provides a shortcut to request multiple related fields (field groups) in the response. Instead of specifying each field using the `fields` parameter, users can specify a field group name using the `expand` parameter to request all fields in the group. Available field groups are listed in [Appendix A](Appendix_A_Available_Fields.md#field-group-listing-by-endpoint); the list can also be accessed programmatically at the [_mapping endpoint](#95mapping-endpoint). The `fields` and `expand` parameters can be used together to request custom combinations of field groups and individual fields.

#### Example

```Shell
curl 'https://api.gdc.cancer.gov/files/ac2ddebd-5e5e-4aea-a430-5a87c6d9c878?expand=cases.samples&pretty=true'
```
```Response
{
  "data": {
    "data_format": "BAM",
    "cases": [
      {
        "samples": [
          {
            "sample_type_id": "11",
            "tumor_descriptor": null,
            "sample_id": "b4e7558d-898e-4d68-a897-381edde0bbcc",
            "sample_type": "Solid Tissue Normal",
            "created_datetime": null,
            "tumor_code": null,
            "time_between_excision_and_freezing": null,
            "composition": null,
            "updated_datetime": "2018-11-15T21:38:54.195821-06:00",
            "days_to_collection": 5980,
            "state": "released",
            "initial_weight": 810.0,
            "preservation_method": null,
            "intermediate_dimension": null,
            "time_between_clamping_and_freezing": null,
            "freezing_method": null,
            "pathology_report_uuid": null,
            "submitter_id": "TCGA-QQ-A5VA-11A",
            "tumor_code_id": null,
            "shortest_dimension": null,
            "oct_embedded": "false",
            "days_to_sample_procurement": null,
            "longest_dimension": null,
            "current_weight": null,
            "is_ffpe": false,
            "tissue_type": "Not Reported"
          }
        ]
      }
    ],
    "access": "controlled",
    "file_name": "000aa811c15656604161e8f0e3a0aae4_gdc_realn.bam",
    "submitter_id": "32872121-d38a-4128-b96a-698a6f18f29d",
    "data_category": "Sequencing Reads",
    "acl": [
      "phs000178"
    ],
    "type": "aligned_reads",
    "platform": "Illumina",
    "created_datetime": "2016-05-26T18:55:53.506549-05:00",
    "file_size": 12667634731,
    "md5sum": "200475f5f6e42520204e5f6aadfe954f",
    "updated_datetime": "2018-11-15T21:38:44.655215-06:00",
    "file_id": "ac2ddebd-5e5e-4aea-a430-5a87c6d9c878",
    "data_type": "Aligned Reads",
    "state": "released",
    "experimental_strategy": "WXS",
    "version": "1",
    "data_release": "12.0 - 27.0"
  },
  "warnings": {}
}
```

### Size and From

GDC API provides a pagination feature that limits the number of results returned by the API. It is implemented using `size` and `from` query parameters.

The `size` query parameter specifies the maximum number of results to return. Default `size` is 10. If the number of query results is greater than `size`, only some of the results will be returned.

The `from` query parameter specifies the first record to return out of the set of results. For example, if there are 20 cases returned from the `cases` endpoint, then setting `from` to `11` will return results 12 to 20. The `from` parameter can be used in conjunction with the `size` parameter to return a specific subset of results.


#### Example


``` Shell1
curl 'https://api.gdc.cancer.gov/files?fields=file_name&from=0&size=2&pretty=true'
```
``` Python1
import requests
import json

files_endpt = 'https://api.gdc.cancer.gov/files'
params = {'fields':'file_name',
          'from':0, 'size':2}
response = requests.get(files_endpt, params = params)
print json.dumps(response.json(), indent=2)

```
```Response1
{
  "data": {
    "hits": [
      {
        "id": "c2cefea6-74d4-4859-8fe2-822767d6f68d",
        "file_name": "30f53128-5def-4d1c-b203-9717e9cf4401_wxs_gdc_realn.bam"
      },
      {
        "id": "070d2103-4350-477b-8bd7-ee529d9d24fb",
        "file_name": "ca5a3304-2af0-4a5f-9479-881918520921.wxs.varscan2.raw_somatic_mutation.vcf.gz"
      }
    ],
    "pagination": {
      "count": 2,
      "total": 596758,
      "size": 2,
      "from": 0,
      "sort": "",
      "page": 1,
      "pages": 298379
    }
  },
  "warnings": {}
}
```
``` Shell2
curl 'https://api.gdc.cancer.gov/files?fields=file_name&from=101&size=5&pretty=true'
```
``` Python2
import requests
import json

files_endpt = 'https://api.gdc.cancer.gov/files'
params = {'fields':'file_name',
          'from':101, 'size':5}
response = requests.get(files_endpt, params = params)
print json.dumps(response.json(), indent=2)
```
``` Output2
{
  "data": {
    "hits": [
      {
        "id": "79c5e6ab-7d33-48fb-8ad3-1a353f1aa8a0",
        "file_name": "bded93b7-da8b-467a-b301-bc0533780b7b.wxs.VarScan2.aliquot.maf.gz"
      },
      {
        "id": "b8b730aa-6f8f-4c7c-ad64-d69f49df56a3",
        "file_name": "5b930358-2132-4c5d-874c-9b94656dcf3b_gdc_realn.bam"
      },
      {
        "id": "bd912d6c-7325-4e90-ab1d-c3a4880f1e84",
        "file_name": "8a737363-742f-40df-ab4e-9a0bdd3adeed.wxs.varscan2.raw_somatic_mutation.vcf.gz"
      },
      {
        "id": "a11d6196-7f01-4c62-8808-0a627250c59c",
        "file_name": "4a030c3f-79e0-4737-92ad-59ba5af89977.wxs.aliquot_ensemble_raw.maf.gz"
      },
      {
        "id": "6fc778a7-6c7d-4aba-b1e4-36c2bb752216",
        "file_name": "2d9f75a7-95fd-418f-9c48-ae981b6853f1.star_fusion.rna_fusion.bedpe"
      }
    ],
    "pagination": {
      "count": 5,
      "total": 596758,
      "size": 5,
      "from": 101,
      "sort": "",
      "page": 21,
      "pages": 119352
    }
  },
  "warnings": {}
}
```

### Sort

The `sort` query parameter sorts the results by a specific field, and with the sort direction specified using the `:asc` (ascending) or `:desc` (descending) prefix, e.g. `sort=field:desc`. A list of all valid _field_ names is available in [Appendix A](Appendix_A_Available_Fields.md); the list can also be accessed programmatically at the [_mapping endpoint](#95mapping-endpoint).

#### Example

Sort cases by `submitter_id` in ascending order:

``` shell
curl  'https://api.gdc.cancer.gov/cases?fields=submitter_id&sort=submitter_id:asc&pretty=true'
```
``` python
import requests
import json

cases_endpt = 'https://api.gdc.cancer.gov/cases'
params = {'fields':'submitter_id',
          'sort':'submitter_id:asc'}
response = requests.get(cases_endpt, params = params)
print json.dumps(response.json(), indent=2)

```
``` Output
{
  "data": {
    "hits": [
      {
        "id": "be37f1f7-2f98-4f74-bc04-6dd2ae2afcad",
        "submitter_id": "01BR001"
      },
      {
        "id": "e6915db0-7c89-484d-8f9f-15cca68b82fc",
        "submitter_id": "01BR008"
      },
      {
        "id": "16614d46-172b-479c-992b-e80a8e9a2c59",
        "submitter_id": "01BR009"
      },
      {
        "id": "567fc9e3-17a6-42b1-a896-5e9a9507d1d8",
        "submitter_id": "01BR010"
      },
      {
        "id": "54e89878-a1bc-4f5a-9d68-4842a469586e",
        "submitter_id": "01BR015"
      },
      {
        "id": "a1c7b7b9-b8c8-48c3-9420-55497f9318fd",
        "submitter_id": "01BR017"
      },
      {
        "id": "ce3c8b98-e275-4cfd-a379-940d675a564b",
        "submitter_id": "01BR018"
      },
      {
        "id": "e4ce89ef-bcaa-418a-8a6b-3602793b9bbf",
        "submitter_id": "01BR020"
      },
      {
        "id": "19d3c861-8a5f-49a2-acc0-b55b25465c35",
        "submitter_id": "01BR023"
      },
      {
        "id": "afae8dce-294a-4108-bb28-376f804ae5c4",
        "submitter_id": "01BR025"
      }
    ],
    "pagination": {
      "count": 10,
      "total": 84392,
      "size": 10,
      "from": 0,
      "sort": "submitter_id:asc",
      "page": 1,
      "pages": 8440
    }
  },
  "warnings": {}
}
```

### Facets
The `facets` parameter provides aggregate information for a specified field. It provides all values that exist for that field, and the number of entities (cases, projects, files, or annotations) that this value. The primary intended use of this parameter is for displaying aggregate information in the GDC Data Portal.

The `facets` parameter can be used in conjunction with the `filters` parameter to get aggregate information for a set of search results. The following limitations apply when using `facets` and `filters` together:

1. The `filters` object's top level operator must be `and`, and the internal filters must be limited to: `=`, `!=`, `in`, `exclude`, `is`, and `not`.
2. The information provided by `facets` for a given field will disregard any filters applied to that same field.

#### Example

This is an example of a request for a count of projects in each program.

```shell
curl  'https://api.gdc.cancer.gov/projects?facets=program.name&from=0&size=0&sort=program.name:asc&pretty=true'
```
```python
import requests
import json

projects_endpt = 'https://api.gdc.cancer.gov/projects'
params = {'facets':'program.name',
          'from':0, 'size':0,
          'sort':'program.name:asc'}
response = requests.get(projects_endpt, params = params)
print json.dumps(response.json(), indent=2)
```
```Response
{
  "data": {
    "hits": [],
    "aggregations": {
      "program.name": {
        "buckets": [
          {
            "doc_count": 33,
            "key": "TCGA"
          },
          {
            "doc_count": 9,
            "key": "TARGET"
          },
          {
            "doc_count": 8,
            "key": "GENIE"
          },
          {
            "doc_count": 2,
            "key": "BEATAML1.0"
          },
          {
            "doc_count": 2,
            "key": "CGCI"
          },
          {
            "doc_count": 2,
            "key": "CMI"
          },
          {
            "doc_count": 2,
            "key": "CPTAC"
          },
          {
            "doc_count": 1,
            "key": "CTSP"
          },
          {
            "doc_count": 1,
            "key": "FM"
          },
          {
            "doc_count": 1,
            "key": "HCMI"
          },
          {
            "doc_count": 1,
            "key": "MMRF"
          },
          {
            "doc_count": 1,
            "key": "NCICCR"
          },
          {
            "doc_count": 1,
            "key": "OHSU"
          },
          {
            "doc_count": 1,
            "key": "ORGANOID"
          },
          {
            "doc_count": 1,
            "key": "VAREPOP"
          },
          {
            "doc_count": 1,
            "key": "WCDT"
          }
        ]
      }
    },
    "pagination": {
      "count": 0,
      "total": 67,
      "size": 0,
      "from": 0,
      "sort": "program.name:asc",
      "page": 1,
      "pages": 67
    }
  },
  "warnings": {}
}
```

#### Example

In this sample POST request, both `filters` and `facets` parameters are used. Note that `facets` ignores the `primary_site` filter.

```Payload
{
    "filters":{
        "op":"and",
        "content":[
            {
                "op":"=",
                "content":{
                    "field":"cases.project.primary_site",
                    "value":"Kidney"
                }
            },
            {
                "op":"=",
                "content":{
                    "field":"project.program.name",
                    "value":"TCGA"
                }
            }
        ]
    },
    "size":"0",
    "facets":"project.primary_site",
    "pretty":"true"
}
```
```Shell
curl --request POST --header "Content-Type: application/json" --data @Payload 'https://api.gdc.cancer.gov/v0/cases'
```
``` Response
{
  "data": {
    "pagination": {
      "count": 0,
      "sort": "",
      "from": 0,
      "page": 1,
      "total": 941,
      "pages": 941,
      "size": 0
    },
    "hits": [],
    "aggregations": {
      "project.primary_site": {
        "buckets": [
          {
            "key": "Brain",
            "doc_count": 1133
          },
          {
            "key": "Breast",
            "doc_count": 1098
          },
          {
            "key": "Lung",
            "doc_count": 1089
          },
          {
            "key": "Kidney",
            "doc_count": 941
          },
          {
            "key": "Colorectal",
            "doc_count": 635
          },
          {
            "key": "Uterus",
            "doc_count": 617
          },
          {
            "key": "Ovary",
            "doc_count": 608
          },
          {
            "key": "Head and Neck",
            "doc_count": 528
          },
          {
            "key": "Thyroid",
            "doc_count": 507
          },
          {
            "key": "Prostate",
            "doc_count": 500
          },
          {
            "key": "Stomach",
            "doc_count": 478
          },
          {
            "key": "Skin",
            "doc_count": 470
          },
          {
            "key": "Bladder",
            "doc_count": 412
          },
          {
            "key": "Liver",
            "doc_count": 377
          },
          {
            "key": "Cervix",
            "doc_count": 308
          },
          {
            "key": "Adrenal Gland",
            "doc_count": 271
          },
          {
            "key": "Soft Tissue",
            "doc_count": 261
          },
          {
            "key": "Bone Marrow",
            "doc_count": 200
          },
          {
            "key": "Esophagus",
            "doc_count": 185
          },
          {
            "key": "Pancreas",
            "doc_count": 185
          },
          {
            "key": "Testis",
            "doc_count": 150
          },
          {
            "key": "Thymus",
            "doc_count": 124
          },
          {
            "key": "Pleura",
            "doc_count": 87
          },
          {
            "key": "Eye",
            "doc_count": 80
          },
          {
            "key": "Lymph Nodes",
            "doc_count": 58
          },
          {
            "key": "Bile Duct",
            "doc_count": 51
          }
        ]
      }
    }
  },
  "warnings": {}
}
```


## Alternative Request Format

The GDC API also supports POST requests with `Content-Type: application/x-www-form-urlencoded` (curl default), which require payloads in the following format:
```
	filters=%7B%0A%20%20%20%20%22op%22%3A%22in%22%2C%0A%20%20%20%20%22content%22%3A%7B%0A%20%20%20%20%20%20%20%20%22field%22%3A%22files.file_id%22%2C%0A%20%20%20%20%20%20%20%20%22value%22%3A%5B%0A%20%20%20%20%20%20%20%20%20%20%20%20%220001801b-54b0-4551-8d7a-d66fb59429bf%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22002c67f2-ff52-4246-9d65-a3f69df6789e%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22003143c8-bbbf-46b9-a96f-f58530f4bb82%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%220043d981-3c6b-463f-b512-ab1d076d3e62%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22004e2a2c-1acc-4873-9379-ef1aa12283b6%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22005239a8-2e63-4ff1-9cd4-714f81837a61%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22006b8839-31e5-4697-b912-8e3f4124dd15%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22006ce9a8-cf38-462e-bb99-7f08499244ab%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22007ce9b5-3268-441e-9ffd-b40d1127a319%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%220084a614-780b-42ec-b85f-7a1b83128cd3%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%2200a5e471-a79f-4d56-8a4c-4847ac037400%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%2200ab2b5a-b59e-4ec9-b297-76f74ff1d3fb%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%2200c5f14e-a398-4076-95d1-25f320ee3a37%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%2200c74a8b-10aa-40cc-991e-3365ea1f3fce%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%2200df5a50-bce3-4edf-a078-641e54800dcb%22%0A%20%20%20%20%20%20%20%20%5D%0A%20%20%20%20%7D%0A%7D&fields=file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id&format=tsv&size=100
```
## Using Wildcards

The GDC API supports the use of the wildcard character, an asterisk (\*), in the `value` fields of a JSON query.  For example, if a user wanted to retrieve information about projects with a disease type that ended in "Adenocarcinoma" a query for `"disease_type": "*Adenocarcinoma"` would be appropriate. See below:

```
{  
   "size":"20000",
   "pretty":"TRUE",
   "fields":"submitter_id,disease_type",
   "format":"TSV",
   "filters":{  
      "op":"=",
      "content":{  
         "field":"disease_type",
         "value":"*Adenocarcinoma"
      }
   }
}
```

## Quicksearch Endpoint

The GDC Portal has a quicksearch functionality that allows for a project, case, or file to be queried from a search box. This function calls the `/v0/all` endpoint, which retrieves the top cases, files, and projects that match to the query. The quicksearch can also be used programmatically through the API.  For example, a search term of 'TCGA' would produce the following query:  

```Shell
curl "https://api.gdc.cancer.gov/v0/all?query=TCGA&size=5"
```
```Response
{"data":{"query":{"hits":[{"disease_type":["Gliomas"],"id":"UHJvamVjdDpUQ0dBLUxHRw==","name":"Brain Lower Grade Glioma","primary_site":["Brain"],"project_id":"TCGA-LGG","project_quicksearch":"Brain Lower Grade Glioma"},{"disease_type":["Myeloid Leukemias"],"id":"UHJvamVjdDpUQ0dBLUxBTUw=","name":"Acute Myeloid Leukemia","primary_site":["Hematopoietic and reticuloendothelial systems"],"project_id":"TCGA-LAML","project_quicksearch":"Acute Myeloid Leukemia"},{"disease_type":["Adenomas and Adenocarcinomas"],"id":"UHJvamVjdDpUQ0dBLUtJUkM=","name":"Kidney Renal Clear Cell Carcinoma","primary_site":["Kidney"],"project_id":"TCGA-KIRC","project_quicksearch":"Kidney Renal Clear Cell Carcinoma"},{"disease_type":["Complex Mixed and Stromal Neoplasms"],"id":"UHJvamVjdDpUQ0dBLVVDUw==","name":"Uterine Carcinosarcoma","primary_site":["Uterus, NOS"],"project_id":"TCGA-UCS","project_quicksearch":"Uterine Carcinosarcoma"},{"disease_type":["Germ Cell Neoplasms"],"id":"UHJvamVjdDpUQ0dBLVRHQ1Q=","name":"Testicular Germ Cell Tumors","primary_site":["Testis"],"project_id":"TCGA-TGCT","project_quicksearch":"Testicular Germ Cell Tumors"}],"total":195221}}}
```

This endpoint can be used to quickly retrieve information about a file.  For example, if a user wanted to know the UUID for `nationwidechildrens.org_biospecimen.TCGA-EL-A4K1.xml`, the following query could be used to quickly retrieve it programmatically:

```Shell
curl "https://api.gdc.cancer.gov/v0/all?query=nationwidechildrens.org_biospecimen.TCGA-EL-A4K1.xml&size=5"
```
```Response
{"data":{"query":{"hits":[{"file_id":"a74abfec-db78-4ed4-9e4b-604b66e30e30","file_name":"nationwidechildrens.org_biospecimen.TCGA-EL-A4K1.xml","id":"RmlsZTphNzRhYmZlYy1kYjc4LTRlZDQtOWU0Yi02MDRiNjZlMzBlMzA=","submitter_id":"nationwidechildrens.org_biospecimen.TCGA-EL-A4K1.xml"}],"total":1}}}
```

## Additional Examples

More examples of API functionality described in this section are provided in [Additional Examples](Additional_Examples.md).
