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
[{"uuid": "1dd28069-5777-4ff9-bd2b-d1ba68e88b06", "version": "1", "file_change": "superseded", "release_date": "2018-08-23", "data_release": "12.0"}, {"uuid": "76b3f4d8-c6b7-4662-ac42-1d27d4684281", "version": "2", "file_change": "released", "release_date": "2022-03-29", "data_release": "32.0"}]

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
        "id": "03974dc9-0162-4de8-9897-09f88693681a",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Nevi and Melanomas",
        "analyte_ids": [
          "9747b614-624b-410a-8b94-854a16cd143a",
          "c8974764-4836-4a34-aeb8-52b491f78d0e",
          "bcea1ed5-b9cb-4a92-ad80-598d8a223fb3"
        ],
        "submitter_id": "HCM-BROD-0334-C43",
        "submitter_analyte_ids": [
          "HCM-BROD-0334-C43-10A-01D",
          "HCM-BROD-0334-C43-85M-01D",
          "HCM-BROD-0334-C43-85M-01R"
        ],
        "days_to_consent": null,
        "aliquot_ids": [
          "dcd74e48-12f3-4a86-a829-c7e055c215b7",
          "ea182abf-041d-474a-bc53-f6fdd05cd999",
          "33f3ba0a-c902-4288-9fa7-5696d959e51d"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0334-C43-85M-01R-A79O-41",
          "HCM-BROD-0334-C43-10A-01D-A79L-36",
          "HCM-BROD-0334-C43-85M-01D-A79L-36"
        ],
        "created_datetime": "2020-05-21T08:55:40.814734-05:00",
        "diagnosis_ids": [
          "3d666f1b-58c2-451f-8ebf-87b5caa02aaf",
          "fedc3533-85f7-4fc6-b996-a1f596e021df"
        ],
        "sample_ids": [
          "cd88baf4-b6eb-4df5-9b42-d55f3aad739c",
          "eb79f8b4-1cc3-4a32-ad51-cfea8cf150f0"
        ],
        "consent_type": null,
        "submitter_sample_ids": [
          "HCM-BROD-0334-C43-10A",
          "HCM-BROD-0334-C43-85M"
        ],
        "primary_site": "Skin",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0334-C43_diagnosis2",
          "HCM-BROD-0334-C43_diagnosis"
        ],
        "updated_datetime": "2021-03-03T15:15:08.075155-06:00",
        "case_id": "03974dc9-0162-4de8-9897-09f88693681a",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "bd0bc175-5b54-47c1-96fc-c6d8afc0c115"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0334-C43-10A-01"
        ]
      },
      {
        "id": "03bfeb7c-cecf-4691-8263-33cdfe391ea9",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Adenomas and Adenocarcinomas",
        "analyte_ids": [
          "db8132c3-47a8-49ec-9b78-bc7d18debf67",
          "f74bf217-dae2-4554-92e4-8707068ea7a7",
          "01764f17-2a97-442e-a08b-8a21303b4770",
          "c9884c81-3c8f-4ad9-a962-42c7459a2276",
          "d4f1a9f8-f748-4f45-aa06-da4d760c4fab"
        ],
        "submitter_id": "HCM-BROD-0124-C25",
        "submitter_analyte_ids": [
          "HCM-BROD-0124-C25-85A-01D",
          "HCM-BROD-0124-C25-01A-01D",
          "HCM-BROD-0124-C25-10A-01D",
          "HCM-BROD-0124-C25-85A-01R",
          "HCM-BROD-0124-C25-01A-01R"
        ],
        "aliquot_ids": [
          "fe8e2565-749a-470a-b843-7afbe95ded81",
          "092656af-b279-46b8-9ccf-b1eabfbd1d6f",
          "677edb2c-fac3-4878-a28d-cf4e0d7873d7",
          "8bb9fff4-24f8-426e-9f2a-4cb30a4ac5c2",
          "f1c6f71d-b125-47bb-91ad-90d7cbff0012"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0124-C25-01A-01D-A78W-36",
          "HCM-BROD-0124-C25-85A-01D-A786-36",
          "HCM-BROD-0124-C25-01A-01R-A78X-41",
          "HCM-BROD-0124-C25-10A-01D-A78W-36",
          "HCM-BROD-0124-C25-85A-01R-A787-41"
        ],
        "created_datetime": "2019-04-16T11:21:56.471158-05:00",
        "diagnosis_ids": [
          "00184ed8-780a-4acf-b5f1-b1fcd6b08dcf"
        ],
        "sample_ids": [
          "e6bc6b9d-553f-4c78-bc11-ebcf7b0d4f27",
          "539593d1-bd9b-4379-8c86-16cdf607cd4e",
          "b0c0b5b0-cf6d-4281-8f4d-43dc77e88bc6"
        ],
        "submitter_sample_ids": [
          "HCM-BROD-0124-C25-85A",
          "HCM-BROD-0124-C25-01A",
          "HCM-BROD-0124-C25-10A"
        ],
        "primary_site": "Pancreas",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0124-C25_diagnosis"
        ],
        "updated_datetime": "2021-07-12T12:25:55.528644-05:00",
        "case_id": "03bfeb7c-cecf-4691-8263-33cdfe391ea9",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "303512c0-382b-4442-a5a7-2699ca8b1384",
          "f07a4dae-2878-452e-9836-6f39c594d38d"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0124-C25-10A-01",
          "HCM-BROD-0124-C25-01A-01"
        ]
      },
      {
        "id": "05f41641-ee22-4d41-bb87-2bfa47cd983f",
        "lost_to_followup": null,
        "slide_ids": [
          "775be57f-6df8-40c3-9c4e-c06dd900a237",
          "7bd1dea3-7819-43e5-a9e4-5fe0a189cc87"
        ],
        "submitter_slide_ids": [
          "HCM-BROD-0095-C15-06A-01-S2-HE",
          "HCM-BROD-0095-C15-06A-01-S1-HE"
        ],
        "days_to_lost_to_followup": null,
        "disease_type": "Adenomas and Adenocarcinomas",
        "analyte_ids": [
          "3cd4c8b4-8c23-4c20-bffb-97e0f5e5ac0a",
          "76ed5a0c-8129-4a73-bbd2-08d9b36bee62",
          "d2acb3a5-7b51-4d86-8c2f-18b3f886a001",
          "e6f06e9e-bf13-44fb-990c-d64b1096cd7c",
          "331fc610-2b43-4b0a-a9b2-3a8a665cb000"
        ],
        "submitter_id": "HCM-BROD-0095-C15",
        "submitter_analyte_ids": [
          "HCM-BROD-0095-C15-06A-11R",
          "HCM-BROD-0095-C15-85A-01R",
          "HCM-BROD-0095-C15-85A-01D",
          "HCM-BROD-0095-C15-10B-01D",
          "HCM-BROD-0095-C15-06A-11D"
        ],
        "days_to_consent": null,
        "aliquot_ids": [
          "24b11e4f-5e04-4d28-875e-a242515f6d07",
          "4af130bc-16ef-42c0-9f01-04dc601f4165",
          "9bd07008-27df-4e7e-be83-2d0fbeb2db94",
          "11efd698-fb14-4a52-aaec-0084afc5bbe0",
          "3f7b05f1-a8cd-48f2-89f1-4a1e4eab92c5"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0095-C15-06A-11R-A79D-41",
          "HCM-BROD-0095-C15-85A-01R-A79D-41",
          "HCM-BROD-0095-C15-10B-01D-A79C-36",
          "HCM-BROD-0095-C15-06A-11D-A79C-36",
          "HCM-BROD-0095-C15-85A-01D-A79C-36"
        ],
        "created_datetime": "2019-04-04T14:07:27.780827-05:00",
        "diagnosis_ids": [
          "cd0da3ad-189b-4be0-a8c6-55e0294e7c73"
        ],
        "sample_ids": [
          "5fcb1b3c-6711-4caa-9c6c-31ed2c0fc238",
          "71568fd7-b545-4979-a907-ef8bf41e76db",
          "d9db85d9-5ea2-41cd-96d2-d50462b5d4b6"
        ],
        "consent_type": null,
        "submitter_sample_ids": [
          "HCM-BROD-0095-C15-85A",
          "HCM-BROD-0095-C15-06A",
          "HCM-BROD-0095-C15-10B"
        ],
        "primary_site": "Esophagus",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0095-C15_diagnosis"
        ],
        "updated_datetime": "2021-01-06T22:55:10.531130-06:00",
        "case_id": "05f41641-ee22-4d41-bb87-2bfa47cd983f",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "0ef29dec-13e5-4faf-997c-a9c9502a353b",
          "c7660224-cef0-4d7d-b532-a4e9d4a7fb7c"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0095-C15-06A-11",
          "HCM-BROD-0095-C15-10B-01"
        ]
      },
      {
        "id": "07a067d0-7dfc-4817-b4c5-9200da20a59f",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Adenomas and Adenocarcinomas",
        "analyte_ids": [
          "3a6a1f63-f9d3-44f0-8659-3b3cd19a42a0",
          "6671382a-20f0-41d7-a35d-c09c7275d08a",
          "d5034c0e-2790-4fb5-82b0-e36c0626b57c"
        ],
        "submitter_id": "HCM-SANG-0268-C18",
        "submitter_analyte_ids": [
          "HCM-SANG-0268-C18-85A-01D",
          "HCM-SANG-0268-C18-85A-01R",
          "HCM-SANG-0268-C18-10A-01D"
        ],
        "days_to_consent": null,
        "aliquot_ids": [
          "e60dec0b-a2db-4b66-9e77-e8f0b8856a53",
          "2a558468-cca8-4fe9-a745-8c76436ce6cb",
          "ae8d0e5c-340b-488e-ba9e-d0912869ee8d",
          "745432c8-ea7f-48e0-bd53-28ead61a7ec0",
          "e2d32e2a-3161-40ab-88ef-0a4dfd3a72b4",
          "c2328729-0a45-4740-a9cb-3b8ff31807b3",
          "6b3085c9-63b8-4b82-8413-f775dbc6d993"
        ],
        "submitter_aliquot_ids": [
          "HCM-SANG-0268-C18-85A-01D-A80T-32-aliquot",
          "HCM-SANG-0268-C18-85A-01R-A80W-32-aliquot",
          "HCM-SANG-0268-C18-10A-01D-A79M-36",
          "HCM-SANG-0268-C18-10A-01D-A80T-32-aliquot",
          "HCM-SANG-0268-C18-85A-01R-A79O-41",
          "HCM-SANG-0268-C18-01A-01D-A80T-32-aliquot",
          "HCM-SANG-0268-C18-85A-01D-A79M-36"
        ],
        "created_datetime": "2020-05-29T12:07:57.849637-05:00",
        "diagnosis_ids": [
          "78f3dc9e-9b51-42b2-b809-405e438f7f68"
        ],
        "sample_ids": [
          "20b5f03e-26b9-4902-9bdd-667168727e6d",
          "17766175-a251-454b-9263-76af20b77290",
          "6a530ae1-79fa-4536-bbf5-86b14a80563e",
          "9eb44d2a-80c0-43d2-822a-0598a7d8e68c",
          "a3817d45-dbaa-4295-8d17-712c0c2438e4",
          "475ac2f5-327f-4d94-bebf-e03058911b59"
        ],
        "consent_type": null,
        "submitter_sample_ids": [
          "HCM-SANG-0268-C18-10A-01D-A80T-32",
          "HCM-SANG-0268-C18-85A-01D-A80T-32",
          "HCM-SANG-0268-C18-01A-01D-A80T-32",
          "HCM-SANG-0268-C18-10A",
          "HCM-SANG-0268-C18-85A-01R-A80W-32",
          "HCM-SANG-0268-C18-85A"
        ],
        "primary_site": "Colon",
        "submitter_diagnosis_ids": [
          "HCM-SANG-0268-C18_diagnosis"
        ],
        "updated_datetime": "2023-02-22T07:39:25.979291-06:00",
        "case_id": "07a067d0-7dfc-4817-b4c5-9200da20a59f",
        "index_date": "Sample Procurement",
        "state": "released",
        "portion_ids": [
          "4d2aeabd-3a5a-42fc-8bcd-301829a32883"
        ],
        "submitter_portion_ids": [
          "HCM-SANG-0268-C18-10A-01"
        ]
      },
      {
        "id": "0a6a14db-ca5c-4bf9-9125-611d672bc67b",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Adenomas and Adenocarcinomas",
        "analyte_ids": [
          "4c154bcd-3cd5-4e3a-bb90-948cf4c92965",
          "e8c8e41f-0a24-4d16-bcdd-3d7422150c4d",
          "7ce86045-0051-465e-8744-f6299a007bce",
          "f32ed211-100b-47ca-b0b0-774a15d58e60",
          "a70c036b-ea29-46d5-ba7f-8152947146f8",
          "117f4f00-2978-4a67-8a18-3acbd77fe930"
        ],
        "submitter_id": "HCM-SANG-0271-D12",
        "submitter_analyte_ids": [
          "HCM-SANG-0271-D12-86A-01R",
          "HCM-SANG-0271-D12-86A-01D",
          "HCM-SANG-0271-D12-85A-01R",
          "HCM-SANG-0271-D12-31B-01D",
          "HCM-SANG-0271-D12-10A-01D",
          "HCM-SANG-0271-D12-85B-01D"
        ],
        "days_to_consent": null,
        "aliquot_ids": [
          "d59f00c0-696c-4444-8e7f-c0140dbcd8e3",
          "d1848dd2-dccb-43c6-bd29-bb7f331089ef",
          "f7125c62-3864-412e-9afa-dc4970e43d05",
          "3bfc0b7c-58c7-47c9-9664-5c2d087f4485",
          "3109b83d-36ed-4850-9847-9710f3921413",
          "8dda870c-374e-49ae-8f0b-44b50ef567cd",
          "80d9347a-3684-4348-b1de-14bdca1471fd",
          "0bc14c0a-9fb9-480d-bfb5-3b7d89c83efd",
          "fe6747dd-12a0-4bb4-b2f3-9dd6f047df02",
          "98284402-4f76-4aa8-893c-95c72b5fb395"
        ],
        "submitter_aliquot_ids": [
          "HCM-SANG-0271-D12-85A-01R-A80W-32-aliquot",
          "HCM-SANG-0271-D12-31B-01D-A80U-36",
          "HCM-SANG-0271-D12-86A-01D-A85C-36",
          "HCM-SANG-0271-D12-85A-01D-A80T-32-aliquot",
          "HCM-SANG-0271-D12-01A-01D-A80T-32-aliquot",
          "HCM-SANG-0271-D12-85B-01D-A80U-36",
          "HCM-SANG-0271-D12-10A-01D-A80U-36",
          "HCM-SANG-0271-D12-10A-01D-A80T-32-aliquot",
          "HCM-SANG-0271-D12-85A-01R-A80V-41",
          "HCM-SANG-0271-D12-86A-01R-A85D-41"
        ],
        "created_datetime": "2020-07-08T11:54:16.081928-05:00",
        "diagnosis_ids": [
          "17a41ef7-aafc-4b70-8d8f-484dc5cd27bd"
        ],
        "sample_ids": [
          "ae2b48ea-7255-4d3d-ba00-1def687c3606",
          "b16a3a6b-b0b1-4749-afbc-aea8f4eb3a5d",
          "dc021e65-03cd-4210-a9b7-cdc971c22223",
          "bc55eaf9-56fb-4f98-9fd8-2e80796b2873",
          "5d4ac7e0-5a91-4e32-bf47-1febc9cc37d7",
          "371a0350-a942-4a63-ac59-873da0cd1e86",
          "4d8a7679-067e-4af7-9b11-d7722adc35ba",
          "2ccc3dc4-042a-477a-8b3f-b8b7ff838762",
          "fa202911-ce39-476e-9631-e81ffb46a402"
        ],
        "consent_type": null,
        "submitter_sample_ids": [
          "HCM-SANG-0271-D12-85A-01R-A80W-32",
          "HCM-SANG-0271-D12-85B",
          "HCM-SANG-0271-D12-01A-01D-A80T-32",
          "HCM-SANG-0271-D12-85A",
          "HCM-SANG-0271-D12-10A-01D-A80T-32",
          "HCM-SANG-0271-D12-10A",
          "HCM-SANG-0271-D12-31B",
          "HCM-SANG-0271-D12-86A",
          "HCM-SANG-0271-D12-85A-01D-A80T-32"
        ],
        "primary_site": "Colon",
        "submitter_diagnosis_ids": [
          "HCM-SANG-0271-D12_diagnosis"
        ],
        "updated_datetime": "2023-02-22T07:39:25.979291-06:00",
        "case_id": "0a6a14db-ca5c-4bf9-9125-611d672bc67b",
        "index_date": "Sample Procurement",
        "state": "released",
        "portion_ids": [
          "bcb558ac-e5bf-4b6c-af72-9b78b7b76722",
          "a9aec562-c531-458d-95f1-768e2610aae7"
        ],
        "submitter_portion_ids": [
          "HCM-SANG-0271-D12-10A-01",
          "HCM-SANG-0271-D12-31B-01"
        ]
      },
      {
        "id": "0cf7d1fe-e9c7-4e84-9497-df13ca2ed2c9",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Adenomas and Adenocarcinomas",
        "analyte_ids": [
          "6ef5c20b-50c3-4fda-bd1e-de4876ccab7b",
          "4a18990e-a0dc-4466-aaed-0053ffa0656a",
          "0a62151e-6c10-4e23-80ff-93970fa20f7a",
          "9b91bd41-f1f6-4712-84c9-54afa62d089b",
          "038fc540-f8e2-44b7-944c-ac6f25cac665"
        ],
        "submitter_id": "HCM-SANG-0287-C20",
        "submitter_analyte_ids": [
          "HCM-SANG-0287-C20-85A-01D",
          "HCM-SANG-0287-C20-01A-01D",
          "HCM-SANG-0287-C20-85A-01R",
          "HCM-SANG-0287-C20-10A-01D",
          "HCM-SANG-0287-C20-01A-01R"
        ],
        "aliquot_ids": [
          "43bb9903-2b44-494f-9275-62dd88803739",
          "c0b86e70-1df5-4acb-9667-cfada31b7c5e",
          "feb5f601-cd08-46bf-a8b3-4553d6ce24ba",
          "1e9d473d-c7b7-4eca-82c5-4d8b40e482a8",
          "22bc24d9-719b-4f3c-8587-e394cf409117",
          "77b7d070-0fd0-4d33-9204-b6f8cd95c415",
          "89ab0d1b-b692-4248-abc6-7dc698b54e67",
          "1870408c-6d85-4e00-9c64-178350135992",
          "48f6c991-9d58-4d1f-96e4-c62dcf0c6bf8"
        ],
        "submitter_aliquot_ids": [
          "HCM-SANG-0287-C20-01A-01D-A78U-36",
          "HCM-SANG-0287-C20-85A-01R-A80W-32-aliquot",
          "HCM-SANG-0287-C20-85A-01R-A78V-41",
          "HCM-SANG-0287-C20-10A-01D-A80T-32-aliquot",
          "HCM-SANG-0287-C20-85A-01D-A80T-32-aliquot",
          "HCM-SANG-0287-C20-01A-01D-A80T-32-aliquot",
          "HCM-SANG-0287-C20-01A-01R-A78V-41",
          "HCM-SANG-0287-C20-10A-01D-A78U-36",
          "HCM-SANG-0287-C20-85A-01D-A78U-36"
        ],
        "created_datetime": "2019-10-14T10:45:59.013881-05:00",
        "diagnosis_ids": [
          "c91c84dc-0e49-412f-9c97-8c69b195cf02"
        ],
        "sample_ids": [
          "a716e59b-6876-4c5e-85c9-a53869482d95",
          "06c2d962-ede5-4f2f-81dd-dc9253a9ddf5",
          "459eaa9f-e97f-4008-8cf6-f6671985fd30",
          "31866d02-8351-497c-8c73-1fc5fef584aa",
          "38243f28-93e3-4a8a-bf92-8019a07b3bec",
          "62ec2e41-acac-4449-b9fe-fd2c938cc811",
          "fa630ad7-f709-43ca-bf5c-d3e83dd51779"
        ],
        "submitter_sample_ids": [
          "HCM-SANG-0287-C20-01A",
          "HCM-SANG-0287-C20-85A-01D-A80T-32",
          "HCM-SANG-0287-C20-85A-01R-A80W-32",
          "HCM-SANG-0287-C20-10A",
          "HCM-SANG-0287-C20-01A-01D-A80T-32",
          "HCM-SANG-0287-C20-85A",
          "HCM-SANG-0287-C20-10A-01D-A80T-32"
        ],
        "primary_site": "Rectum",
        "submitter_diagnosis_ids": [
          "HCM-SANG-0287-C20_diagnosis"
        ],
        "updated_datetime": "2023-02-22T07:39:25.979291-06:00",
        "case_id": "0cf7d1fe-e9c7-4e84-9497-df13ca2ed2c9",
        "index_date": "Sample Procurement",
        "state": "released",
        "portion_ids": [
          "003cfb21-0ffb-44e2-9961-ebdd5c39f361",
          "a7236e0c-decf-48c3-a977-d069436420b7"
        ],
        "submitter_portion_ids": [
          "HCM-SANG-0287-C20-01A-01",
          "HCM-SANG-0287-C20-10A-01"
        ]
      },
      {
        "id": "0e9a9e97-f0bf-4f4a-84cc-73eccfc627b1",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Adenomas and Adenocarcinomas",
        "analyte_ids": [
          "78c3017c-6de1-4ccb-bf4d-27922b0b1f38",
          "f877946f-8a20-40dd-b5e0-5d4a350c2528",
          "c70fdba0-7324-4ea1-b15c-8de62a459a50",
          "9c7ccff1-5150-46a4-94dd-a928df9db3e7",
          "9bc12019-74f2-4784-9322-0398a2c1b3d1"
        ],
        "submitter_id": "HCM-CSHL-0376-D37",
        "submitter_analyte_ids": [
          "HCM-CSHL-0376-D37-31A-11D",
          "HCM-CSHL-0376-D37-10A-01D",
          "HCM-CSHL-0376-D37-31A-11R",
          "HCM-CSHL-0376-D37-85P-01D",
          "HCM-CSHL-0376-D37-85P-01R"
        ],
        "aliquot_ids": [
          "790b5310-35ac-4382-b062-db3078e8b20a",
          "33019b27-f718-4436-95fd-ff4d8a28c923",
          "6c69fc6a-6606-4f79-8ad8-4bf6bfaab56d",
          "aa16d6e6-6a3d-487f-82b1-9ffa258b851e",
          "5facfb50-5793-43a6-be50-5c1afe446dcc"
        ],
        "submitter_aliquot_ids": [
          "HCM-CSHL-0376-D37-31A-11R-A78V-41",
          "HCM-CSHL-0376-D37-85P-01D-A78T-36",
          "HCM-CSHL-0376-D37-31A-11D-A78T-36",
          "HCM-CSHL-0376-D37-85P-01R-A78V-41",
          "HCM-CSHL-0376-D37-10A-01D-A78T-36"
        ],
        "created_datetime": "2019-10-14T13:24:10.078043-05:00",
        "diagnosis_ids": [
          "9ce37f22-bbde-4447-bbf0-28e85f4f2837"
        ],
        "sample_ids": [
          "ac11b2da-38a9-444e-a922-3b78f15be942",
          "cbfc2432-c011-41ee-81fd-5efd8c0dac79",
          "76b129b9-954b-4471-949b-e118dc778e2d"
        ],
        "submitter_sample_ids": [
          "HCM-CSHL-0376-D37-10A",
          "HCM-CSHL-0376-D37-85P",
          "HCM-CSHL-0376-D37-31A"
        ],
        "primary_site": "Colon",
        "submitter_diagnosis_ids": [
          "HCM-CSHL-0376-D37_diagnosis"
        ],
        "updated_datetime": "2023-02-22T07:39:25.979291-06:00",
        "case_id": "0e9a9e97-f0bf-4f4a-84cc-73eccfc627b1",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "d0c90c2f-d6d9-4601-b2f8-6ca30b48f405",
          "f2b394b2-35f1-492b-969e-df863a2714cc"
        ],
        "submitter_portion_ids": [
          "HCM-CSHL-0376-D37-10A-01",
          "HCM-CSHL-0376-D37-31A-11"
        ]
      },
      {
        "id": "149a8565-e0c5-4474-a693-d44f1b445c0c",
        "lost_to_followup": "Yes",
        "slide_ids": [
          "846c933d-4bca-4995-ac45-1caf38ee481b",
          "4be18133-d1ff-400f-b57b-dded404249c0",
          "8f24d92e-cfee-4cce-8eec-e4f302766255",
          "08650e03-96e9-4213-8892-b3da6928171f"
        ],
        "submitter_slide_ids": [
          "HCM-BROD-0199-C71-01A-01-S2-HE",
          "HCM-BROD-0199-C71-02A-01-S2-HE",
          "HCM-BROD-0199-C71-01A-01-S1-HE",
          "HCM-BROD-0199-C71-02A-01-S1-HE"
        ],
        "days_to_lost_to_followup": null,
        "disease_type": "Gliomas",
        "analyte_ids": [
          "3b8a2f07-a6d0-4998-8fe8-c2138372e191",
          "ae4eb402-7c13-46f0-bdf8-361c7fdc9430",
          "b720a33e-c858-47f4-9e51-e58748004e96",
          "fdf06863-7d14-4f79-ad2c-b46d443d235c",
          "a4f88ff3-7c9e-4edc-808f-aa01378ec68d",
          "cca66f05-89e0-48a7-8948-53e556fee5be",
          "e658c087-c37d-47ac-9472-3b6b6eed7188",
          "6be6aff1-9694-4c7f-8f43-f0f902b95849",
          "39f9ce44-944c-49b8-b215-05b117d5e62d"
        ],
        "submitter_id": "HCM-BROD-0199-C71",
        "submitter_analyte_ids": [
          "HCM-BROD-0199-C71-85A-01D",
          "HCM-BROD-0199-C71-01A-11D",
          "HCM-BROD-0199-C71-85S-01R",
          "HCM-BROD-0199-C71-85A-01R",
          "HCM-BROD-0199-C71-02A-11R",
          "HCM-BROD-0199-C71-02A-11D",
          "HCM-BROD-0199-C71-01A-11R",
          "HCM-BROD-0199-C71-10A-01D",
          "HCM-BROD-0199-C71-85R-01D"
        ],
        "aliquot_ids": [
          "9f18daaa-cd67-436e-85a1-1e9b3e1e2135",
          "eff144c4-4e2a-497d-998a-5ae6bccf2576",
          "7ec23cac-61f4-4f8f-afa7-ed8a3ca1493f",
          "6ea15978-902f-45ba-b148-cc4247341882",
          "529f3cbf-f65d-4a60-a562-b3fcfbd7d4c9",
          "22542d68-476b-47fd-91da-03db887756d6",
          "7f4e7f6b-33b0-49da-95bd-88643d5e14ff",
          "fc0121b9-5de8-4af6-90db-b36dd8207ebf",
          "c9dd10db-bf23-4bb2-a6db-fec1b45549c8"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0199-C71-85R-01D-A80U-36",
          "HCM-BROD-0199-C71-02A-11R-A80V-41",
          "HCM-BROD-0199-C71-02A-11D-A80U-36",
          "HCM-BROD-0199-C71-01A-11D-A786-36",
          "HCM-BROD-0199-C71-10A-01D-A786-36",
          "HCM-BROD-0199-C71-01A-11R-A787-41",
          "HCM-BROD-0199-C71-85A-01D-A786-36",
          "HCM-BROD-0199-C71-85S-01R-A80V-41",
          "HCM-BROD-0199-C71-85A-01R-A787-41"
        ],
        "created_datetime": "2019-04-04T15:00:32.807421-05:00",
        "diagnosis_ids": [
          "df0230e6-d07c-4814-9aee-5be560d1ce58",
          "2d72ecab-4038-4ef7-b921-79b88ad62722"
        ],
        "sample_ids": [
          "b5919dd1-039e-4ab4-b6d5-37b9d483893f",
          "769c5e1b-89e3-431c-9420-31de5efe5a22",
          "5e0faf77-f1b3-4b1a-8b6a-8066270eddb4",
          "ee0e94ad-b5b0-4ad6-912c-c4010f1a1d26",
          "600a5b42-f871-47a2-8124-68184edf10bd",
          "a1dbbc0c-173c-4455-a780-8682dd2e258a"
        ],
        "submitter_sample_ids": [
          "HCM-BROD-0199-C71-85R",
          "HCM-BROD-0199-C71-02A",
          "HCM-BROD-0199-C71-85S",
          "HCM-BROD-0199-C71-85A",
          "HCM-BROD-0199-C71-10A",
          "HCM-BROD-0199-C71-01A"
        ],
        "primary_site": "Brain",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0199-C71_diagnosis",
          "HCM-BROD-0199-C71_diagnosis2"
        ],
        "updated_datetime": "2021-01-06T22:55:10.531130-06:00",
        "case_id": "149a8565-e0c5-4474-a693-d44f1b445c0c",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "9dccf7ac-58a0-4c9e-9fb6-4ead8c37d48b",
          "5afaa056-b09e-4328-ae5b-feb70faa3595",
          "77827fa6-18a9-4e1a-b49f-800cec351fa8"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0199-C71-10A-01",
          "HCM-BROD-0199-C71-02A-11",
          "HCM-BROD-0199-C71-01A-11"
        ]
      },
      {
        "id": "19b1e69a-355a-4dd7-9c56-d701f6c2c5a0",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Adenomas and Adenocarcinomas",
        "analyte_ids": [
          "181d5e6e-026d-4983-97e8-f4d9e28a0cfe",
          "293683aa-f2ef-4d14-8ec7-681c870d3b71",
          "c1f6530c-004c-4b10-bae5-a571542aabd2",
          "2140e379-dd8b-4440-a626-bb0e01f8fc00",
          "4d337e44-6468-43b8-bf2f-4301750dab99"
        ],
        "submitter_id": "HCM-SANG-0299-C15",
        "submitter_analyte_ids": [
          "HCM-SANG-0299-C15-10B-01D",
          "HCM-SANG-0299-C15-85X-01R",
          "HCM-SANG-0299-C15-85A-01R",
          "HCM-SANG-0299-C15-85X-01D",
          "HCM-SANG-0299-C15-85A-01D"
        ],
        "aliquot_ids": [
          "37fee31c-8669-4057-ae95-424586fa2a05",
          "d17fd7af-fd14-4cb7-a6c1-c23339589288",
          "65b7ffa3-0f1f-4918-a75f-2343720fe40c",
          "a331201c-3aef-4eec-83d5-d38f4211c1b1",
          "3ef40a4d-21cd-4be4-8c57-cba30c8e0778",
          "b8af2a87-5b55-4c73-88f1-90ce5e1d05f5",
          "a1988273-c9f8-4cd1-b1a6-5daa2c7a3e51",
          "b535f87b-0a61-49bd-828b-bb5e3c11f2f7",
          "3bb67aa8-5155-425e-9afc-9aca5120b0f0",
          "fa20b306-a835-46d8-ab03-f8ea0b585381",
          "65ec4e30-e782-4b6f-b985-33e9c7e72a0a"
        ],
        "submitter_aliquot_ids": [
          "HCM-SANG-0299-C15-85B-01D-A80T-32-aliquot",
          "HCM-SANG-0299-C15-10A-01D-A80T-32-aliquot",
          "HCM-SANG-0299-C15-85A-01R-A78V-41",
          "HCM-SANG-0299-C15-85A-01D-A78U-36",
          "HCM-SANG-0299-C15-85A-01R-A80W-32-aliquot",
          "HCM-SANG-0299-C15-85X-01D-A78U-36",
          "HCM-SANG-0299-C15-85X-01R-A78V-41",
          "HCM-SANG-0299-C15-85A-01D-A80T-32-aliquot",
          "HCM-SANG-0299-C15-01A-01D-A80T-32-aliquot",
          "HCM-SANG-0299-C15-85B-01R-A80W-32-aliquot",
          "HCM-SANG-0299-C15-10B-01D-A78U-36"
        ],
        "created_datetime": "2019-10-14T10:46:36.257369-05:00",
        "diagnosis_ids": [
          "7c6aa4ce-6661-4491-a827-c0a8045743a6",
          "c1a6f70f-b871-4b4e-a292-aef67c2d4776"
        ],
        "sample_ids": [
          "c56f3f94-deda-4d21-8f8e-658108995dfa",
          "ce46b9f4-5244-432d-aa87-026e0a27d71a",
          "cf134dfb-5126-4bee-bcab-39d584335a21",
          "34e0d6ee-97f5-420f-9b34-4784098125f7",
          "ac8476c1-53bf-43b5-8695-db441ed1a720",
          "9000b366-a14b-44bc-a782-484e09765b2d",
          "8541bc01-57e6-4bf7-a42e-da1c8e790633",
          "277d7ed9-4e6a-429b-8f4e-438ff0d2ba7a",
          "d12d9184-5317-4730-af2e-a8428456a2a7"
        ],
        "submitter_sample_ids": [
          "HCM-SANG-0299-C15-85X",
          "HCM-SANG-0299-C15-85B-01D-A80T-32",
          "HCM-SANG-0299-C15-85A-01R-A80W-32",
          "HCM-SANG-0299-C15-10B",
          "HCM-SANG-0299-C15-01A-01D-A80T-32",
          "HCM-SANG-0299-C15-85A",
          "HCM-SANG-0299-C15-85A-01D-A80T-32",
          "HCM-SANG-0299-C15-85B-01R-A80W-32",
          "HCM-SANG-0299-C15-10A-01D-A80T-32"
        ],
        "primary_site": "Esophagus",
        "submitter_diagnosis_ids": [
          "HCM-SANG-0299-C15_diagnosis2",
          "HCM-SANG-0299-C15_diagnosis"
        ],
        "updated_datetime": "2023-02-22T07:39:25.979291-06:00",
        "case_id": "19b1e69a-355a-4dd7-9c56-d701f6c2c5a0",
        "index_date": "Sample Procurement",
        "state": "released",
        "portion_ids": [
          "b8192c57-9cda-4dca-a590-1e13beadf2a0"
        ],
        "submitter_portion_ids": [
          "HCM-SANG-0299-C15-10B-01"
        ]
      },
      {
        "id": "19f1d344-4c14-4733-abbd-c2db6737e210",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Ductal and Lobular Neoplasms",
        "analyte_ids": [
          "a21d5f05-af57-41bf-ab60-32d2c869611b",
          "696f4b03-5f0d-4eaa-975d-a9feb64dae07",
          "c8d4fb23-c55d-4bb6-bd3d-fc1f159d7a33",
          "5db09eb8-fecc-433f-aa34-a12a7c9333dd",
          "51598aba-be8e-42a0-bf6c-aca34776fc1f"
        ],
        "submitter_id": "HCM-CSHL-0081-C25",
        "submitter_analyte_ids": [
          "HCM-CSHL-0081-C25-11A-11D",
          "HCM-CSHL-0081-C25-85A-01R",
          "HCM-CSHL-0081-C25-85B-01D",
          "HCM-CSHL-0081-C25-01A-11R",
          "HCM-CSHL-0081-C25-01A-11D"
        ],
        "aliquot_ids": [
          "0a65e1bd-6af3-44ba-924c-193eb8e099d6",
          "61112c39-c838-4f84-89c7-a33bfe5dea88",
          "25cd08c3-09a9-406d-8ec3-ab4946224cf1",
          "973d4fa3-4fa7-4ff9-8176-cc51c53b7079",
          "9aa21acf-1e24-4b7a-a3c1-73354bdd81b6"
        ],
        "submitter_aliquot_ids": [
          "HCM-CSHL-0081-C25-11A-11D-A78M-36",
          "HCM-CSHL-0081-C25-01A-11R-A78N-41",
          "HCM-CSHL-0081-C25-85B-01D-A78M-36",
          "HCM-CSHL-0081-C25-01A-11D-A78M-36",
          "HCM-CSHL-0081-C25-85A-01R-A78N-41"
        ],
        "created_datetime": "2019-09-19T08:58:31.776805-05:00",
        "diagnosis_ids": [
          "c1cca3a7-e0ac-40a1-9db7-6902b48d3c62"
        ],
        "sample_ids": [
          "d9f23187-9a29-426a-9ead-4bb3a2ce6cf9",
          "3709004e-b04d-4473-aa29-8dd84176d17d",
          "adcc54e3-074b-4ca6-b179-0a5df8efeb36",
          "05478d15-885a-4c44-a46a-81bbe6c9ee11"
        ],
        "submitter_sample_ids": [
          "HCM-CSHL-0081-C25-85B",
          "HCM-CSHL-0081-C25-01A",
          "HCM-CSHL-0081-C25-85A",
          "HCM-CSHL-0081-C25-11A"
        ],
        "primary_site": "Pancreas",
        "submitter_diagnosis_ids": [
          "HCM-CSHL-0081-C25_diagnosis"
        ],
        "updated_datetime": "2023-02-22T07:39:25.979291-06:00",
        "case_id": "19f1d344-4c14-4733-abbd-c2db6737e210",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "d9588542-d7ef-413e-900c-3f816b583525",
          "b77573f1-e6a3-43c1-a56d-a207c39e18c4"
        ],
        "submitter_portion_ids": [
          "HCM-CSHL-0081-C25-01A-11",
          "HCM-CSHL-0081-C25-11A-11"
        ]
      }
    ],
    "pagination": {
      "count": 10,
      "total": 40232,
      "size": 10,
      "from": 0,
      "sort": "",
      "page": 1,
      "pages": 4024
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
0286c31b-a704-4d7d-99e3-0bc4e8975b8b	HCM-CSHL-0084-C25
02f6d684-b6b5-419a-b0e1-b74d0a384a30	HCM-BROD-0408-C71
03974dc9-0162-4de8-9897-09f88693681a	HCM-BROD-0334-C43
03bfeb7c-cecf-4691-8263-33cdfe391ea9	HCM-BROD-0124-C25
04cbceab-f945-482b-956b-840756a17a4a	HCM-BROD-0421-C71

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
				<id>0286c31b-a704-4d7d-99e3-0bc4e8975b8b</id>
				<submitter_id>HCM-CSHL-0084-C25</submitter_id>
			</item>
			<item>
				<id>02f6d684-b6b5-419a-b0e1-b74d0a384a30</id>
				<submitter_id>HCM-BROD-0408-C71</submitter_id>
			</item>
			<item>
				<id>03974dc9-0162-4de8-9897-09f88693681a</id>
				<submitter_id>HCM-BROD-0334-C43</submitter_id>
			</item>
			<item>
				<id>03bfeb7c-cecf-4691-8263-33cdfe391ea9</id>
				<submitter_id>HCM-BROD-0124-C25</submitter_id>
			</item>
			<item>
				<id>04cbceab-f945-482b-956b-840756a17a4a</id>
				<submitter_id>HCM-BROD-0421-C71</submitter_id>
			</item>
		</hits>
		<pagination>
			<count>5</count>
			<total>86962</total>
			<size>5</size>
			<from>0</from>
			<sort/>
			<page>1</page>
			<pages>17393</pages>
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
{"data": {"hits": [{"id": "be37f1f7-2f98-4f74-bc04-6dd2ae2afcad", "submitter_id": "01BR001"}, {"id": "e6915db0-7c89-484d-8f9f-15cca68b82fc", "submitter_id": "01BR008"}, {"id": "16614d46-172b-479c-992b-e80a8e9a2c59", "submitter_id": "01BR009"}, {"id": "567fc9e3-17a6-42b1-a896-5e9a9507d1d8", "submitter_id": "01BR010"}, {"id": "54e89878-a1bc-4f5a-9d68-4842a469586e", "submitter_id": "01BR015"}], "pagination": {"count": 5, "total": 86962, "size": 5, "from": 0, "sort": "None", "page": 1, "pages": 17393}}, "warnings": {}}
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
      "total": 86962,
      "size": 5,
      "from": 0,
      "sort": "None",
      "page": 1,
      "pages": 17393
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
        "id": "d570eccc-3c1c-4c4f-ae04-96be71fbe016",
        "cases": [
          {
            "submitter_id": "TCGA-AN-A0FL"
          }
        ],
        "file_name": "TCGA-AN-A0FL-01Z-00-DX1.20A041C6-A306-4599-A7D1-65032A252AA9.svs",
        "file_id": "d570eccc-3c1c-4c4f-ae04-96be71fbe016",
        "file_size": 1055798681
      },
      {
        "id": "0f8d8202-a1ca-4ea1-98b2-c20a6b08479a",
        "cases": [
          {
            "submitter_id": "TCGA-AN-A0FL"
          }
        ],
        "file_name": "nationwidechildrens.org_ssf.TCGA-AN-A0FL.xml",
        "file_id": "0f8d8202-a1ca-4ea1-98b2-c20a6b08479a",
        "file_size": 15519
      },
      {
        "id": "b76f87b3-99c5-4297-b2df-8cbea8ecaf61",
        "cases": [
          {
            "submitter_id": "TCGA-BH-A18F"
          }
        ],
        "file_name": "7c4e4c2a-a0b1-424f-97d8-359825674429.wxs.aliquot_ensemble_masked.maf.gz",
        "file_id": "b76f87b3-99c5-4297-b2df-8cbea8ecaf61",
        "file_size": 21571
      },
      {
        "id": "be6d269d-4305-4643-b98e-af703a067761",
        "cases": [
          {
            "submitter_id": "TCGA-BH-A18F"
          }
        ],
        "file_name": "HITCH_p_TCGASNP_b93_N_GenomeWideSNP_6_E11_741424.CEL",
        "file_id": "be6d269d-4305-4643-b98e-af703a067761",
        "file_size": 69084893
      },
      {
        "id": "fed73119-1d5e-4f7e-9713-183d1916422b",
        "cases": [
          {
            "submitter_id": "TCGA-BH-A18F"
          }
        ],
        "file_name": "3b928f83-14a7-4bd6-a9b0-744b414d4495.wxs.varscan2.raw_somatic_mutation.vcf.gz",
        "file_id": "fed73119-1d5e-4f7e-9713-183d1916422b",
        "file_size": 35903
      },
      {
        "id": "6877b045-91f1-4030-82ff-b90507e11e17",
        "cases": [
          {
            "submitter_id": "TCGA-BH-A18F"
          }
        ],
        "file_name": "5057e3cb-25cd-4a67-8d31-6ac8508ba3c7.methylation_array.sesame.level3betas.txt",
        "file_id": "6877b045-91f1-4030-82ff-b90507e11e17",
        "file_size": 770500
      },
      {
        "id": "07e8cdc7-d228-4752-ad19-800abd507277",
        "cases": [
          {
            "submitter_id": "TCGA-BH-A0BM"
          }
        ],
        "file_name": "TCGA-BRCA.28dcad29-448e-4bcb-911d-556c6f4a5573.star_fusion.rna_fusion.tsv",
        "file_id": "07e8cdc7-d228-4752-ad19-800abd507277",
        "file_size": 234
      },
      {
        "id": "fef57b45-ede1-49b0-b60d-957a55a15e0e",
        "cases": [
          {
            "submitter_id": "TCGA-BH-A0BM"
          }
        ],
        "file_name": "nationwidechildrens.org_biospecimen.TCGA-BH-A0BM.xml",
        "file_id": "fef57b45-ede1-49b0-b60d-957a55a15e0e",
        "file_size": 127218
      },
      {
        "id": "81a1b323-88b6-4837-bccf-ac84a79828b6",
        "cases": [
          {
            "submitter_id": "TCGA-BH-A0BM"
          }
        ],
        "file_name": "TCGA-BRCA.4570b87f-8116-48bf-86d3-b993536c88db.gene_level_copy_number.v36.tsv",
        "file_id": "81a1b323-88b6-4837-bccf-ac84a79828b6",
        "file_size": 3446816
      },
      {
        "id": "c6bf94a6-9940-4155-86b4-bbb10875dbdb",
        "cases": [
          {
            "submitter_id": "TCGA-BH-A18F"
          }
        ],
        "file_name": "TCGA-BRCA.88cae21a-4890-4fdd-a678-c4864620942c.star_fusion.rna_fusion.bedpe",
        "file_id": "c6bf94a6-9940-4155-86b4-bbb10875dbdb",
        "file_size": 229
      }
    ],
    "pagination": {
      "count": 10,
      "total": 931947,
      "size": 10,
      "from": 0,
      "sort": "",
      "page": 1,
      "pages": 93195
    }
  },
  "warnings": {}
}
```

### Expand

The `expand` parameter provides a shortcut to request multiple related fields (field groups) in the response. Instead of specifying each field using the `fields` parameter, users can specify a field group name using the `expand` parameter to request all fields in the group. Available field groups are listed in [Appendix A](Appendix_A_Available_Fields.md#field-group-listing-by-endpoint); the list can also be accessed programmatically at the [_mapping endpoint](#95mapping-endpoint). The `fields` and `expand` parameters can be used together to request custom combinations of field groups and individual fields.

#### Example

```Shell
curl 'https://api.gdc.cancer.gov/files/573ee7e9-b8bd-419e-808b-a027c4311731?expand=cases.samples&pretty=true'
```
```Response
{
  "data": {
    "proportion_reads_mapped": 0.9648433596149857,
    "access": "controlled",
    "proportion_base_mismatch": 0.004117986,
    "contamination_error": 0,
    "acl": [
      "phs000178"
    ],
    "type": "aligned_reads",
    "platform": "Illumina",
    "created_datetime": "2022-05-12T14:42:10.014925-05:00",
    "md5sum": "25e89d79c47d12b13f4eb8e8f39f7f64",
    "updated_datetime": "2022-11-01T11:52:54.136033-05:00",
    "pairs_on_diff_chr": 1170013,
    "state": "released",
    "data_format": "BAM",
    "total_reads": 379313036,
    "proportion_coverage_30x": 0.000109,
    "cases": [
      {
        "samples": [
          {
            "sample_type_id": "10",
            "tumor_descriptor": "Not Reported",
            "sample_id": "4e128a37-be58-477a-a01f-448179360b7c",
            "sample_type": "Blood Derived Normal",
            "tumor_code": null,
            "created_datetime": null,
            "time_between_excision_and_freezing": null,
            "composition": "Not Reported",
            "updated_datetime": "2022-04-28T22:05:09.013808-05:00",
            "days_to_collection": 6755,
            "state": "released",
            "initial_weight": null,
            "preservation_method": null,
            "intermediate_dimension": null,
            "time_between_clamping_and_freezing": null,
            "freezing_method": null,
            "pathology_report_uuid": null,
            "submitter_id": "TCGA-B6-A0RI-10A",
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
    "file_name": "c9478f7d-bfe3-4e80-8161-39b3d440fa16_wgs_gdc_realn.bam",
    "mean_coverage": 5.452655,
    "proportion_reads_duplicated": 0.009253781617987946,
    "submitter_id": "a4e380e5-420e-49af-986d-e721601065fb",
    "data_category": "Sequencing Reads",
    "proportion_coverage_10x": 0.07674,
    "file_size": 42958286722,
    "contamination": 0,
    "average_base_quality": 32,
    "file_id": "573ee7e9-b8bd-419e-808b-a027c4311731",
    "data_type": "Aligned Reads",
    "average_insert_size": 207,
    "average_read_length": 51,
    "experimental_strategy": "WGS",
    "version": "1",
    "data_release": "36.0 - 37.0"
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
        "id": "d570eccc-3c1c-4c4f-ae04-96be71fbe016",
        "file_name": "TCGA-AN-A0FL-01Z-00-DX1.20A041C6-A306-4599-A7D1-65032A252AA9.svs"
      },
      {
        "id": "0f8d8202-a1ca-4ea1-98b2-c20a6b08479a",
        "file_name": "nationwidechildrens.org_ssf.TCGA-AN-A0FL.xml"
      }
    ],
    "pagination": {
      "count": 2,
      "total": 931947,
      "size": 2,
      "from": 0,
      "sort": "",
      "page": 1,
      "pages": 465974
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
        "id": "297933f5-1316-4cb6-b53f-9dbfa7f3d7ed",
        "file_name": "TCGA-B6-A0RH-01A-02-TSB.ea83f31e-defb-4436-8a58-5b66b18d13b5.svs"
      },
      {
        "id": "2f31e897-b3e8-49f1-a400-ccf9f00f294a",
        "file_name": "URAEI_p_TCGASNP_b85_N_GenomeWideSNP_6_F01_735050.grch38.seg.v2.txt"
      },
      {
        "id": "ebd6cf90-4f6b-4193-887a-22fdb5645fbc",
        "file_name": "TCGA-BRCA.5994c06d-ee9b-4ead-b3d1-2e1f286f7d6d.ascat2.allelic_specific.seg.txt"
      },
      {
        "id": "aebd6b5a-e676-4357-93df-523b31b55ea0",
        "file_name": "TCGA-BRCA.c737131c-636f-4e1b-89b8-bb2d6ddd8164.star_fusion.rna_fusion.bedpe"
      },
      {
        "id": "aa83a7e7-e9cc-4330-a7be-ca750cffb74c",
        "file_name": "URAEI_p_TCGASNP_b85_N_GenomeWideSNP_6_F01_735050.birdseed.data.txt"
      }
    ],
    "pagination": {
      "count": 5,
      "total": 931947,
      "size": 5,
      "from": 101,
      "sort": "",
      "page": 21,
      "pages": 186390
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
      "total": 86962,
      "size": 10,
      "from": 0,
      "sort": "None",
      "page": 1,
      "pages": 8697
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
            "doc_count": 3,
            "key": "CGCI"
          },
          {
            "doc_count": 3,
            "key": "CMI"
          },
          {
            "doc_count": 3,
            "key": "MATCH"
          },
          {
            "doc_count": 2,
            "key": "BEATAML1.0"
          },
          {
            "doc_count": 2,
            "key": "CPTAC"
          },
          {
            "doc_count": 1,
            "key": "APOLLO"
          },
          {
            "doc_count": 1,
            "key": "CDDP_EAGLE"
          },
          {
            "doc_count": 1,
            "key": "CTSP"
          },
          {
            "doc_count": 1,
            "key": "EXCEPTIONAL_RESPONDERS"
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
            "key": "MP2PRT"
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
            "key": "REBC"
          },
          {
            "doc_count": 1,
            "key": "TRIO"
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
      "total": 78,
      "size": 0,
      "from": 0,
      "sort": "None",
      "page": 1,
      "pages": 78
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
    "hits": [],
    "aggregations": {
      "project.primary_site": {
        "buckets": [
          {
            "doc_count": 1202,
            "key": "kidney"
          },
          {
            "doc_count": 1191,
            "key": "brain"
          },
          {
            "doc_count": 1176,
            "key": "bronchus and lung"
          },
          {
            "doc_count": 1156,
            "key": "breast"
          },
          {
            "doc_count": 952,
            "key": "colon"
          },
          {
            "doc_count": 947,
            "key": "stomach"
          },
          {
            "doc_count": 878,
            "key": "uterus, nos"
          },
          {
            "doc_count": 869,
            "key": "ovary"
          },
          {
            "doc_count": 821,
            "key": "corpus uteri"
          },
          {
            "doc_count": 789,
            "key": "other and unspecified parts of tongue"
          },
          {
            "doc_count": 670,
            "key": "connective, subcutaneous and other soft tissues"
          },
          {
            "doc_count": 633,
            "key": "rectosigmoid junction"
          },
          {
            "doc_count": 586,
            "key": "bones, joints and articular cartilage of other and unspecified sites"
          },
          {
            "doc_count": 565,
            "key": "thyroid gland"
          },
          {
            "doc_count": 528,
            "key": "base of tongue"
          },
          {
            "doc_count": 528,
            "key": "floor of mouth"
          },
          {
            "doc_count": 528,
            "key": "gum"
          },
          {
            "doc_count": 528,
            "key": "hypopharynx"
          },
          {
            "doc_count": 528,
            "key": "larynx"
          },
          {
            "doc_count": 528,
            "key": "lip"
          },
          {
            "doc_count": 528,
            "key": "oropharynx"
          },
          {
            "doc_count": 528,
            "key": "other and ill-defined sites in lip, oral cavity and pharynx"
          },
          {
            "doc_count": 528,
            "key": "other and unspecified parts of mouth"
          },
          {
            "doc_count": 528,
            "key": "palate"
          },
          {
            "doc_count": 528,
            "key": "tonsil"
          },
          {
            "doc_count": 500,
            "key": "prostate gland"
          },
          {
            "doc_count": 498,
            "key": "retroperitoneum and peritoneum"
          },
          {
            "doc_count": 470,
            "key": "skin"
          },
          {
            "doc_count": 448,
            "key": "heart, mediastinum, and pleura"
          },
          {
            "doc_count": 428,
            "key": "liver and intrahepatic bile ducts"
          },
          {
            "doc_count": 412,
            "key": "bladder"
          },
          {
            "doc_count": 307,
            "key": "cervix uteri"
          },
          {
            "doc_count": 271,
            "key": "adrenal gland"
          },
          {
            "doc_count": 261,
            "key": "bones, joints and articular cartilage of limbs"
          },
          {
            "doc_count": 261,
            "key": "meninges"
          },
          {
            "doc_count": 261,
            "key": "other and unspecified male genital organs"
          },
          {
            "doc_count": 261,
            "key": "peripheral nerves and autonomic nervous system"
          },
          {
            "doc_count": 258,
            "key": "hematopoietic and reticuloendothelial systems"
          },
          {
            "doc_count": 208,
            "key": "testis"
          },
          {
            "doc_count": 185,
            "key": "esophagus"
          },
          {
            "doc_count": 185,
            "key": "pancreas"
          },
          {
            "doc_count": 179,
            "key": "other and ill-defined sites"
          },
          {
            "doc_count": 179,
            "key": "other endocrine glands and related structures"
          },
          {
            "doc_count": 179,
            "key": "spinal cord, cranial nerves, and other parts of central nervous system"
          },
          {
            "doc_count": 172,
            "key": "rectum"
          },
          {
            "doc_count": 172,
            "key": "unknown"
          },
          {
            "doc_count": 124,
            "key": "thymus"
          },
          {
            "doc_count": 80,
            "key": "eye and adnexa"
          },
          {
            "doc_count": 58,
            "key": "lymph nodes"
          },
          {
            "doc_count": 58,
            "key": "other and unspecified major salivary glands"
          },
          {
            "doc_count": 58,
            "key": "small intestine"
          },
          {
            "doc_count": 51,
            "key": "gallbladder"
          },
          {
            "doc_count": 51,
            "key": "other and unspecified parts of biliary tract"
          }
        ]
      }
    },
    "pagination": {
      "count": 0,
      "total": 1202,
      "size": 0,
      "from": 0,
      "sort": "",
      "page": 1,
      "pages": 1202
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
{"data":{"query":{"hits":[{"disease_type":["Adenomas and Adenocarcinomas"],"id":"UHJvamVjdDpUQ0dBLUFDQw==","name":"Adrenocortical Carcinoma","primary_site":["Adrenal gland"],"project_id":"TCGA-ACC","project_quicksearch":"Adrenocortical Carcinoma"},{"disease_type":["Adenomas and Adenocarcinomas"],"id":"UHJvamVjdDpUQ0dBLUtJQ0g=","name":"Kidney Chromophobe","primary_site":["Kidney"],"project_id":"TCGA-KICH","project_quicksearch":"Kidney Chromophobe"},{"disease_type":["Adenomas and Adenocarcinomas"],"id":"UHJvamVjdDpUQ0dBLUxJSEM=","name":"Liver Hepatocellular Carcinoma","primary_site":["Liver and intrahepatic bile ducts"],"project_id":"TCGA-LIHC","project_quicksearch":"Liver Hepatocellular Carcinoma"},{"disease_type":["Myeloid Leukemias"],"id":"UHJvamVjdDpUQ0dBLUxBTUw=","name":"Acute Myeloid Leukemia","primary_site":["Hematopoietic and reticuloendothelial systems"],"project_id":"TCGA-LAML","project_quicksearch":"Acute Myeloid Leukemia"},{"disease_type":["Adenomas and Adenocarcinomas"],"id":"UHJvamVjdDpUQ0dBLUtJUlA=","name":"Kidney Renal Papillary Cell Carcinoma","primary_site":["Kidney"],"project_id":"TCGA-KIRP","project_quicksearch":"Kidney Renal Papillary Cell Carcinoma"}],"total":183550}}}
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
