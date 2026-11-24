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
| [\_mapping](/API/Users_Guide/Search_and_Retrieval/#_mapping-endpoint) | Information about elements that can be used to query other endpoints |

The choice of endpoint determines what is listed in the search results. The `files` endpoint will generate a list of files, whereas the `cases` endpoint will generate a list of cases. Each of the above endpoints, other than `_mapping`, can query and return any of the related fields in the [GDC Data Model](../../Data/Data_Model/GDC_Data_Model.md). So the `cases` endpoint can be queried for file fields (e.g. to look for cases that have certain types of experimental data), and the `files` endpoint can be queried for clinical metadata associated with a case (e.g. to look for files from cases diagnosed with a specific cancer type).

### `Project` Endpoint
The `projects` endpoint provides access to project records, the highest level of data organization in the GDC.

#### Example
This example is a query for projects contained in the GDC. It uses the [from](#size-and-from), [size](#size-and-from), [sort](#sort), and [pretty](#pretty) parameters, and returns the first two projects sorted by project id.

=== "Shell"

    ```shell
    curl 'https://api.gdc.cancer.gov/projects?from=0&size=2&sort=project_id:asc&pretty=true'
    ```

=== "Output"

    ```json
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

=== "Shell"

    ```shell
    curl 'https://api.gdc.cancer.gov/projects/TARGET-NBL?expand=summary,summary.experimental_strategies,summary.data_categories&pretty=true'
    ```

=== "Response"    

    ```json
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

This example is a query for files contained in the GDC. It uses the [from](#size-and-from), [size](#size-and-from), [sort](#sort), and [pretty](#pretty) parameters, and returns only the first two files, sorted by file size, from smallest to largest.

=== "Shell"

    ```shell
    curl 'https://api.gdc.cancer.gov/files?from=0&size=2&sort=file_size:asc&pretty=true'
    ```
=== "Output"

    ```json 
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

=== "Shell"

    ```Shell
    curl 'https://api.gdc.cancer.gov/files/20f45e04-3c10-4f11-b57b-719880eab69e?pretty=true'
    ```

=== "Output"

    ```json
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

=== "Shell1"

    ```shell
    curl 'https://api.gdc.cancer.gov/files/versions/1dd28069-5777-4ff9-bd2b-d1ba68e88b06,2a03abac-f1a2-49a9-a57c-7543739dd862?pretty=true'
    ```

=== "Output1"

    ```json 
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

=== "Shell2"

    ```shell
    curl --request POST --header "Content-Type: text/tsv"  https://api.gdc.cancer.gov/files/versions/manifest?pretty=true --data-binary @gdc_manifest_20180809_154816.txt
    ```
    
=== "Output2"

    ```json
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

This example is a query for files contained in GDC. It returns case where submitter id is `TCGA-BH-A0EA`, using the [pretty](#pretty) and [filters](#filters-specifying-the-query) parameters and the following [filtering operators](#query-format):

	{"op":"and","content":[{"op":"in","content":{"field":"submitter_id","value":["TCGA-BH-A0EA"]}}]}

Command:

=== "Shell"

    ```shell
    curl 'https://api.gdc.cancer.gov/cases?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22submitter_id%22%2C%22value%22%3A%5B%22TCGA-BH-A0EA%22%5D%7D%7D%5D%7D%0A%0A&pretty=true'
    ```

=== "Output"

    ```json 
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

=== "Shell"

    ```shell
    curl 'https://api.gdc.cancer.gov/cases/1f601832-eee3-48fb-acf5-80c4a454f26e?pretty=true&expand=diagnoses'
    ```

=== "Response"

    ```json
    
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
            "tumor_grade": "Not Reported",
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

The query uses the [filters](#filters-specifying-the-query) parameter to specify entity UUIDs. Code samples below include the bare and percent-encoded filter JSON.

=== "Filter-json"

    ```json
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

=== "Filter-JSON-percent-encoded"

    ```
    %7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22entity_id%22%2C%22value%22%3A%5B%22e0d36cc0-652c-4224-bb10-09d15c7bd8f1%22%2C%2225ebc29a-7598-4ae4-ba7f-618d448882cc%22%2C%22fe660d7c-2746-4b50-ab93-b2ed99960553%22%5D%7D%7D
    ```

=== "Shell"

    ```shell
    curl 'https://api.gdc.cancer.gov/annotations?filters=%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22entity_id%22%2C%22value%22%3A%5B%22e0d36cc0-652c-4224-bb10-09d15c7bd8f1%22%2C%2225ebc29a-7598-4ae4-ba7f-618d448882cc%22%2C%22fe660d7c-2746-4b50-ab93-b2ed99960553%22%5D%7D%7D&pretty=true'
    ```

=== "Output"

    ```json 
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

=== "Shell"

    ```shell
    curl 'https://api.gdc.cancer.gov/history/1dd28069-5777-4ff9-bd2b-d1ba68e88b06'
    ```
=== "Output"

    ```json 
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

=== "Shell"

    ```shell
    curl 'https://api.gdc.cancer.gov/projects/_mapping'
    ```

=== "Output"

    ```json
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
| =        | equals (string or number)                        | one                | sex_at_birth = "female"                                      |
| !=       | does not equal (string or number)                | one                | project_id != "TARGET-AML"                                   |
| <        | less than (number)                               | one                | age at diagnosis < 90y                                       |
| <=       | less than or equal (number)                      | one                | age at diagnosis <= 17                                       |
| >        | greater than (number)                            | one                | age at diagnosis > 50                                        |
| >=       | greater than or equal (number)                   | one                | age at diagnosis >= 18                                       |
| is       | is (missing)                                     | one                | sex_at_birth is missing                                      |
| not      | not (missing)                                    | one                | race not missing                                             |
| in       | matches a string or number in (a list)           | multiple           | primary_site in [Brain, Lung]                                |
| exclude  | removes results that only match the query value in the specified field | multiple           | experimental_strategy exclude [WXS, WGS, "Genotyping array"] |
|excludeifany| removes results that match at least one instance of the query value in the specified field | multiple           | experimental_strategy excludeifany [WXS, WGS, "Genotyping array"] |
| and      | (operation1) and (operation2)                    | multiple           | (primary_site in [Brain, Lung]) and (sex_at_birth = "female") |
| or       | (operation1) or (operation2)                     | multiple           | (project_id != "TARGET-AML") or (age at diagnosis < 90y)     |

The `field` operand specifies a field that corresponds to a property defined in the [GDC Data Dictionary](../../Data_Dictionary/viewer.md). A list of supported fields is provided in [Appendix A](Appendix_A_Available_Fields.md); the list can also be accessed programmatically at the [_mapping endpoint](#_mapping-endpoint).

The `value` operand specifies the search terms. Users can get a list of available values for a specific property by making a call to the appropriate API endpoint using the `facets` parameter, e.g. `https://api.gdc.cancer.gov/v0/cases?facets=demographic.sex_at_birth&size=0&pretty=true`. See [Facets](#facets) for details.

A simple query with a single operator looks like this:

	{
	    "op":"=",
	    "content":{
	        "field":"cases.demographic.sex_at_birth",
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

#### Example: `exclude` and `excludeifany` Operators

This example demonstrates the behavior of the `exclude` and `excludeifany` operators, showcasing how each operator behaves when cases contain multiple values for the `diagnoses.classification_of_tumor` field.

Note: When a field contains only a single value for a given result (such as `file_name` for a file), the `exclude` and `excludeifany` operators produce identical results.

##### No `exclude` or `excludeifany` Filters Example

Querying a single `case_id` where at least one diagnosis contains `classification_of_tumor` equal to `metastasis` may return multiple diagnoses within the case. In this example, the returned case includes two diagnoses: one with `classification_of_tumor` set to `metastasis` and another set to `primary`.

=== "Shell"

    ```shell
    curl --location 'https://api.gdc.cancer.gov/cases' \
    --header 'Content-Type: application/json' \
    --data '{
        "size": 10,
        "filters": {
            "op": "and",
            "content": [
                {
                    "op": "=",
                    "content": {
                        "field": "case_id",
                        "value": [
                            "33165ed4-5f38-4566-8d0d-82e2a6d2743f"
                        ]
                    }
                },            
                {
                    "op": "=",
                    "content": {
                        "field": "diagnoses.classification_of_tumor",
                        "value": [
                            "metastasis"
                        ]
                    }
                }
            ]
        },
        "expand": "diagnoses",
        "fields": "diagnoses.classification_of_tumor"
    } '
    ```

=== "Output"

    ```json 
    {
        "data": {
            "hits": [
                {
                    "id": "33165ed4-5f38-4566-8d0d-82e2a6d2743f",
                    "diagnoses": [
                        {
                            "morphology": "Not Reported",
                            "submitter_id": "HCM-BROD-0871-C49_diagnosis2",
                            "days_to_diagnosis": 3041,
                            "created_datetime": "2021-12-21T13:56:21.205280-06:00",
                            "tissue_or_organ_of_origin": "Uterus, NOS",
                            "age_at_diagnosis": 25208,
                            "primary_diagnosis": "Not Reported",
                            "classification_of_tumor": "metastasis",
                            "updated_datetime": "2025-04-07T17:13:09.336794-05:00",
                            "diagnosis_id": "d245f226-ed78-44d4-a452-9f219612c88c",
                            "icd_10_code": "C79.8",
                            "site_of_resection_or_biopsy": "Pelvis, NOS",
                            "state": "released",
                            "diagnosis_is_primary_disease": false
                        },
                        {
                            "ajcc_pathologic_t": "T0",
                            "morphology": "8890/3",
                            "ajcc_pathologic_stage": "Stage III",
                            "ajcc_pathologic_n": "N0",
                            "ajcc_pathologic_m": "M0",
                            "submitter_id": "HCM-BROD-0871-C49_diagnosis",
                            "created_datetime": "2021-12-21T13:56:21.205280-06:00",
                            "tissue_or_organ_of_origin": "Uterus, NOS",
                            "days_to_last_follow_up": 3065.0,
                            "age_at_diagnosis": 22167,
                            "primary_diagnosis": "Leiomyosarcoma, NOS",
                            "ajcc_clinical_stage": "Stage III",
                            "classification_of_tumor": "primary",
                            "updated_datetime": "2025-04-07T17:13:09.336794-05:00",
                            "metastasis_at_diagnosis": "No Metastasis",
                            "diagnosis_id": "f297b2b5-29a7-4f2d-a094-e2666fb04bf5",
                            "prior_malignancy": "yes",
                            "icd_10_code": "C49.9",
                            "site_of_resection_or_biopsy": "Not Reported",
                            "state": "released",
                            "prior_treatment": "No",
                            "diagnosis_is_primary_disease": true,
                            "ajcc_staging_system_edition": "8th"
                        }
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

##### Exclude Example

In this example, the query uses the `exclude` operator, which excludes cases only if *all* diagnoses within a case have a `classification_of_tumor` value matching the specified exclusion. The returned result includes the same case, as it contains both a diagnosis with `primary` and one with `metastasis`. This demonstrates that `exclude` will not remove the case unless every diagnosis matches the exclusion criteria.

=== "Shell"

    ```shell
    curl --location 'https://api.gdc.cancer.gov/cases' \
    --header 'Content-Type: application/json' \
    --data '{
        "size": 10,
        "filters": {
            "op": "and",
            "content": [
                {
                    "op": "=",
                    "content": {
                        "field": "case_id",
                        "value": [
                            "33165ed4-5f38-4566-8d0d-82e2a6d2743f"
                        ]
                    }
                },            
                {
                    "op": "exclude",
                    "content": {
                        "field": "diagnoses.classification_of_tumor",
                        "value": [
                            "metastasis"
                        ]
                    }
                }
            ]
        },
        "expand": "diagnoses",
        "fields": "diagnoses.classification_of_tumor"
    } '
    ```

=== "Output"

    ```json 
    {
        "data": {
            "hits": [
                {
                    "id": "33165ed4-5f38-4566-8d0d-82e2a6d2743f",
                    "diagnoses": [
                        {
                            "morphology": "Not Reported",
                            "submitter_id": "HCM-BROD-0871-C49_diagnosis2",
                            "days_to_diagnosis": 3041,
                            "created_datetime": "2021-12-21T13:56:21.205280-06:00",
                            "tissue_or_organ_of_origin": "Uterus, NOS",
                            "age_at_diagnosis": 25208,
                            "primary_diagnosis": "Not Reported",
                            "classification_of_tumor": "metastasis",
                            "updated_datetime": "2025-04-07T17:13:09.336794-05:00",
                            "diagnosis_id": "d245f226-ed78-44d4-a452-9f219612c88c",
                            "icd_10_code": "C79.8",
                            "site_of_resection_or_biopsy": "Pelvis, NOS",
                            "state": "released",
                            "diagnosis_is_primary_disease": false
                        },
                        {
                            "ajcc_pathologic_t": "T0",
                            "morphology": "8890/3",
                            "ajcc_pathologic_stage": "Stage III",
                            "ajcc_pathologic_n": "N0",
                            "ajcc_pathologic_m": "M0",
                            "submitter_id": "HCM-BROD-0871-C49_diagnosis",
                            "created_datetime": "2021-12-21T13:56:21.205280-06:00",
                            "tissue_or_organ_of_origin": "Uterus, NOS",
                            "days_to_last_follow_up": 3065.0,
                            "age_at_diagnosis": 22167,
                            "primary_diagnosis": "Leiomyosarcoma, NOS",
                            "ajcc_clinical_stage": "Stage III",
                            "classification_of_tumor": "primary",
                            "updated_datetime": "2025-04-07T17:13:09.336794-05:00",
                            "metastasis_at_diagnosis": "No Metastasis",
                            "diagnosis_id": "f297b2b5-29a7-4f2d-a094-e2666fb04bf5",
                            "prior_malignancy": "yes",
                            "icd_10_code": "C49.9",
                            "site_of_resection_or_biopsy": "Not Reported",
                            "state": "released",
                            "prior_treatment": "No",
                            "diagnosis_is_primary_disease": true,
                            "ajcc_staging_system_edition": "8th"
                        }
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

##### Excludeifany Example

When the query is updated to use the `excludeifany` operator, no results are returned because at least one diagnosis within the case has a `classification_of_tumor` value of `metastasis`. The `excludeifany` operator excludes any case where any element in the array matches the exclusion value.

=== "Shell"

    ```shell
    curl --location 'https://api.gdc.cancer.gov/cases' \
    --header 'Content-Type: application/json' \
    --data '{
        "size": 10,
        "filters": {
            "op": "and",
            "content": [
                {
                    "op": "=",
                    "content": {
                        "field": "case_id",
                        "value": [
                            "33165ed4-5f38-4566-8d0d-82e2a6d2743f"
                        ]
                    }
                },            
                {
                    "op": "excludeifany",
                    "content": {
                        "field": "diagnoses.classification_of_tumor",
                        "value": [
                            "metastasis"
                        ]
                    }
                }
            ]
        },
        "expand": "diagnoses",
        "fields": "diagnoses.classification_of_tumor"
    } ' 
    ```

=== "Output"

    ```json 
    {
        "data": {
            "hits": [],
            "pagination": {
                "count": 0,
                "total": 0,
                "size": 10,
                "from": 0,
                "sort": "",
                "page": 0,
                "pages": 0
            }
        },
        "warnings": {}
    }
    ```


#### Example: HTTP GET Request

This example requests `male` cases using HTTP GET.

The JSON object to be passed to the GDC API looks like:

    {
        "op": "=",
        "content": {
            "field": "cases.demographic.sex_at_birth",
            "value": [
                "male"
           ]
        }
    }

URL-encoding the above JSON object using [Percent-(URL)-encoding tool](https://www.freeformatter.com/url-encoder.html) results in the following string:

    %7B%0D%0A++++%22op%22%3A+%22%3D%22%2C%0D%0A++++%22content%22%3A+%7B%0D%0A++++++++%22field%22%3A+%22cases.demographic.sex_at_birth%22%2C%0D%0A++++++++%22value%22%3A+%5B%0D%0A++++++++++++%22male%22%0D%0A+++++++%5D%0D%0A++++%7D%0D%0A%7D

The above string can now be passed to the GDC API using the `filters` parameter:

=== "Shell"

    ```shell
    curl 'https://api.gdc.cancer.gov/cases?filters=%7B%0D%0A++++%22op%22%3A+%22%3D%22%2C%0D%0A++++%22content%22%3A+%7B%0D%0A++++++++%22field%22%3A+%22cases.demographic.sex_at_birth%22%2C%0D%0A++++++++%22value%22%3A+%5B%0D%0A++++++++++++%22male%22%0D%0A+++++++%5D%0D%0A++++%7D%0D%0A%7D&pretty=true'
    ```

=== "Python"

    ```python
    import requests
    import json
    cases_endpt = 'https://api.gdc.cancer.gov/cases'
    filt = {"op":"=",
            "content":{
                "field": "cases.demographic.sex_at_birth",
                "value": ["male"]
            }
    }
    params = {'filters':json.dumps(filt), 'sort':'demographic.sex_at_birth:asc'}
    # requests URL-encodes automatically
    response = requests.get(cases_endpt, params = params)
    print(json.dumps(response.json(), indent=2))
    ```

=== "Output"

    ```json 
    {
      "data": {
        "hits": [
          {
            "id": "0e9262d1-5aa8-4528-9aee-2815afcd23cd",
            "slide_ids": [
              "3e754cbb-ce5e-41d2-aef6-51a96b39b8ea",
              "5f4b2950-9427-4bf4-96db-239da84e5224"
            ],
            "submitter_slide_ids": [
              "HCM-WCMC-0950-C67-01A-01-S2-HE",
              "HCM-WCMC-0950-C67-01A-01-S1-HE"
            ],
            "disease_type": "Transitional Cell Papillomas and Carcinomas",
            "analyte_ids": [
              "f4c8bb09-0593-47e1-b134-69992529e9bb",
              "f153b565-462c-4649-9f77-42aea65b35a6",
              "c777e68c-9a38-4d6b-ac88-d357c06a8796",
              "e83734f5-607b-43e9-b5f7-c0003891d547",
              "e664a9f4-1e80-48cf-af47-6b4a9b612f25"
            ],
            "submitter_id": "HCM-WCMC-0950-C67",
            "submitter_analyte_ids": [
              "HCM-WCMC-0950-C67-85A-01R",
              "HCM-WCMC-0950-C67-10A-01D",
              "HCM-WCMC-0950-C67-85A-01D",
              "HCM-WCMC-0950-C67-01A-01R",
              "HCM-WCMC-0950-C67-01A-01D"
            ],
            "aliquot_ids": [
              "91911835-c196-4df8-84e4-41151e72f571",
              "318d0505-345e-40cd-b6cf-fcf2413bd828",
              "7d9ee243-b34f-4b77-81af-a75ee45a8558",
              "fa2bc724-43cd-4a01-ab7d-77f2d16783e7",
              "a31597ac-b537-46e4-8d7d-5cfee78eadd6"
            ],
            "submitter_aliquot_ids": [
              "HCM-WCMC-0950-C67-01A-01R-A88K-41",
              "HCM-WCMC-0950-C67-85A-01D-A88H-36",
              "HCM-WCMC-0950-C67-10A-01D-A88H-36",
              "HCM-WCMC-0950-C67-85A-01R-A88J-41",
              "HCM-WCMC-0950-C67-01A-01D-A88H-36"
            ],
            "created_datetime": "2021-12-09T14:30:39.067286-06:00",
            "diagnosis_ids": [
              "82806bbe-6770-4697-93aa-c8b87aa73044"
            ],
            "sample_ids": [
              "00380420-c546-4d7b-ae1a-4cd90b4f729e",
              "40b7f970-85c1-41fc-ba25-f247806646f7",
              "a80e45ae-b144-40dc-9ffc-e8b24052c277"
            ],
            "submitter_sample_ids": [
              "HCM-WCMC-0950-C67-10A",
              "HCM-WCMC-0950-C67-85A",
              "HCM-WCMC-0950-C67-01A"
            ],
            "primary_site": "Bladder",
            "submitter_diagnosis_ids": [
              "HCM-WCMC-0950-C67_diagnosis"
            ],
            "updated_datetime": "2025-04-07T17:13:09.336794-05:00",
            "case_id": "0e9262d1-5aa8-4528-9aee-2815afcd23cd",
            "index_date": "Diagnosis",
            "state": "released",
            "portion_ids": [
              "a82ae730-818d-4e24-a3b5-04506073043e",
              "f3d0cd1f-a7c3-49ce-81e2-931bda3f63c3"
            ],
            "submitter_portion_ids": [
              "HCM-WCMC-0950-C67-01A-01",
              "HCM-WCMC-0950-C67-10A-01"
            ]
          },
          {
            "id": "0ee53efd-b992-4aeb-a091-8dd1bd32da6e",
            "disease_type": "Adenomas and Adenocarcinomas",
            "analyte_ids": [
              "c2a69b6f-6e97-4765-8f05-2ef5b37be4df",
              "440b278b-03c2-4a16-bdc3-c78677046166",
              "960b4dc1-061d-4459-b59d-be0d06f4674c"
            ],
            "submitter_id": "HCM-BROD-0682-C64",
            "submitter_analyte_ids": [
              "HCM-BROD-0682-C64-85A-01D",
              "HCM-BROD-0682-C64-85A-01R",
              "HCM-BROD-0682-C64-10A-01D"
            ],
            "aliquot_ids": [
              "7ff361d1-a60e-477d-b6cb-cf1a5c61894d",
              "ae0ab148-5d9a-4e61-a043-e2a84bce1d7d",
              "6e93d944-bf0c-46ce-b558-151c96da4e42"
            ],
            "submitter_aliquot_ids": [
              "HCM-BROD-0682-C64-85A-01D-A85C-36",
              "HCM-BROD-0682-C64-10A-01D-A85C-36",
              "HCM-BROD-0682-C64-85A-01R-A85D-41"
            ],
            "created_datetime": "2021-06-03T11:50:05.387795-05:00",
            "diagnosis_ids": [
              "3c95e3a3-c388-452c-88e8-13704668f55a"
            ],
            "sample_ids": [
              "76590128-9b0f-4aa2-b7b9-62acb3aa9cb0",
              "333bca08-ddef-4d31-8eaf-0999514c88d1"
            ],
            "submitter_sample_ids": [
              "HCM-BROD-0682-C64-10A",
              "HCM-BROD-0682-C64-85A"
            ],
            "primary_site": "Kidney",
            "submitter_diagnosis_ids": [
              "HCM-BROD-0682-C64_diagnosis"
            ],
            "updated_datetime": "2025-04-07T17:13:09.336794-05:00",
            "case_id": "0ee53efd-b992-4aeb-a091-8dd1bd32da6e",
            "index_date": "Diagnosis",
            "state": "released",
            "portion_ids": [
              "a3f89e93-a6ed-4dcb-a380-4b422ffcc6d7"
            ],
            "submitter_portion_ids": [
              "HCM-BROD-0682-C64-10A-01"
            ]
          },
          {
            "id": "5868f125-ca3b-4b98-85ea-64f0ea5b45da",
            "slide_ids": [
              "a2f8961c-02ae-40d4-8701-e31692fcc096",
              "87173226-109b-4f62-9689-2c987bec00dc"
            ],
            "submitter_slide_ids": [
              "HCM-CSHL-0582-C18-06A-01-S1-HE",
              "HCM-CSHL-0582-C18-06A-01-S2-HE"
            ],
            "disease_type": "Adenomas and Adenocarcinomas",
            "analyte_ids": [
              "cac05fae-8042-4d82-a516-aa8cabbdc735",
              "cf74ccd3-dae0-4b01-b29f-fc3c0a364ca5",
              "dc4396e0-f380-4c99-bae7-4cefed5eeea6",
              "465cb1b3-fd6c-4f7e-86b1-e0dd0843c321",
              "8e0c3202-15af-4fbd-ad93-e0801df03000"
            ],
            "submitter_id": "HCM-CSHL-0582-C18",
            "submitter_analyte_ids": [
              "HCM-CSHL-0582-C18-06A-11R",
              "HCM-CSHL-0582-C18-85M-01D",
              "HCM-CSHL-0582-C18-06A-11D",
              "HCM-CSHL-0582-C18-10A-01D",
              "HCM-CSHL-0582-C18-85M-01R"
            ],
            "aliquot_ids": [
              "12053654-c327-4a01-9e6c-46409abd5f20",
              "9b662aa0-89a5-413a-9601-3c6cba2406d5",
              "bc113dd7-84e8-4b7b-912e-564e0ad959a9",
              "61c389ec-9771-489c-9bf0-401dca5da2cd",
              "f253ed43-a6e0-430b-a2fa-a108f5a06749"
            ],
            "submitter_aliquot_ids": [
              "HCM-CSHL-0582-C18-06A-11D-A82H-36",
              "HCM-CSHL-0582-C18-06A-11R-A82J-41",
              "HCM-CSHL-0582-C18-10A-01D-A82H-36",
              "HCM-CSHL-0582-C18-85M-01R-A82J-41",
              "HCM-CSHL-0582-C18-85M-01D-A82H-36"
            ],
            "created_datetime": "2020-10-30T18:29:27.459298-05:00",
            "diagnosis_ids": [
              "a3a2c172-34bc-4879-98ce-b520afe288d9",
              "21126e95-4518-417c-bd7c-e7c00c88a0e7"
            ],
            "sample_ids": [
              "adfc9442-eefc-4639-9b03-452ed10b3ef4",
              "ce06a937-9e13-4b43-a793-915bb64274e3",
              "ff88a67e-f0da-4848-9c25-fd5f42c6ea1e"
            ],
            "submitter_sample_ids": [
              "HCM-CSHL-0582-C18-85M",
              "HCM-CSHL-0582-C18-06A",
              "HCM-CSHL-0582-C18-10A"
            ],
            "primary_site": "Colon",
            "submitter_diagnosis_ids": [
              "HCM-CSHL-0582-C18_diagnosis2",
              "HCM-CSHL-0582-C18_diagnosis"
            ],
            "updated_datetime": "2025-03-31T16:34:59.906497-05:00",
            "case_id": "5868f125-ca3b-4b98-85ea-64f0ea5b45da",
            "index_date": "Diagnosis",
            "state": "released",
            "portion_ids": [
              "cbb65334-6871-4b04-a7b3-0e952d44a299",
              "a15dddff-e4e6-49e2-a5b7-2234fa0ffab4"
            ],
            "submitter_portion_ids": [
              "HCM-CSHL-0582-C18-06A-11",
              "HCM-CSHL-0582-C18-10A-01"
            ]
          },
          {
            "id": "dd7c58fc-0765-4d02-9a60-515a721cd369",
            "lost_to_followup": null,
            "days_to_lost_to_followup": null,
            "disease_type": "Adenomas and Adenocarcinomas",
            "analyte_ids": [
              "9a62892d-21a9-415c-b790-4c28805afeb1",
              "0908a123-760d-4de7-bbfc-5623fd593aa4",
              "51943fc2-631a-452d-a95d-80e204824578"
            ],
            "submitter_id": "HCM-BROD-0328-C15",
            "submitter_analyte_ids": [
              "HCM-BROD-0328-C15-85A-01D",
              "HCM-BROD-0328-C15-85A-01R",
              "HCM-BROD-0328-C15-10B-01D"
            ],
            "aliquot_ids": [
              "5e7a0ec0-5d47-4080-83b8-ece859036acb",
              "67285f61-6429-47ae-a385-2848c75fbe18",
              "286c9b98-534c-4f5c-839d-1a3e2c74f5a6"
            ],
            "submitter_aliquot_ids": [
              "HCM-BROD-0328-C15-85A-01D-A79L-36",
              "HCM-BROD-0328-C15-10B-01D-A79L-36",
              "HCM-BROD-0328-C15-85A-01R-A79O-41"
            ],
            "created_datetime": "2019-04-04T15:56:24.085531-05:00",
            "diagnosis_ids": [
              "adb35634-5e40-4885-a3ac-9b33f83e2e90"
            ],
            "sample_ids": [
              "e9957167-8692-48fe-9c24-43fc57f5bc55",
              "fc8cc23b-a948-45e6-a518-6109c78ce276"
            ],
            "submitter_sample_ids": [
              "HCM-BROD-0328-C15-10B",
              "HCM-BROD-0328-C15-85A"
            ],
            "primary_site": "Esophagus",
            "submitter_diagnosis_ids": [
              "HCM-BROD-0328-C15_diagnosis"
            ],
            "updated_datetime": "2025-03-31T16:34:59.906497-05:00",
            "case_id": "dd7c58fc-0765-4d02-9a60-515a721cd369",
            "index_date": "Diagnosis",
            "state": "released",
            "portion_ids": [
              "2d2a22c9-2438-4614-85f3-74898a697f69"
            ],
            "submitter_portion_ids": [
              "HCM-BROD-0328-C15-10B-01"
            ]
          },
          {
            "id": "12603388-e5f0-45ac-ab7b-3a839de7a9e8",
            "disease_type": "Adenomas and Adenocarcinomas",
            "submitter_id": "HCM-SANG-1335-C18",
            "aliquot_ids": [
              "28d6f03e-eb51-4146-b846-b25f7a113c53",
              "50066afd-9d47-427b-a36c-569b1fd56000",
              "5496d096-d085-426c-af20-02f22bacfc8b"
            ],
            "submitter_aliquot_ids": [
              "HCM-SANG-1335-C18-85A-01D-A80T-32-aliquot",
              "HCM-SANG-1335-C18-85A-01R-A80W-32-aliquot",
              "HCM-SANG-1335-C18-10A-01D-A80T-32-aliquot"
            ],
            "created_datetime": "2022-11-07T08:35:49.832199-06:00",
            "diagnosis_ids": [
              "05c6c397-b06f-4fc3-806d-584e01bb976b"
            ],
            "sample_ids": [
              "92af994f-026c-4928-a1ab-95bb6a97b6f7",
              "a0c5e36b-b7ff-4cdb-89f6-1e84e752dc1a",
              "bfc0d425-2539-4fb1-970a-d5b7ca32aafb"
            ],
            "submitter_sample_ids": [
              "HCM-SANG-1335-C18-85A-01R-A80W-32",
              "HCM-SANG-1335-C18-85A-01D-A80T-32",
              "HCM-SANG-1335-C18-10A-01D-A80T-32"
            ],
            "primary_site": "Colon",
            "submitter_diagnosis_ids": [
              "HCM-SANG-1335-C18_diagnosis"
            ],
            "updated_datetime": "2025-04-07T17:13:09.336794-05:00",
            "case_id": "12603388-e5f0-45ac-ab7b-3a839de7a9e8",
            "index_date": "Sample Procurement",
            "state": "released"
          },
          {
            "id": "1312ba61-8d33-4b61-89cf-0c6e5c6d06a7",
            "disease_type": "Adenomas and Adenocarcinomas",
            "submitter_id": "HCM-SANG-1322-C15",
            "aliquot_ids": [
              "5bf6d83f-c7dd-4d4a-8e08-dd002340b190",
              "653d94fa-07fa-4a1b-9a00-6469027dec40",
              "9582d868-547e-4c37-863d-21cff2ad3afa",
              "077cd901-d3d4-4c35-8b28-a4c0c425500b"
            ],
            "submitter_aliquot_ids": [
              "HCM-SANG-1322-C15-85A-01D-A80T-32-aliquot",
              "HCM-SANG-1322-C15-85A-01R-A80W-32-aliquot",
              "HCM-SANG-1322-C15-10A-01D-A80T-32-aliquot",
              "HCM-SANG-1322-C15-08A-01D-A80T-32-aliquot"
            ],
            "created_datetime": "2022-11-07T08:35:49.832199-06:00",
            "diagnosis_ids": [
              "1b0aa02b-19d9-4e30-98d7-86d641180351"
            ],
            "sample_ids": [
              "85e0c249-4bd8-4311-8695-86fc7ab45212",
              "475c6bcb-3287-44ff-9128-7deda95902de",
              "b282ba81-bcfb-4e48-b847-1c5c897f3b39",
              "3d4e94f9-211e-4f08-9ba9-b0330ca96b7b"
            ],
            "submitter_sample_ids": [
              "HCM-SANG-1322-C15-08A-01D-A80T-32",
              "HCM-SANG-1322-C15-85A-01D-A80T-32",
              "HCM-SANG-1322-C15-85A-01R-A80W-32",
              "HCM-SANG-1322-C15-10A-01D-A80T-32"
            ],
            "primary_site": "Esophagus",
            "submitter_diagnosis_ids": [
              "HCM-SANG-1322-C15_diagnosis"
            ],
            "updated_datetime": "2025-04-07T17:13:09.336794-05:00",
            "case_id": "1312ba61-8d33-4b61-89cf-0c6e5c6d06a7",
            "index_date": "Sample Procurement",
            "state": "released"
          },
          {
            "id": "35a8dada-9a8b-4dd7-8565-e5e2c1e69958",
            "lost_to_followup": null,
            "slide_ids": [
              "8bbddca9-7db2-4514-94c8-40c98aecff02",
              "b4825a86-47ae-443e-9293-7892c982f4d1"
            ],
            "submitter_slide_ids": [
              "HCM-BROD-0046-C71-02A-01-S2-HE",
              "HCM-BROD-0046-C71-02A-01-S1-HE"
            ],
            "days_to_lost_to_followup": null,
            "disease_type": "Gliomas",
            "analyte_ids": [
              "5e713f2d-b3b6-40f9-80f1-fb53f87f7421",
              "7625dfef-8aba-453b-b830-237e6110ec26",
              "b964ed74-a27c-433c-9ad8-252ddee02b1b",
              "0ec84cd1-1957-45b1-a4ae-fdae0615de26"
            ],
            "submitter_id": "HCM-BROD-0046-C71",
            "submitter_analyte_ids": [
              "HCM-BROD-0046-C71-02A-11R",
              "HCM-BROD-0046-C71-85A-01R",
              "HCM-BROD-0046-C71-02A-11D",
              "HCM-BROD-0046-C71-85B-01D"
            ],
            "days_to_consent": null,
            "aliquot_ids": [
              "bf831b26-3782-410a-8634-f8401279d122",
              "6349f05b-3c22-4246-9cac-91e7ce73e203",
              "f95a9257-51af-4a94-a892-6a29a91d8528",
              "5c0f2139-0285-4ee4-9fb6-f7d696e9c372"
            ],
            "submitter_aliquot_ids": [
              "HCM-BROD-0046-C71-02A-11D-A78T-36",
              "HCM-BROD-0046-C71-02A-11R-A78V-41",
              "HCM-BROD-0046-C71-85A-01R-A78V-41",
              "HCM-BROD-0046-C71-85B-01D-A78T-36"
            ],
            "created_datetime": "2019-04-04T14:43:01.620260-05:00",
            "diagnosis_ids": [
              "bbe94910-fbe5-49bf-b984-2d2021bc1b30"
            ],
            "sample_ids": [
              "7cd03d20-6709-4f9e-9e02-1b68b8259050",
              "05f85ea2-6b34-40e6-8350-b344c4ac1369",
              "02312ea6-0a03-429b-a4bb-29c636e89397"
            ],
            "consent_type": null,
            "submitter_sample_ids": [
              "HCM-BROD-0046-C71-85A",
              "HCM-BROD-0046-C71-02A",
              "HCM-BROD-0046-C71-85B"
            ],
            "primary_site": "Brain",
            "submitter_diagnosis_ids": [
              "HCM-BROD-0046-C71_diagnosis"
            ],
            "updated_datetime": "2025-03-31T16:34:59.906497-05:00",
            "case_id": "35a8dada-9a8b-4dd7-8565-e5e2c1e69958",
            "index_date": "Diagnosis",
            "state": "released",
            "portion_ids": [
              "a9c741e5-6b80-4cdc-8da4-5d7f89aec18a"
            ],
            "submitter_portion_ids": [
              "HCM-BROD-0046-C71-02A-11"
            ]
          },
          {
            "id": "14d24b35-141b-4743-b4a9-321e2393bb54",
            "slide_ids": [
              "df3d99b3-27fd-44b3-a4dd-c0a7a3125b66",
              "7a634416-9cf2-45a9-8910-6a3db8d0a11a"
            ],
            "submitter_slide_ids": [
              "HCM-CSHL-0625-C18-06A-01-S2-HE",
              "HCM-CSHL-0625-C18-06A-01-S1-HE"
            ],
            "disease_type": "Adenomas and Adenocarcinomas",
            "analyte_ids": [
              "36d574ad-b36c-4691-8550-ffd332daf891",
              "5634ce55-23f6-4956-afe4-ff8b7def6f29",
              "dbc44a11-fcb6-481b-8d37-7c514c9738b8",
              "9467a0c7-96f7-4ebd-a617-94d06502d0a5",
              "7b547352-46f8-45b9-833c-df7144a99dbc"
            ],
            "submitter_id": "HCM-CSHL-0625-C18",
            "submitter_analyte_ids": [
              "HCM-CSHL-0625-C18-06A-11R",
              "HCM-CSHL-0625-C18-10A-01D",
              "HCM-CSHL-0625-C18-85M-01R",
              "HCM-CSHL-0625-C18-85M-01D",
              "HCM-CSHL-0625-C18-06A-11D"
            ],
            "aliquot_ids": [
              "7600c987-7246-4777-bc96-4eeb308b0c60",
              "325f86ae-e738-4d93-92b1-53b00a428e5b",
              "ee960355-6dc7-4807-aabd-7a130000262d",
              "c26e7443-e3a5-4018-863a-994e4ebeb132",
              "8a222af8-5fc8-4eb8-a50e-dd3aa42b9b40"
            ],
            "submitter_aliquot_ids": [
              "HCM-CSHL-0625-C18-85M-01D-A85K-36",
              "HCM-CSHL-0625-C18-10A-01D-A85K-36",
              "HCM-CSHL-0625-C18-06A-11R-A85M-41",
              "HCM-CSHL-0625-C18-06A-11D-A85K-36",
              "HCM-CSHL-0625-C18-85M-01R-A85M-41"
            ],
            "created_datetime": "2021-07-06T12:40:12.582209-05:00",
            "diagnosis_ids": [
              "dceb4b03-d5e1-4097-a188-e1fa3199f844",
              "ddc2cd05-9276-45b7-a629-974b02869272"
            ],
            "sample_ids": [
              "e2e51bb0-fa8c-4c50-b235-318c528436a3",
              "66df9961-1439-4f0c-ba39-b49652dcdeae",
              "1cc415d6-250f-4edf-952c-aa9fbfa5a52f"
            ],
            "submitter_sample_ids": [
              "HCM-CSHL-0625-C18-10A",
              "HCM-CSHL-0625-C18-06A",
              "HCM-CSHL-0625-C18-85M"
            ],
            "primary_site": "Colon",
            "submitter_diagnosis_ids": [
              "HCM-CSHL-0625-C18_diagnosis",
              "HCM-CSHL-0625-C18_diagnosis2"
            ],
            "updated_datetime": "2025-04-07T17:13:09.336794-05:00",
            "case_id": "14d24b35-141b-4743-b4a9-321e2393bb54",
            "index_date": "Diagnosis",
            "state": "released",
            "portion_ids": [
              "9b52cede-7745-4540-ab0e-998c1d2976fd",
              "34fe335e-0038-4141-962f-4eb00e69f875"
            ],
            "submitter_portion_ids": [
              "HCM-CSHL-0625-C18-10A-01",
              "HCM-CSHL-0625-C18-06A-11"
            ]
          },
          {
            "id": "157d4813-8461-4593-8232-863c7ace4b35",
            "slide_ids": [
              "76e8580c-a519-426c-a98d-84639a318cc9",
              "a68a9731-d510-4b80-963c-0e53a287218e",
              "d505d837-1e05-4b2e-85f1-232047339330",
              "262bd51e-189b-4250-9838-de97f462f762"
            ],
            "submitter_slide_ids": [
              "HCM-CSHL-0860-C18-06A-01-S2-HE",
              "HCM-CSHL-0860-C18-11A-S1-HE",
              "HCM-CSHL-0860-C18-11A-S2-HE",
              "HCM-CSHL-0860-C18-06A-01-S1-HE"
            ],
            "disease_type": "Adenomas and Adenocarcinomas",
            "analyte_ids": [
              "486be753-7b36-4b12-ab09-859006718773",
              "4f49f92d-aaaa-48d0-b609-83c9a294784d",
              "84834dac-99ee-4803-8709-fcfb79a49be4",
              "e08c3742-2db9-4278-bbbf-ef78cb979047",
              "70d65521-7e4e-4001-8806-a7f0415a656e"
            ],
            "submitter_id": "HCM-CSHL-0860-C18",
            "submitter_analyte_ids": [
              "HCM-CSHL-0860-C18-85M-01D",
              "HCM-CSHL-0860-C18-06A-11R",
              "HCM-CSHL-0860-C18-06A-11D",
              "HCM-CSHL-0860-C18-85M-01R",
              "HCM-CSHL-0860-C18-11A-01D"
            ],
            "aliquot_ids": [
              "eb040f56-135a-430b-bdbd-c45372f4d7fc",
              "20ac2691-66ad-4978-b56a-0ffb7484a4de",
              "f2ec41de-5a0b-40d0-a5e9-ce624102976b",
              "2027c333-dc19-4a63-bf94-8e7b9078ec60",
              "c4a846b4-2194-4280-9364-d33e8c6fb38e"
            ],
            "submitter_aliquot_ids": [
              "HCM-CSHL-0860-C18-06A-11R-A885-41",
              "HCM-CSHL-0860-C18-85M-01R-A885-41",
              "HCM-CSHL-0860-C18-11A-01D-A881-36",
              "HCM-CSHL-0860-C18-06A-11D-A881-36",
              "HCM-CSHL-0860-C18-85M-01D-A881-36"
            ],
            "created_datetime": "2021-11-10T16:06:04.214162-06:00",
            "diagnosis_ids": [
              "43034527-795f-4d93-8ce0-777dd3947050",
              "de5b2d58-3b24-4e7a-bc6f-8558fc7ec068"
            ],
            "sample_ids": [
              "205029b5-42b6-4adf-b4f0-857b7fc22fe0",
              "3f3e66ab-449d-445b-9b0d-3efe4407f2a7",
              "e68fe9c2-e70b-4c34-89e9-fd11bf608ac7"
            ],
            "submitter_sample_ids": [
              "HCM-CSHL-0860-C18-06A",
              "HCM-CSHL-0860-C18-11A",
              "HCM-CSHL-0860-C18-85M"
            ],
            "primary_site": "Colon",
            "submitter_diagnosis_ids": [
              "HCM-CSHL-0860-C18_diagnosis2",
              "HCM-CSHL-0860-C18_diagnosis"
            ],
            "updated_datetime": "2025-04-07T17:13:09.336794-05:00",
            "case_id": "157d4813-8461-4593-8232-863c7ace4b35",
            "index_date": "Diagnosis",
            "state": "released",
            "portion_ids": [
              "a2757cc9-313c-43b4-94b6-2d5b1f253672",
              "91801435-f54e-43d5-8c35-1f0560bfc49e"
            ],
            "submitter_portion_ids": [
              "HCM-CSHL-0860-C18-11A-01",
              "HCM-CSHL-0860-C18-06A-11"
            ]
          },
          {
            "id": "745d26a1-1072-4ac0-b168-03f4287c32b0",
            "slide_ids": [
              "939f1bd4-fc3b-4963-9306-5edbe819c855",
              "f435a03e-0517-40df-b540-17af191a7622"
            ],
            "submitter_slide_ids": [
              "HCM-BROD-0644-C25-06B-01-S1-HE",
              "HCM-BROD-0644-C25-06B-01-S2-HE"
            ],
            "disease_type": "Epithelial Neoplasms, NOS",
            "analyte_ids": [
              "5f8f96b4-14bf-44b1-8899-6623bdef942f",
              "3cf3e3eb-0577-46aa-a82d-890522427021",
              "469bec52-e0c3-4333-948f-9a2d6a784599",
              "ec0d7a33-f4ef-4b31-8485-264bcf4677af",
              "007e99c8-d719-4b0e-84e0-2d8482a31159"
            ],
            "submitter_id": "HCM-BROD-0644-C25",
            "submitter_analyte_ids": [
              "HCM-BROD-0644-C25-85M-01R",
              "HCM-BROD-0644-C25-06B-01D",
              "HCM-BROD-0644-C25-06B-01R",
              "HCM-BROD-0644-C25-10A-01D",
              "HCM-BROD-0644-C25-85M-01D"
            ],
            "aliquot_ids": [
              "e232849a-087b-4e98-883b-c8346527d430",
              "b922a932-004e-4c64-bb45-17c5a4f82784",
              "f356b6fb-9902-4420-a64a-1aa6a8eda705",
              "1ec464bb-2abd-4931-92fb-e33716f9d128",
              "3f0a435a-7475-44e5-8779-6777a44eae60"
            ],
            "submitter_aliquot_ids": [
              "HCM-BROD-0644-C25-85M-01D-A83U-36",
              "HCM-BROD-0644-C25-06B-01R-A83V-41",
              "HCM-BROD-0644-C25-85M-01R-A83V-41",
              "HCM-BROD-0644-C25-10A-01D-A83U-36",
              "HCM-BROD-0644-C25-06B-01D-A83U-36"
            ],
            "created_datetime": "2020-10-08T14:49:03.103433-05:00",
            "diagnosis_ids": [
              "87084624-73bb-47ea-9ae3-380158c67ac0",
              "c4383dd9-ba1a-480d-8fd1-6cda99501eb1"
            ],
            "sample_ids": [
              "23af2ae5-4b4c-4641-a182-e737c7a70ab9",
              "db97e33c-8e35-4905-ab27-b209383fbfc2",
              "517daee3-c532-483a-bd8a-e6c1099fe71c"
            ],
            "submitter_sample_ids": [
              "HCM-BROD-0644-C25-10A",
              "HCM-BROD-0644-C25-06B",
              "HCM-BROD-0644-C25-85M"
            ],
            "primary_site": "Pancreas",
            "submitter_diagnosis_ids": [
              "HCM-BROD-0644-C25_diagnosis",
              "HCM-BROD-0644-C25_diagnosis2"
            ],
            "updated_datetime": "2025-03-31T16:34:59.906497-05:00",
            "case_id": "745d26a1-1072-4ac0-b168-03f4287c32b0",
            "index_date": "Diagnosis",
            "state": "released",
            "portion_ids": [
              "0d8a3c18-7455-4a89-94af-f136c4b0e581",
              "2863a9c7-06fc-4333-9057-c559bc6ceb8b"
            ],
            "submitter_portion_ids": [
              "HCM-BROD-0644-C25-06B-01",
              "HCM-BROD-0644-C25-10A-01"
            ]
          }
        ],
        "pagination": {
          "count": 10,
          "total": 23793,
          "size": 10,
          "from": 0,
          "sort": "None",
          "page": 1,
          "pages": 2380
        }
      },
      "warnings": {}
    }
    ```



#### Example: HTTP POST Request

This example demonstrates how to obtain metadata in TSV format for a set of files using their UUIDs (e.g. UUIDs obtained from a [download manifest file generated by the GDC Data Portal](/Data_Portal/Users_Guide/Repository/#generating-a-manifest-file-for-the-data-transfer-tool)).

The first step is to construct a JSON query object, including `filters`, `fields`, `format`, and `size` parameters. The object is then submitted as HTTP POST payload to the GDC API using curl, in order to retrieve a TSV file with the requested metadata.

=== "Payload_txt"

    ```json
    {
        "filters":{
            "op":"in",
            "content":{
                "field":"files.file_id",
                "value":[
                    "b4bce3ff-7fdc-4849-880b-56f2b348ceac",
                    "5ca9fa79-53bc-4e91-82cd-5715038ee23e",
                    "b7c3e5ad-4ffc-4fc4-acbf-1dfcbd2e5382",
                    "1bc55036-2c7a-4333-87fd-b336056a8f06",
                    "6c00906c-7eb7-484d-93a9-d5b2075e7b50",
                    "b246b593-e1d1-4711-8aa5-d8e9eead9e2b",
                    "a2ee0837-f3a9-42be-ba26-f2bb1d6a50c0",
                    "72c1126d-384f-4f8b-becf-92e0779525b7",
                    "53984642-d431-490e-8a8a-c59b872ace66",
                    "b2c51c85-f0fd-4125-a4c4-2f8c3af35bf0",
                    "eeab9e3f-6158-4335-a5a0-7b4b79ac6056",
                    "f045d839-19e8-4d72-9ff3-809f00487934",
                    "e6ddb2b6-f137-488b-b171-c1748c089e15",
                    "8f69d768-5170-4245-b3dd-d57a14ef8bf9",
                    "4e1669f6-9f9c-4d0e-afd8-f4e29a7602af"
                ]
            }
        },
        "format":"TSV",
        "fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id",
        "size":"100"
    }
    ```

=== "Shell"

    ```shell
    curl --request POST --header "Content-Type: application/json" --data @Payload.txt 'https://api.gdc.cancer.gov/files' > File_metadata.txt
    ```

=== "File_metadata_txt"

    ```
    cases.0.case_id	cases.0.samples.0.portions.0.analytes.0.aliquots.0.aliquot_id	cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id	cases.0.samples.0.sample_id	cases.0.samples.0.sample_type	cases.0.samples.0.submitter_id	cases.0.samples.0.tissue_type	cases.0.samples.0.tumor_descriptor	cases.0.submitter_id	data_category	data_type	file_id	file_name	id
    f0d8a1fe-e313-44f1-99cc-b965cbeeff0e	922226ba-6244-4953-ad42-f4daa474c288	TCGA-C8-A27B-10A-01D-A167-09	31139082-7978-45aa-9d8f-ac4789ac5cec	Blood Derived Normal	TCGA-C8-A27B-10A	Normal	Not Applicable	TCGA-C8-A27B	Sequencing Reads	Aligned Reads	a2ee0837-f3a9-42be-ba26-f2bb1d6a50c0	922226ba-6244-4953-ad42-f4daa474c288_wxs_gdc_realn.bam	a2ee0837-f3a9-42be-ba26-f2bb1d6a50c0
    e62a728d-390f-428a-bea1-fc8c9814fb11	641a0220-6eec-434a-b606-e256113b65da	TCGA-CE-A484-10A-01D-A23U-08	27a8008e-044a-4966-b518-cc6905e292ca	Blood Derived Normal	TCGA-CE-A484-10A	Normal	Not Applicable	TCGA-CE-A484	Sequencing Reads	Aligned Reads	4e1669f6-9f9c-4d0e-afd8-f4e29a7602af	641a0220-6eec-434a-b606-e256113b65da_wxs_gdc_realn.bam	4e1669f6-9f9c-4d0e-afd8-f4e29a7602af
    29c8f468-5ac1-4d6c-8376-e36e6d246926	02e65074-ffda-4795-b8f5-1bfd20bd1019	TCGA-B5-A0K0-10A-01W-A062-09	1df69e2e-f392-465f-8e61-4671ba2fcd35	Blood Derived Normal	TCGA-B5-A0K0-10A	Normal	Not Applicable	TCGA-B5-A0K0	Sequencing Reads	Aligned Reads	53984642-d431-490e-8a8a-c59b872ace66	02e65074-ffda-4795-b8f5-1bfd20bd1019_wxs_gdc_realn.bam	53984642-d431-490e-8a8a-c59b872ace66
    b5c1e511-baf2-45b3-9919-110e8941e3c2	2a8cb8fe-b64f-453e-8139-7ede12f3fc51	TCGA-EC-A24G-10A-01D-A16D-09	61cf2e54-1b8d-40a0-9c73-a7449cbd570a	Blood Derived Normal	TCGA-EC-A24G-10A	Normal	Not Applicable	TCGA-EC-A24G	Sequencing Reads	Aligned Reads	eeab9e3f-6158-4335-a5a0-7b4b79ac6056	2a8cb8fe-b64f-453e-8139-7ede12f3fc51_wxs_gdc_realn.bam	eeab9e3f-6158-4335-a5a0-7b4b79ac6056
    24faa36a-268d-4a13-b3ae-eacd431a2bcc	c3feacc2-5a26-4bb2-a312-8b2ee53ccad1	TCGA-EE-A2GU-10A-01D-A198-08	cc4a5ed8-376a-4842-a25d-ffb07d8e1ca0	Blood Derived Normal	TCGA-EE-A2GU-10A	Normal	Not Applicable	TCGA-EE-A2GU	Sequencing Reads	Aligned Reads	8f69d768-5170-4245-b3dd-d57a14ef8bf9	c3feacc2-5a26-4bb2-a312-8b2ee53ccad1_wxs_gdc_realn.bam	8f69d768-5170-4245-b3dd-d57a14ef8bf9
    8fc9cc74-f388-49f0-b957-debb62638634	30919a1a-df9f-4604-835e-f66ac7bcacdf	TCGA-DA-A1IB-10A-01D-A198-08	432952c5-6505-4220-a581-f65270a45281	Blood Derived Normal	TCGA-DA-A1IB-10A	Normal	Not Applicable	TCGA-DA-A1IB	Sequencing Reads	Aligned Reads	f045d839-19e8-4d72-9ff3-809f00487934	30919a1a-df9f-4604-835e-f66ac7bcacdf_wxs_gdc_realn.bam	f045d839-19e8-4d72-9ff3-809f00487934
    a9255dcb-b236-4777-ac43-555e3a5386c3	9f4ffc2f-d006-4d86-b3b1-b25020481893	TCGA-EB-A3XB-10B-01D-A23B-08	0e1d4c7c-204d-4765-b090-68ed4cd83835	Blood Derived Normal	TCGA-EB-A3XB-10B	Normal	Not Applicable	TCGA-EB-A3XB	Sequencing Reads	Aligned Reads	b246b593-e1d1-4711-8aa5-d8e9eead9e2b	9f4ffc2f-d006-4d86-b3b1-b25020481893_wxs_gdc_realn.bam	b246b593-e1d1-4711-8aa5-d8e9eead9e2b
    a4225cb2-7b4b-4122-b6b9-629c26e3ea56	f4799bdc-b207-4053-9a4b-5a26ebf8ab91	TCGA-B9-A69E-10A-01D-A31X-10	5d6d6cd4-6a7b-499d-936a-1be9bf74b07f	Blood Derived Normal	TCGA-B9-A69E-10A	Normal	Not Applicable	TCGA-B9-A69E	Sequencing Reads	Aligned Reads	b2c51c85-f0fd-4125-a4c4-2f8c3af35bf0	f4799bdc-b207-4053-9a4b-5a26ebf8ab91_wxs_gdc_realn.bam	b2c51c85-f0fd-4125-a4c4-2f8c3af35bf0
    7a2cf5ce-8317-4fff-946e-b9937afab815	8c34ffe2-9012-4b4a-b610-a42a9c6a9780	TCGA-AX-A2HG-10A-01D-A17D-09	ef4b80ec-b453-48ec-8ad8-ccac83e1e4db	Blood Derived Normal	TCGA-AX-A2HG-10A	Normal	Not Applicable	TCGA-AX-A2HG	Sequencing Reads	Aligned Reads	e6ddb2b6-f137-488b-b171-c1748c089e15	8c34ffe2-9012-4b4a-b610-a42a9c6a9780_wxs_gdc_realn.bam	e6ddb2b6-f137-488b-b171-c1748c089e15
    e7a1cbe2-793c-4747-8412-8be794f2382b	66cbb40f-14b3-40c0-a332-e8a8e21bca11	TCGA-G7-6790-10A-01D-1962-08	4be83d0f-8b09-4e9e-8318-358371d34332	Blood Derived Normal	TCGA-G7-6790-10A	Normal	Not Applicable	TCGA-G7-6790	Sequencing Reads	Aligned Reads	6c00906c-7eb7-484d-93a9-d5b2075e7b50	66cbb40f-14b3-40c0-a332-e8a8e21bca11_wxs_gdc_realn.bam	6c00906c-7eb7-484d-93a9-d5b2075e7b50
    c787c4da-c564-44f1-89eb-dd9da107acb1	c723584a-c404-4c88-bfea-e40f5dbba542	TCGA-EB-A44O-10A-01D-A25O-08	5b738547-1825-4684-81bd-864bf2eb43ef	Blood Derived Normal	TCGA-EB-A44O-10A	Normal	Not Applicable	TCGA-EB-A44O	Sequencing Reads	Aligned Reads	72c1126d-384f-4f8b-becf-92e0779525b7	c723584a-c404-4c88-bfea-e40f5dbba542_wxs_gdc_realn.bam	72c1126d-384f-4f8b-becf-92e0779525b7
    fec0da58-1047-44d2-b6d1-c18cceed43dc	cd761feb-9a20-4495-8943-c6243532a5cf	TCGA-E9-A295-10A-01D-A16D-09	e74183e1-f0b4-412a-8dac-a62d404add78	Blood Derived Normal	TCGA-E9-A295-10A	Normal	Not Applicable	TCGA-E9-A295	Sequencing Reads	Aligned Reads	5ca9fa79-53bc-4e91-82cd-5715038ee23e	cd761feb-9a20-4495-8943-c6243532a5cf_wxs_gdc_realn.bam	5ca9fa79-53bc-4e91-82cd-5715038ee23e
    53886143-c1c6-40e9-88e6-e4e5e0271fc8	e96d5811-4736-40dd-966d-e0e172aeb0af	TCGA-A2-A3XX-10A-01D-A23C-09	c6eb6218-ad71-40a6-88b7-a4f1a015b816	Blood Derived Normal	TCGA-A2-A3XX-10A	Normal	Not Applicable	TCGA-A2-A3XX	Sequencing Reads	Aligned Reads	1bc55036-2c7a-4333-87fd-b336056a8f06	e96d5811-4736-40dd-966d-e0e172aeb0af_wxs_gdc_realn.bam	1bc55036-2c7a-4333-87fd-b336056a8f06
    8aaa4e25-5c12-4ace-96dc-91aaa0c4457c	b4e4630a-b38c-4b62-b0e8-d73f0e3b4e47	TCGA-B0-5094-11A-01D-1421-08	7519d7a8-c3ee-417b-9cfc-111bc5ad0637	Solid Tissue Normal	TCGA-B0-5094-11A	Normal	Not Applicable	TCGA-B0-5094	Sequencing Reads	Aligned Reads	b4bce3ff-7fdc-4849-880b-56f2b348ceac	b4e4630a-b38c-4b62-b0e8-d73f0e3b4e47_wxs_gdc_realn.bam	b4bce3ff-7fdc-4849-880b-56f2b348ceac
    ae55b2d3-62a1-419e-9f9a-5ddfac356db4	45c68b6b-0bed-424d-9a77-4f87bbaa3649	TCGA-B0-5117-11A-01D-1421-08	b1116541-bece-4df3-b3dd-cec50aeb277b	Solid Tissue Normal	TCGA-B0-5117-11A	Normal	Not Applicable	TCGA-B0-5117	Sequencing Reads	Aligned Reads	b7c3e5ad-4ffc-4fc4-acbf-1dfcbd2e5382	45c68b6b-0bed-424d-9a77-4f87bbaa3649_wxs_gdc_realn.bam	b7c3e5ad-4ffc-4fc4-acbf-1dfcbd2e5382

    ```


### Format

Specifies the format of the API response: JSON (default), `TSV` or `XML`.

#### Examples

=== "Shell1"

    ```shell
    curl  'https://api.gdc.cancer.gov/cases?fields=submitter_id&size=5&format=TSV'
    ```

=== "Python1"

    ```python
    import requests

    cases_endpt = 'https://api.gdc.cancer.gov/cases'
    params = {'fields':'submitter_id',
              'format':'TSV'}
    response = requests.get(cases_endpt, params = params)
    print(response.content)
    ```

=== "Response1"

    ```
    id	submitter_id
    0286c31b-a704-4d7d-99e3-0bc4e8975b8b	HCM-CSHL-0084-C25
    02f6d684-b6b5-419a-b0e1-b74d0a384a30	HCM-BROD-0408-C71
    03974dc9-0162-4de8-9897-09f88693681a	HCM-BROD-0334-C43
    03bfeb7c-cecf-4691-8263-33cdfe391ea9	HCM-BROD-0124-C25
    04cbceab-f945-482b-956b-840756a17a4a	HCM-BROD-0421-C71

    ```

=== "Shell 2"

    ```shell
    curl  'https://api.gdc.cancer.gov/cases?fields=submitter_id&size=5&format=XML&pretty=true'
    ```

=== "Python2"

    ```python
    import requests

    cases_endpt = 'https://api.gdc.cancer.gov/cases'
    params = {'fields':'submitter_id',
              'format':'XML',
              'pretty':'true'}
    response = requests.get(cases_endpt, params = params)
    print(response.content)
    ```

=== "Output2"

    ```xml
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

=== "Request1"

    ```shell
    curl  'https://api.gdc.cancer.gov/cases?fields=submitter_id&sort=submitter_id:asc&size=5'
    ```

=== "Response1"

    ```json
    {"data": {"hits": [{"id": "be37f1f7-2f98-4f74-bc04-6dd2ae2afcad", "submitter_id": "01BR001"}, {"id": "e6915db0-7c89-484d-8f9f-15cca68b82fc", "submitter_id": "01BR008"}, {"id": "16614d46-172b-479c-992b-e80a8e9a2c59", "submitter_id": "01BR009"}, {"id": "567fc9e3-17a6-42b1-a896-5e9a9507d1d8", "submitter_id": "01BR010"}, {"id": "54e89878-a1bc-4f5a-9d68-4842a469586e", "submitter_id": "01BR015"}], "pagination": {"count": 5, "total": 86962, "size": 5, "from": 0, "sort": "None", "page": 1, "pages": 17393}}, "warnings": {}}
    ```

=== "Request2"

    ```shell
    curl  'https://api.gdc.cancer.gov/cases?fields=submitter_id&sort=submitter_id:asc&size=5&pretty=true'
    ```

=== "Response2"

    ```json
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

=== "Shell"

    ```shell
    curl 'https://api.gdc.cancer.gov/files?fields=cases.submitter_id,file_id,file_name,file_size&pretty=true'
    ```

=== "Python"    

    ```python
    import requests
    import json
    
    files_endpt = 'https://api.gdc.cancer.gov/files'
    params = {'fields':'cases.submitter_id,file_id,file_name,file_size'}
    response = requests.get(files_endpt, params = params)
    print(json.dumps(response.json(), indent=2))
    ```

=== "Response"

    ```json
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

The `expand` parameter provides a shortcut to request multiple related fields (field groups) in the response. Instead of specifying each field using the `fields` parameter, users can specify a field group name using the `expand` parameter to request all fields in the group. Available field groups are listed in [Appendix A](Appendix_A_Available_Fields.md#field-group-listing-by-endpoint); the list can also be accessed programmatically at the [_mapping endpoint](#_mapping-endpoint). The `fields` and `expand` parameters can be used together to request custom combinations of field groups and individual fields.

#### Example

=== "Shell"

    ```Shell
    curl 'https://api.gdc.cancer.gov/files/573ee7e9-b8bd-419e-808b-a027c4311731?expand=cases.samples&pretty=true'
    ```

=== "Response"

    ```json
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

=== "Shell1"

    ``` Shell
    curl 'https://api.gdc.cancer.gov/files?fields=file_name&from=0&size=2&pretty=true'
    ```

===  "Python1"

    ``` Python
    import requests
    import json

    files_endpt = 'https://api.gdc.cancer.gov/files'
    params = {'fields':'file_name',
              'from':0, 'size':2}
    response = requests.get(files_endpt, params = params)
    print(json.dumps(response.json(), indent=2))

    ```

=== "Response1"

    ```json
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

=== "Shell2"

    ``` Shell
    curl 'https://api.gdc.cancer.gov/files?fields=file_name&from=101&size=5&pretty=true'
    ```

=== "Python2"

    ``` Python
    import requests
    import json

    files_endpt = 'https://api.gdc.cancer.gov/files'
    params = {'fields':'file_name',
              'from':101, 'size':5}
    response = requests.get(files_endpt, params = params)
    print(json.dumps(response.json(), indent=2))
    ```

=== "Output2"

    ```json
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

The `sort` query parameter sorts the results by a specific field, and with the sort direction specified using the `:asc` (ascending) or `:desc` (descending) prefix, e.g. `sort=field:desc`. A list of all valid _field_ names is available in [Appendix A](Appendix_A_Available_Fields.md); the list can also be accessed programmatically at the [_mapping endpoint](#_mapping-endpoint).

#### Example

Sort cases by `submitter_id` in ascending order:

=== "Shell"

    ``` shell
    curl  'https://api.gdc.cancer.gov/cases?fields=submitter_id&sort=submitter_id:asc&pretty=true'
    ```

=== "Python"

    ``` python
    import requests
    import json

    cases_endpt = 'https://api.gdc.cancer.gov/cases'
    params = {'fields':'submitter_id',
              'sort':'submitter_id:asc'}
    response = requests.get(cases_endpt, params = params)
    print(json.dumps(response.json(), indent=2))

    ```

=== "Output"
    
    ``` json
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

=== "Shell"

    ```shell
    curl  'https://api.gdc.cancer.gov/projects?facets=program.name&from=0&size=0&sort=program.name:asc&pretty=true'
    ```

=== "Python"

    ```python
    import requests
    import json

    projects_endpt = 'https://api.gdc.cancer.gov/projects'
    params = {'facets':'program.name',
              'from':0, 'size':0,
              'sort':'program.name:asc'}
    response = requests.get(projects_endpt, params = params)
    print(json.dumps(response.json(), indent=2))
    ```

=== "Response"

    ```
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
                "doc_count": 10,
                "key": "MATCH"
              },
              {
                "doc_count": 9,
                "key": "TARGET"
              },
              {
                "doc_count": 4,
                "key": "CGCI"
              },
              {
                "doc_count": 3,
                "key": "CMI"
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
                "doc_count": 2,
                "key": "MP2PRT"
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
          "total": 79,
          "size": 0,
          "from": 0,
          "sort": "None",
          "page": 1,
          "pages": 79
        }
      },
      "warnings": {}
    }
    ```

#### Example

In this sample POST request, both `filters` and `facets` parameters are used. Note that `facets` ignores the `primary_site` filter.

=== "Payload"

    ```json
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

=== "Shell"

    ```Shell
    curl --request POST --header "Content-Type: application/json" --data @Payload 'https://api.gdc.cancer.gov/v0/cases'
    ```

=== "Response"

    ```json
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
	filters=%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.file_id%22%2C%22value%22%3A%5B%22b4bce3ff-7fdc-4849-880b-56f2b348ceac%22%2C%225ca9fa79-53bc-4e91-82cd-5715038ee23e%22%2C%22b7c3e5ad-4ffc-4fc4-acbf-1dfcbd2e5382%22%2C%221bc55036-2c7a-4333-87fd-b336056a8f06%22%2C%226c00906c-7eb7-484d-93a9-d5b2075e7b50%22%2C%22b246b593-e1d1-4711-8aa5-d8e9eead9e2b%22%2C%22a2ee0837-f3a9-42be-ba26-f2bb1d6a50c0%22%2C%2272c1126d-384f-4f8b-becf-92e0779525b7%22%2C%2253984642-d431-490e-8a8a-c59b872ace66%22%2C%22b2c51c85-f0fd-4125-a4c4-2f8c3af35bf0%22%2C%22eeab9e3f-6158-4335-a5a0-7b4b79ac6056%22%2C%22f045d839-19e8-4d72-9ff3-809f00487934%22%2C%22e6ddb2b6-f137-488b-b171-c1748c089e15%22%2C%228f69d768-5170-4245-b3dd-d57a14ef8bf9%22%2C%224e1669f6-9f9c-4d0e-afd8-f4e29a7602af%22%5D%7D%7D&fields=file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id&format=tsv&size=100
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

The GDC Portal has a quicksearch functionality that allows for a project, case, or file to be queried from a search box. This function calls the `/v0/all` endpoint, which retrieves the top cases, files, projects, genes, mutations, and annotations that match to the query. The quicksearch can also be used programmatically through the API.  For example, a search term of 'TCGA' would produce the following query:  

=== "Shell"

    ```Shell
    curl "https://api.gdc.cancer.gov/v0/all?query=TCGA&size=5"
    ```

=== "Response"   

    ```json
    {"data":{"query":{"hits":[{"disease_type":["Adenomas and Adenocarcinomas"],"id":"UHJvamVjdDpUQ0dBLUFDQw==","name":"Adrenocortical Carcinoma","primary_site":["Adrenal gland"],"project_id":"TCGA-ACC","project_quicksearch":"Adrenocortical Carcinoma"},{"disease_type":["Adenomas and Adenocarcinomas"],"id":"UHJvamVjdDpUQ0dBLUtJQ0g=","name":"Kidney Chromophobe","primary_site":["Kidney"],"project_id":"TCGA-KICH","project_quicksearch":"Kidney Chromophobe"},{"disease_type":["Adenomas and Adenocarcinomas"],"id":"UHJvamVjdDpUQ0dBLUxJSEM=","name":"Liver Hepatocellular Carcinoma","primary_site":["Liver and intrahepatic bile ducts"],"project_id":"TCGA-LIHC","project_quicksearch":"Liver Hepatocellular Carcinoma"},{"disease_type":["Myeloid Leukemias"],"id":"UHJvamVjdDpUQ0dBLUxBTUw=","name":"Acute Myeloid Leukemia","primary_site":["Hematopoietic and reticuloendothelial systems"],"project_id":"TCGA-LAML","project_quicksearch":"Acute Myeloid Leukemia"},{"disease_type":["Adenomas and Adenocarcinomas"],"id":"UHJvamVjdDpUQ0dBLUtJUlA=","name":"Kidney Renal Papillary Cell Carcinoma","primary_site":["Kidney"],"project_id":"TCGA-KIRP","project_quicksearch":"Kidney Renal Papillary Cell Carcinoma"}],"total":183550}}}
    ```

This endpoint can be used to quickly retrieve information about a file.  For example, if a user wanted to know the UUID for `nationwidechildrens.org_biospecimen.TCGA-EL-A4K1.xml`, the following query could be used to quickly retrieve it programmatically:

=== "Shell"

    ```Shell
    curl "https://api.gdc.cancer.gov/v0/all?query=nationwidechildrens.org_biospecimen.TCGA-EL-A4K1.xml&size=5"
    ```

=== "Response"

    ```json
    {"data":{"query":{"hits":[{"file_id":"a74abfec-db78-4ed4-9e4b-604b66e30e30","file_name":"nationwidechildrens.org_biospecimen.TCGA-EL-A4K1.xml","id":"RmlsZTphNzRhYmZlYy1kYjc4LTRlZDQtOWU0Yi02MDRiNjZlMzBlMzA=","submitter_id":"nationwidechildrens.org_biospecimen.TCGA-EL-A4K1.xml"}],"total":1}}}
    ```

## Additional Examples

More examples of API functionality described in this section are provided in [Additional Examples](Additional_Examples.md).
