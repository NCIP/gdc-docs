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

>**Note:** Requests for information stored in the GDC Legacy Archive must be directed to `legacy/` endpoints. See [Getting Started](Getting_Started.md#gdc-legacy-archive) for details.

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
curl 'https://api.gdc.cancer.gov/projects?from=0&size=2&sort=project.project_id:asc&pretty=true'
```
``` Output
{
  "data": {
    "hits": [
      {
        "dbgap_accession_number": null,
        "disease_type": [
          "Brain Lower Grade Glioma"
        ],
        "released": true,
        "state": "legacy",
        "primary_site": [
          "Brain"
        ],
        "project_id": "TCGA-LGG",
        "id": "TCGA-LGG",
        "name": "Brain Lower Grade Glioma"
      },
      {
        "dbgap_accession_number": null,
        "disease_type": [
          "Thyroid Carcinoma"
        ],
        "released": true,
        "state": "legacy",
        "primary_site": [
          "Thyroid"
        ],
        "project_id": "TCGA-THCA",
        "id": "TCGA-THCA",
        "name": "Thyroid Carcinoma"
      }
    ],
    "pagination": {
      "count": 2,
      "sort": "project.project_id:asc",
      "from": 0,
      "page": 1,
      "total": 39,
      "pages": 20,
      "size": 2
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
    "dbgap_accession_number": "phs000467",
    "disease_type": [
      "Neuroblastoma"
    ],
    "summary": {
      "data_categories": [
        {
          "case_count": 151,
          "file_count": 471,
          "data_category": "Transcriptome Profiling"
        },
        {
          "case_count": 1127,
          "file_count": 3,
          "data_category": "Biospecimen"
        },
        {
          "case_count": 216,
          "file_count": 1732,
          "data_category": "Simple Nucleotide Variation"
        },
        {
          "case_count": 7,
          "file_count": 1,
          "data_category": "Clinical"
        },
        {
          "case_count": 270,
          "file_count": 599,
          "data_category": "Raw Sequencing Data"
        }
      ],
      "case_count": 1127,
      "file_count": 2806,
      "experimental_strategies": [
        {
          "case_count": 221,
          "file_count": 2174,
          "experimental_strategy": "WXS"
        },
        {
          "case_count": 151,
          "file_count": 628,
          "experimental_strategy": "RNA-Seq"
        }
      ],
      "file_size": 8157614402888
    },
    "released": true,
    "state": "legacy",
    "primary_site": [
      "Nervous System"
    ],
    "project_id": "TARGET-NBL",
    "name": "Neuroblastoma"
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
        "data_release": "13.0",
        "data_type": "Raw Simple Somatic Mutation",
        "updated_datetime": "2018-07-20T22:27:55.342974+00:00",
        "file_name": "333193d5-ca9a-4262-81f5-e9f3b44358fe.vcf.gz",
        "submitter_id": "AD19_SimpleSomaticMutation",
        "file_id": "333193d5-ca9a-4262-81f5-e9f3b44358fe",
        "file_size": 866,
        "id": "333193d5-ca9a-4262-81f5-e9f3b44358fe",
        "created_datetime": "2017-09-10T19:16:02.549312-05:00",
        "md5sum": "e33e95edb778fe67643162ef0ae3297e",
        "data_format": "VCF",
        "acl": [
          "phs001179"
        ],
        "access": "controlled",
        "state": "released",
        "version": "1",
        "data_category": "Simple Nucleotide Variation",
        "type": "simple_somatic_mutation",
        "experimental_strategy": "Targeted Sequencing"
      },
      {
        "data_release": "13.0",
        "data_type": "Raw Simple Somatic Mutation",
        "updated_datetime": "2018-07-20T22:27:55.342974+00:00",
        "file_name": "d9114e23-0f62-4979-aefc-0dd4d5eb891b.vcf.gz",
        "submitter_id": "AD116_SimpleSomaticMutation",
        "file_id": "d9114e23-0f62-4979-aefc-0dd4d5eb891b",
        "file_size": 866,
        "id": "d9114e23-0f62-4979-aefc-0dd4d5eb891b",
        "created_datetime": "2017-09-10T21:53:02.376246-05:00",
        "md5sum": "95bbfd0586d3c284e9f88edf3bf26065",
        "data_format": "VCF",
        "acl": [
          "phs001179"
        ],
        "access": "controlled",
        "state": "released",
        "version": "1",
        "data_category": "Simple Nucleotide Variation",
        "type": "simple_somatic_mutation",
        "experimental_strategy": "Targeted Sequencing"
      }
    ],
    "pagination": {
      "count": 2,
      "sort": "file_size:asc",
      "from": 0,
      "page": 1,
      "total": 356381,
      "pages": 178191,
      "size": 2
    }
  },
  "warnings": {}
}
```

#### Retrieval of file metadata using individual UUIDs:

The `/files` endpoint supports a simple query format that retrieves the metadata of a single file using its UUID.  Note that the `/files` endpoint is inactive when querying for earlier file versions.  In that case, the `/history` or `/files/versions` endpoints should be used instead.

```Shell
curl 'https://api.gdc.cancer.gov/files/874e71e0-83dd-4d3e-8014-10141b49f12c?pretty=true'
```
``` Output
{
  "data": {
    "data_release": "13.0",
    "data_type": "Raw Simple Somatic Mutation",
    "updated_datetime": "2018-07-20T22:27:55.342974+00:00",
    "created_datetime": "2016-06-03T17:03:06.608739-05:00",
    "file_name": "874e71e0-83dd-4d3e-8014-10141b49f12c.vcf.gz",
    "md5sum": "acf2929b1b825bcd1377023e8b8767ec",
    "data_format": "VCF",
    "acl": [
      "phs000178"
    ],
    "access": "controlled",
    "state": "live",
    "version": "1",
    "file_id": "874e71e0-83dd-4d3e-8014-10141b49f12c",
    "data_category": "Simple Nucleotide Variation",
    "file_size": 122293,
    "submitter_id": "TCGA-V4-A9EZ-01A-11D-A39W-08_TCGA-V4-A9EZ-10A-01D-A39Z-08_mutect",
    "type": "simple_somatic_mutation",
    "experimental_strategy": "WXS"
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
    "latest_size": 332092,
    "latest_id": "1dd28069-5777-4ff9-bd2b-d1ba68e88b06",
    "latest_version": "1",
    "filename": "1dd28069-5777-4ff9-bd2b-d1ba68e88b06.vcf.gz",
    "state": "validated",
    "version": "1",
    "latest_filename": "1dd28069-5777-4ff9-bd2b-d1ba68e88b06.vcf.gz",
    "latest_release": [
      "13.0"
    ],
    "latest_state": "validated",
    "release": "13.0",
    "latest_md5": "c2f9b196e154906a70c7ec46492a859d",
    "size": 332092,
    "id": "1dd28069-5777-4ff9-bd2b-d1ba68e88b06",
    "md5": "c2f9b196e154906a70c7ec46492a859d"
  },
  {
    "latest_size": 6653119038,
    "latest_id": "2a03abac-f1a2-49a9-a57c-7543739dd862",
    "latest_version": "1",
    "filename": "a5d86cde-32ca-4ed6-b1a5-5a47575f2ac6_gdc_realn_rehead.bam",
    "state": "validated",
    "version": "1",
    "latest_filename": "a5d86cde-32ca-4ed6-b1a5-5a47575f2ac6_gdc_realn_rehead.bam",
    "latest_release": [
      "13.0"
    ],
    "latest_state": "validated",
    "release": "13.0",
    "latest_md5": "48686fcd84ac713d44261ca9e26b89fb",
    "size": 6653119038,
    "id": "2a03abac-f1a2-49a9-a57c-7543739dd862",
    "md5": "48686fcd84ac713d44261ca9e26b89fb"
  }
]
```
```Shell2
curl --request POST --header "Content-Type: text/tsv"  https://api.gdc.cancer.gov/files/versions/manifest?pretty=true --data-binary @gdc_manifest_20180809_154816.txt
```

``` Output2
[{
  "latest_size": 44857,
  "state": "validated",
  "latest_version": "1",
  "filename": "nationwidechildrens.org_clinical.TCGA-13-1500.xml",
  "latest_id": "0b20e27c-9a09-4f15-923f-d5b4f185dc22",
  "version": "1",
  "latest_filename": "nationwidechildrens.org_clinical.TCGA-13-1500.xml",
  "latest_release": [
    "12.0"
  ],
  "latest_state": "validated",
  "release": "12.0",
  "latest_md5": "597aa4df24c4d544b6c25cbd8b25a33e",
  "md5": "597aa4df24c4d544b6c25cbd8b25a33e",
  "id": "0b20e27c-9a09-4f15-923f-d5b4f185dc22",
  "size": 44857
},{
  "latest_size": 27620,
  "state": "validated",
  "latest_version": "1",
  "filename": "BUCKS_p_TCGA_272_273_N_GenomeWideSNP_6_G05_1320676.grch38.seg.v2.txt",
  "latest_id": "3edc7084-013c-4493-8507-c00b0e9962d8",
  "version": "1",
  "latest_filename": "BUCKS_p_TCGA_272_273_N_GenomeWideSNP_6_G05_1320676.grch38.seg.v2.txt",
  "latest_release": [
    "12.0"
  ],
  "latest_state": "validated",
  "release": "12.0",
  "latest_md5": "35a18d990a05eedfaf96e753bee0b96d",
  "md5": "35a18d990a05eedfaf96e753bee0b96d",
  "id": "3edc7084-013c-4493-8507-c00b0e9962d8",
  "size": 27620
},{
  "latest_size": 2346,
  "state": "validated",
  "latest_version": "1",
  "filename": "a22f5e32-b16e-458f-a412-7e438056ece6.vep.vcf.gz",
  "latest_id": "a22f5e32-b16e-458f-a412-7e438056ece6",
  "version": "1",
  "latest_filename": "a22f5e32-b16e-458f-a412-7e438056ece6.vep.vcf.gz",
  "latest_release": [
    "12.0"
  ],
  "latest_state": "validated",
  "release": "12.0",
  "latest_md5": "68b2433b31679bbbc6681919a1b81762",
  "md5": "68b2433b31679bbbc6681919a1b81762",
  "id": "a22f5e32-b16e-458f-a412-7e438056ece6",
  "size": 2346
},{
  "latest_size": 35411,
  "state": "validated",
  "latest_version": "1",
  "filename": "CYANS_p_TCGAb_422_423_424_NSP_GenomeWideSNP_6_G12_1513758.nocnv_grch38.seg.v2.txt",
  "latest_id": "ac7d2078-bd6b-446e-b30a-d889da5624b6",
  "version": "1",
  "latest_filename": "CYANS_p_TCGAb_422_423_424_NSP_GenomeWideSNP_6_G12_1513758.nocnv_grch38.seg.v2.txt",
  "latest_release": [
    "12.0"
  ],
  "latest_state": "validated",
  "release": "12.0",
  "latest_md5": "6338826b620773062232830fad51ae64",
  "md5": "6338826b620773062232830fad51ae64",
  "id": "ac7d2078-bd6b-446e-b30a-d889da5624b6",
  "size": 35411
}]
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
	{
	  "data": {
	    "hits": [
	      {
	        "updated_datetime": "2017-03-04T16:39:19.244769-06:00",
	        "submitter_analyte_ids": [
	          "TCGA-BH-A0EA-01A-11R",
	          "TCGA-BH-A0EA-10A-01W",
	          "TCGA-BH-A0EA-01A-11W",
	          "TCGA-BH-A0EA-01A-11D",
	          "TCGA-BH-A0EA-10A-01D"
	        ],
	        "analyte_ids": [
	          "fe678556-acf4-4bde-a95e-860bb0150a95",
	          "66ed0f86-5ca5-4dec-ba76-7ee4dcf31831",
	          "f19f408a-815f-43d9-8032-e9482b796371",
	          "69ddc092-88a0-4839-a2bb-9f1c9e760409",
	          "30cb470f-66d4-4085-8c30-83a42e8453d4"
	        ],
	        "submitter_id": "TCGA-BH-A0EA",
	        "case_id": "1f601832-eee3-48fb-acf5-80c4a454f26e",
	        "id": "1f601832-eee3-48fb-acf5-80c4a454f26e",
	        "disease_type": "Breast Invasive Carcinoma",
	        "sample_ids": [
	          "9a6c71a6-82cd-42b1-a93f-f569370848d6",
	          "7f791228-dd77-4ab0-8227-d784a4c7fea1"
	        ],
	        "portion_ids": [
	          "cb6086d1-3416-4310-b109-e8fa6e8b72d4",
	          "8629bf5a-cdaf-4f6a-90bb-27dd4a7565c5",
	          "ae4f5816-f97a-4605-9b05-9ab820467dee"
	        ],
	        "submitter_portion_ids": [
	          "TCGA-BH-A0EA-01A-21-A13C-20",
	          "TCGA-BH-A0EA-01A-11",
	          "TCGA-BH-A0EA-10A-01"
	        ],
	        "created_datetime": null,
	        "slide_ids": [
	          "90154ea1-6b76-4445-870e-d531d6fa1239",
	          "a0826f0d-986a-491b-8c6f-b34f8929f3ee"
	        ],
	        "state": "live",
	        "aliquot_ids": [
	          "eef9dce1-6ba6-432b-bbe2-53c7dbe64fe7",
	          "cde982b7-3b0a-49eb-8710-a599cb0e44c1",
	          "b1a3739d-d554-4202-b96f-f25a444e2042",
	          "97c64d6a-7dce-4d0f-9cb3-b3e4eb4719c5",
	          "561b8777-801a-49ed-a306-e7dafeb044b6",
	          "42d050e4-e8ee-4442-b9c0-0ee14706b138",
	          "ca71ca96-cbb7-4eab-9487-251dda34e107",
	          "cfbd5476-e83a-401d-9f9a-639c73a0e35b",
	          "edad5bd3-efe0-4c5f-b05c-2c0c2951c45a",
	          "262715e1-835c-4f16-8ee7-6900e26f7cf5",
	          "2beb34c4-d493-4a73-b21e-de77d43251ff",
	          "bcb7fc6d-60a0-48b7-aa81-14c0dda72d76"
	        ],
	        "primary_site": "Breast",
	        "submitter_aliquot_ids": [
	          "TCGA-BH-A0EA-10A-01D-A113-01",
	          "TCGA-BH-A0EA-01A-11R-A115-07",
	          "TCGA-BH-A0EA-01A-11D-A10Y-09",
	          "TCGA-BH-A0EA-01A-11D-A314-09",
	          "TCGA-BH-A0EA-01A-11R-A114-13",
	          "TCGA-BH-A0EA-01A-11D-A111-01",
	          "TCGA-BH-A0EA-01A-11D-A112-05",
	          "TCGA-BH-A0EA-01A-11D-A10X-02",
	          "TCGA-BH-A0EA-10A-01D-A110-09",
	          "TCGA-BH-A0EA-10A-01W-A12U-09",
	          "TCGA-BH-A0EA-10A-01D-A10Z-02",
	          "TCGA-BH-A0EA-01A-11W-A12T-09"
	        ],
	        "submitter_sample_ids": [
	          "TCGA-BH-A0EA-10A",
	          "TCGA-BH-A0EA-01A"
	        ],
	        "submitter_slide_ids": [
	          "TCGA-BH-A0EA-01A-01-MSA",
	          "TCGA-BH-A0EA-01A-01-TSA"
	        ]
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

#### Retrieval of case metadata using individual UUIDs:

The `cases` endpoint supports a simple query format that retrieves the metadata of a single case using its UUID:

```shell
curl 'https://api.gdc.cancer.gov/cases/1f601832-eee3-48fb-acf5-80c4a454f26e?pretty=true&expand=diagnoses'
```
```Response
{
  "data": {
    "diagnoses": [
      {
        "classification_of_tumor": "not reported",
        "last_known_disease_status": "not reported",
        "updated_datetime": "2016-05-16T10:59:16.740358-05:00",
        "primary_diagnosis": "c50.9",
        "submitter_id": "TCGA-BH-A0EA_diagnosis",
        "tumor_stage": "stage iia",
        "age_at_diagnosis": 26548.0,
        "vital_status": "dead",
        "morphology": "8500/3",
        "days_to_death": 991.0,
        "days_to_last_known_disease_status": null,
        "days_to_last_follow_up": null,
        "state": null,
        "days_to_recurrence": null,
        "diagnosis_id": "84654ad5-2a2c-5c3b-8340-ecac6a5550fe",
        "tumor_grade": "not reported",
        "tissue_or_organ_of_origin": "c50.9",
        "days_to_birth": -26548.0,
        "progression_or_recurrence": "not reported",
        "prior_malignancy": "not reported",
        "site_of_resection_or_biopsy": "c50.9",
        "created_datetime": null
      }
    ],
    "sample_ids": [
      "7f791228-dd77-4ab0-8227-d784a4c7fea1",
      "9a6c71a6-82cd-42b1-a93f-f569370848d6"
    ],
    "portion_ids": [
      "cb6086d1-3416-4310-b109-e8fa6e8b72d4",
      "8629bf5a-cdaf-4f6a-90bb-27dd4a7565c5",
      "ae4f5816-f97a-4605-9b05-9ab820467dee"
    ],
    "submitter_portion_ids": [
      "TCGA-BH-A0EA-01A-11",
      "TCGA-BH-A0EA-01A-21-A13C-20",
      "TCGA-BH-A0EA-10A-01"
    ],
    "created_datetime": null,
    "submitter_aliquot_ids": [
      "TCGA-BH-A0EA-01A-11R-A114-13",
      "TCGA-BH-A0EA-01A-11D-A111-01",
      "TCGA-BH-A0EA-01A-11W-A12T-09",
      "TCGA-BH-A0EA-01A-11R-A114-13",
      "TCGA-BH-A0EA-01A-11R-A115-07",
      "TCGA-BH-A0EA-01A-11D-A111-01",
      "TCGA-BH-A0EA-01A-11D-A314-09",
      "TCGA-BH-A0EA-01A-11D-A112-05",
      "TCGA-BH-A0EA-01A-11D-A10Y-09",
      "TCGA-BH-A0EA-01A-11D-A10X-02",
      "TCGA-BH-A0EA-01A-11W-A12T-09",
      "TCGA-BH-A0EA-01A-11D-A10X-02",
      "TCGA-BH-A0EA-01A-11D-A10Y-09",
      "TCGA-BH-A0EA-01A-11D-A314-09",
      "TCGA-BH-A0EA-01A-11R-A115-07",
      "TCGA-BH-A0EA-01A-11D-A112-05",
      "TCGA-BH-A0EA-10A-01D-A110-09",
      "TCGA-BH-A0EA-10A-01D-A113-01",
      "TCGA-BH-A0EA-10A-01W-A12U-09",
      "TCGA-BH-A0EA-10A-01D-A10Z-02",
      "TCGA-BH-A0EA-10A-01D-A113-01",
      "TCGA-BH-A0EA-10A-01D-A110-09",
      "TCGA-BH-A0EA-10A-01W-A12U-09",
      "TCGA-BH-A0EA-10A-01D-A10Z-02"
    ],
    "updated_datetime": "2016-05-02T14:37:43.619198-05:00",
    "submitter_analyte_ids": [
      "TCGA-BH-A0EA-01A-11R",
      "TCGA-BH-A0EA-01A-11D",
      "TCGA-BH-A0EA-01A-11W",
      "TCGA-BH-A0EA-10A-01W",
      "TCGA-BH-A0EA-10A-01D"
    ],
    "analyte_ids": [
      "30cb470f-66d4-4085-8c30-83a42e8453d4",
      "66ed0f86-5ca5-4dec-ba76-7ee4dcf31831",
      "f19f408a-815f-43d9-8032-e9482b796371",
      "69ddc092-88a0-4839-a2bb-9f1c9e760409",
      "fe678556-acf4-4bde-a95e-860bb0150a95"
    ],
    "submitter_id": "TCGA-BH-A0EA",
    "case_id": "1f601832-eee3-48fb-acf5-80c4a454f26e",
    "state": null,
    "aliquot_ids": [
      "bcb7fc6d-60a0-48b7-aa81-14c0dda72d76",
      "97c64d6a-7dce-4d0f-9cb3-b3e4eb4719c5",
      "edad5bd3-efe0-4c5f-b05c-2c0c2951c45a",
      "bcb7fc6d-60a0-48b7-aa81-14c0dda72d76",
      "ca71ca96-cbb7-4eab-9487-251dda34e107",
      "97c64d6a-7dce-4d0f-9cb3-b3e4eb4719c5",
      "eef9dce1-6ba6-432b-bbe2-53c7dbe64fe7",
      "42d050e4-e8ee-4442-b9c0-0ee14706b138",
      "561b8777-801a-49ed-a306-e7dafeb044b6",
      "262715e1-835c-4f16-8ee7-6900e26f7cf5",
      "edad5bd3-efe0-4c5f-b05c-2c0c2951c45a",
      "262715e1-835c-4f16-8ee7-6900e26f7cf5",
      "561b8777-801a-49ed-a306-e7dafeb044b6",
      "eef9dce1-6ba6-432b-bbe2-53c7dbe64fe7",
      "ca71ca96-cbb7-4eab-9487-251dda34e107",
      "42d050e4-e8ee-4442-b9c0-0ee14706b138",
      "cfbd5476-e83a-401d-9f9a-639c73a0e35b",
      "2beb34c4-d493-4a73-b21e-de77d43251ff",
      "b1a3739d-d554-4202-b96f-f25a444e2042",
      "cde982b7-3b0a-49eb-8710-a599cb0e44c1",
      "2beb34c4-d493-4a73-b21e-de77d43251ff",
      "cfbd5476-e83a-401d-9f9a-639c73a0e35b",
      "b1a3739d-d554-4202-b96f-f25a444e2042",
      "cde982b7-3b0a-49eb-8710-a599cb0e44c1"
    ],
    "slide_ids": [
      "90154ea1-6b76-4445-870e-d531d6fa1239",
      "a0826f0d-986a-491b-8c6f-b34f8929f3ee"
    ],
    "submitter_sample_ids": [
      "TCGA-BH-A0EA-01A",
      "TCGA-BH-A0EA-10A"
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
[{"release_date": "2018-07-23", "version": "1", "uuid": "1dd28069-5777-4ff9-bd2b-d1ba68e88b06", "file_change": "released", "data_release": "13.0"}]
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
        "sample_ids": [
          "1d014bf1-95ae-42e3-ae39-97ff4841d8ca",
          "6b685bfc-651b-48d1-8e68-32c8096ea205"
        ],
        "portion_ids": [
          "c061217a-266a-496d-8a96-3489191afa87",
          "0d3a6a58-0e00-4889-bc73-5ddb5a387738",
          "e858ee92-0438-48e9-a70d-80ef2c0ad539"
        ],
        "submitter_portion_ids": [
          "TCGA-66-2770-01A-21-2193-20",
          "TCGA-66-2770-01A-01",
          "TCGA-66-2770-11A-01"
        ],
        "created_datetime": null,
        "submitter_aliquot_ids": [
          "TCGA-66-2770-01A-01D-1522-08",
          "TCGA-66-2770-01A-01D-0848-05",
          "TCGA-66-2770-01A-01W-0879-09",
          "TCGA-66-2770-11A-01W-0878-08",
          "TCGA-66-2770-01A-01R-0849-01",
          "TCGA-66-2770-01A-01W-0877-08",
          "TCGA-66-2770-01A-01D-0846-06",
          "TCGA-66-2770-11A-01W-0880-09",
          "TCGA-66-2770-01A-01D-0964-09",
          "TCGA-66-2770-11A-01D-0846-06",
          "TCGA-66-2770-01A-01D-0845-04",
          "TCGA-66-2770-01A-01W-0881-10",
          "TCGA-66-2770-11A-01D-0963-08",
          "TCGA-66-2770-11A-01D-0844-01",
          "TCGA-66-2770-01A-01R-0851-07",
          "TCGA-66-2770-11A-01W-0882-10",
          "TCGA-66-2770-11A-01D-1522-08",
          "TCGA-66-2770-01A-01T-1557-13",
          "TCGA-66-2770-01A-01D-0847-02",
          "TCGA-66-2770-01A-01D-0844-01",
          "TCGA-66-2770-11A-01D-0847-02",
          "TCGA-66-2770-11A-01D-0964-09",
          "TCGA-66-2770-01A-01D-0963-08",
          "TCGA-66-2770-01A-01R-0850-03",
          "TCGA-66-2770-11A-01D-0845-04",
          "TCGA-66-2770-01A-01T-0852-07"
        ],
        "updated_datetime": "2016-05-02T15:57:03.730994-05:00",
        "submitter_analyte_ids": [
          "TCGA-66-2770-01A-01D",
          "TCGA-66-2770-11A-01W",
          "TCGA-66-2770-01A-01T",
          "TCGA-66-2770-01A-01W",
          "TCGA-66-2770-01A-01R",
          "TCGA-66-2770-11A-01D"
        ],
        "analyte_ids": [
          "385807d3-78de-4558-8d93-702d93fc835a",
          "247acc7a-b4f5-47e9-86da-5ea9b04ad444",
          "151b8cb9-6b0a-4db9-9b0e-62aa501b35d9",
          "e549aebd-4dda-4ea8-8ccf-56c03bc8b2be",
          "631ad4eb-845a-4e70-96ad-4b40157218a8",
          "9a75640e-09d4-42b7-8cb4-75d62b39e98a"
        ],
        "submitter_id": "TCGA-66-2770",
        "case_id": "f1b357e4-d67a-42c9-b0b7-12f69fa3da58",
        "state": null,
        "aliquot_ids": [
          "a2d10f8e-6b27-4df0-bd25-ac24992d0bb4",
          "8c1c733a-abed-468f-b4d0-d1ac34ba6d8b",
          "cad8d384-3b7a-4f70-89c2-5584ae75c5eb",
          "42e774cf-3c4a-4efd-9665-378cb6b4afac",
          "3755168b-f5da-422d-847a-566cb112a8d7",
          "cae4d249-ba67-4316-8761-7e71e3813182",
          "aa6e700c-ce01-4cc9-87de-8bf615a8aa1a",
          "ad5c4069-e616-4ab4-9b03-b196f9189b20",
          "07c26ea4-0584-4cb0-8e5a-d057b8fe6c14",
          "f95c2cb5-d20a-4f1f-8f2a-95a2d37fbdc4",
          "817bf327-e583-4704-b294-c3645dcc4adf",
          "2246cb75-38bd-491f-b6ee-99f4781f2564",
          "a81b9090-626d-492d-9baf-7fa3ef70111c",
          "5cd6f026-894e-45f6-bc59-d6f056e63846",
          "e417903d-ab76-44f0-aae9-3a91fa9a8d3c",
          "1d809a56-31ca-49d8-a57b-e773236b24de",
          "df60a743-ef4b-43ea-bc5a-4d75e8befb8a",
          "871350e2-958f-401c-ae86-6bc880a01942",
          "3dc4207d-5671-4c3d-b75a-d39ef69b564c",
          "69b77cc0-d00a-4ea3-9b39-3e3019d9e292",
          "3d035ee8-9523-4771-8738-c8a5a2f91403",
          "775e46bd-e56f-40fa-9891-aaedc1d49395",
          "d1c60049-922a-42d4-bd7e-8cf4ace47f05",
          "5220a53f-f3fc-476c-aa72-65a038eb2fd8",
          "b7e44e6e-ccf9-4b75-a258-159912ab51ca",
          "42750622-28d7-4d32-9262-b139fe77bc01"
        ],
        "slide_ids": [
          "a10196d2-7a81-4e1e-a9a7-62d123c30875",
          "72edc1ba-916d-42a2-9f22-6254c6e54c5c",
          "ff15eeb9-550e-4c78-90cc-a6cce8ccc3df",
          "71ccfb52-169d-4176-94d6-fff5b75f853d"
        ],
        "submitter_sample_ids": [
          "TCGA-66-2770-11A",
          "TCGA-66-2770-01A"
        ]
      },
      {
        "sample_ids": [
          "06889714-2a40-4248-98ee-f690b301e36a",
          "9f43a0c6-ea19-4021-b0ed-026f33ce1c33"
        ],
        "portion_ids": [
          "3a001d28-7cf9-4c61-b155-73938aebaa25",
          "79554cfd-e853-481e-8e37-1e296034094e"
        ],
        "submitter_portion_ids": [
          "TCGA-02-0075-01A-01",
          "TCGA-02-0075-10A-01"
        ],
        "created_datetime": null,
        "submitter_aliquot_ids": [
          "TCGA-02-0075-01A-01W-0204-02",
          "TCGA-02-0075-01A-01R-0194-03",
          "TCGA-02-0075-01A-01D-0198-02",
          "TCGA-02-0075-01A-01R-0202-01",
          "TCGA-02-0075-10A-01W-0207-09",
          "TCGA-02-0075-01A-01R-0676-04",
          "TCGA-02-0075-10A-01D-0198-02",
          "TCGA-02-0075-10A-01D-0197-06",
          "TCGA-02-0075-10A-01D-0193-01",
          "TCGA-02-0075-01A-01W-0207-09",
          "TCGA-02-0075-01A-01W-0206-08",
          "TCGA-02-0075-01A-01D-0193-01",
          "TCGA-02-0075-10A-01W-0205-10",
          "TCGA-02-0075-01A-01R-0201-02",
          "TCGA-02-0075-10A-01W-0204-02",
          "TCGA-02-0075-01A-01D-0199-05",
          "TCGA-02-0075-10A-01W-0206-08",
          "TCGA-02-0075-01A-01D-0196-04",
          "TCGA-02-0075-01A-01T-0195-07",
          "TCGA-02-0075-10A-01D-0196-04",
          "TCGA-02-0075-01A-01D-0197-06",
          "TCGA-02-0075-01A-01D-0888-01",
          "TCGA-02-0075-01A-01R-0195-07",
          "TCGA-02-0075-01A-01W-0205-10"
        ],
        "updated_datetime": "2016-05-02T15:00:01.972331-05:00",
        "submitter_analyte_ids": [
          "TCGA-02-0075-01A-01R",
          "TCGA-02-0075-10A-01D",
          "TCGA-02-0075-01A-01W",
          "TCGA-02-0075-01A-01T",
          "TCGA-02-0075-01A-01D",
          "TCGA-02-0075-10A-01W"
        ],
        "analyte_ids": [
          "fec22de0-a2b9-45df-9854-1ebe76cee84e",
          "b4d11c50-61f1-4d4a-815f-1c0413018d7f",
          "c48673d0-a38d-44e1-8cfd-e91cb23ea2d5",
          "24f1852c-999a-4ea8-917c-fcfd683e2aca",
          "aa431260-a0fc-4924-80ce-61cab8b5e83e",
          "11f21140-d761-44ca-a9b2-b24099df3b15"
        ],
        "submitter_id": "TCGA-02-0075",
        "case_id": "b196f82b-ef3f-4e05-99f7-da5df65e691e",
        "state": null,
        "aliquot_ids": [
          "75531fe0-101e-4220-bd47-98892c90ee70",
          "e5ea38d4-f47c-4c8a-8bab-13631e0a9a7b",
          "d48b7c2c-daac-4496-af8f-1f45ca43f627",
          "bbba08fc-2514-4e15-afb7-41eecc7e876f",
          "0685b37f-a47c-4222-a846-bf9f3c000de3",
          "683986da-3cee-446d-9b7a-83bef25815c9",
          "e6ffdb20-a1be-4664-bcd3-cc7a4de6f40b",
          "5d1f25c0-9e1a-41ad-9735-134f39dbf70e",
          "528b40b9-246f-4ba3-8209-777136638e62",
          "33131479-5d69-4262-a549-ba8864320f3b",
          "5c7822fc-cf4f-4f62-8482-7c0ce1b7ab9a",
          "b95e7659-e3a4-4e96-b98c-f67d26b85322",
          "30c84aca-f9db-4e07-ac34-1a92b1652ca1",
          "d5e3b5cc-06e0-4294-9d3c-8f3b63acae3d",
          "b14b3d09-3a7f-41a6-81df-2757efa67906",
          "513040e2-dc29-4e2c-86fb-57371eede17a",
          "21c3be1b-7c1e-4864-99d1-486cfe5d8f1d",
          "5e28e5dc-6dfa-44a9-8793-9134cb4cdda5",
          "b8c25892-4773-428f-a02c-f930931268e8",
          "266d5260-08e4-4cec-87f3-ca415bd98575",
          "8859a3ae-f85d-4ef2-830b-80f42f98d53e",
          "ac018a8c-a6e2-4291-a4bf-a330ae9c441e",
          "4b022f7f-7549-4d97-9d41-4e5f2e9ec74c",
          "caad3dfa-74a9-4ecc-95c1-86f6fbfd4ab5"
        ],
        "slide_ids": [
          "39f547cd-5dc3-4bf4-99ea-073bb161c23c",
          "5f096267-0cc2-4cc5-a206-7357159633d7"
        ],
        "submitter_sample_ids": [
          "TCGA-02-0075-10A",
          "TCGA-02-0075-01A"
        ]
      },
      {
        "sample_ids": [
          "ba08195b-31cf-4bb1-a470-23740225c99d",
          "929889c4-e474-4104-b69b-fac7e414a59e"
        ],
        "portion_ids": [
          "48a36eb4-79fb-45e7-8bb1-0fa1d5fcda2c",
          "1de5e67a-ac3f-4c18-92c4-27ba1868c7ac",
          "e09fc5e7-e8d2-4bf9-b12b-17b22e0387e4"
        ],
        "submitter_portion_ids": [
          "TCGA-EJ-A8FU-10A-01",
          "TCGA-EJ-A8FU-01A-21-A43L-20",
          "TCGA-EJ-A8FU-01A-11"
        ],
        "created_datetime": null,
        "submitter_aliquot_ids": [
          "TCGA-EJ-A8FU-01A-11R-A36B-13",
          "TCGA-EJ-A8FU-01A-11R-A36G-07",
          "TCGA-EJ-A8FU-01A-11D-A363-01",
          "TCGA-EJ-A8FU-10A-01D-A361-01",
          "TCGA-EJ-A8FU-10A-01D-A362-08",
          "TCGA-EJ-A8FU-01A-11W-A447-08",
          "TCGA-EJ-A8FU-01A-11D-A365-05",
          "TCGA-EJ-A8FU-01A-11D-A364-08",
          "TCGA-EJ-A8FU-10A-01W-A446-08"
        ],
        "updated_datetime": "2016-05-02T15:57:04.948573-05:00",
        "submitter_analyte_ids": [
          "TCGA-EJ-A8FU-01A-11W",
          "TCGA-EJ-A8FU-01A-11D",
          "TCGA-EJ-A8FU-01A-11R",
          "TCGA-EJ-A8FU-10A-01W",
          "TCGA-EJ-A8FU-10A-01D"
        ],
        "analyte_ids": [
          "2d4e4925-6ac8-498f-882b-4bbf319f6b7b",
          "8d09b982-1256-4674-b383-d6ca4b4bb3c8",
          "c74495d9-63bf-4ac0-b10e-04b3b06103c1",
          "b9884d98-af57-4901-8b9d-4fdbf73d2c5a",
          "2f16ac02-13bf-44fd-bbd7-658c1c384928"
        ],
        "submitter_id": "TCGA-EJ-A8FU",
        "case_id": "23e56e08-e11d-4e83-88a8-1254675b3af8",
        "state": null,
        "aliquot_ids": [
          "e77da017-5dc6-4e32-9568-755e4ee9b533",
          "c9b286d1-d500-4bb3-bb3d-5bf40b1b1265",
          "b7867d52-7987-46d4-a595-0ff5b5375a58",
          "5586ad35-94b7-459e-8982-8e7fb25697a1",
          "162a63f7-594f-4669-a06d-b4899c7fe86a",
          "b8b1ab44-ee6e-4ac5-9efd-d5bd07e67b9c",
          "7adcdf73-3ad3-4da7-ab27-2888f1d4f53a",
          "eb498e52-3eae-402f-8cac-ec930f8d938d",
          "293f781c-c2c7-479b-b1a6-5f951a2c5e5a"
        ],
        "slide_ids": [
          "454a95d5-d084-4f36-b1f1-32c6c23ab46e"
        ],
        "submitter_sample_ids": [
          "TCGA-EJ-A8FU-01A",
          "TCGA-EJ-A8FU-10A"
        ]
      },
      {
        "sample_ids": [
          "d43f0112-fe59-4842-9fda-1189e5fb7248",
          "213cbbe5-c382-47a1-b936-bf40c2c99091"
        ],
        "portion_ids": [
          "26441aae-22e5-4e69-b3f5-34ccde356c93",
          "60d7a93c-0634-438e-a72a-ce63630bb890",
          "246a8f01-7ef2-4737-a984-49aa0b41c089"
        ],
        "submitter_portion_ids": [
          "TCGA-F2-6879-10A-01",
          "TCGA-F2-6879-01A-21-A39M-20",
          "TCGA-F2-6879-01A-11"
        ],
        "created_datetime": "2016-05-02T16:23:44.347995-05:00",
        "submitter_aliquot_ids": [
          "TCGA-F2-6879-01A-11R-2155-13",
          "TCGA-F2-6879-10A-01D-2153-01",
          "TCGA-F2-6879-10A-01D-2152-26",
          "TCGA-F2-6879-01A-11D-2157-05",
          "TCGA-F2-6879-10A-01D-2154-08",
          "TCGA-F2-6879-01A-11D-A45X-08",
          "TCGA-F2-6879-01A-11D-2154-08",
          "TCGA-F2-6879-01A-11W-2179-08",
          "TCGA-F2-6879-01A-11D-2153-01",
          "TCGA-F2-6879-01A-11R-2156-07",
          "TCGA-F2-6879-01A-11D-2152-26",
          "TCGA-F2-6879-10A-01D-A45X-08",
          "TCGA-F2-6879-10A-01W-2179-08",
          "TCGA-F2-6879-01A-01D-YYYY-23"
        ],
        "updated_datetime": "2016-05-02T16:23:44.347995-05:00",
        "submitter_analyte_ids": [
          "TCGA-F2-6879-10A-01D",
          "TCGA-F2-6879-01A-11R",
          "TCGA-F2-6879-10A-01W",
          "TCGA-F2-6879-01A-11W",
          "TCGA-F2-6879-01A-11D"
        ],
        "analyte_ids": [
          "e87dde8d-3bf5-42d8-9a77-620d5c4943e0",
          "30ade77d-996b-4031-93ab-6b341d49eb0a",
          "1d94bd70-6621-4a94-8102-d673663e6665",
          "ea65d92e-1597-410d-84d8-abb2a6235b3e",
          "79697034-1cec-4d92-8195-8a35258ab477"
        ],
        "submitter_id": "TCGA-F2-6879",
        "case_id": "8d9bd437-8b4b-4da5-87ba-6b5790f05022",
        "state": null,
        "aliquot_ids": [
          "e7533585-b062-4d74-b511-05dc806a1357",
          "e107952a-cc2b-4410-b0f9-62e7115430a0",
          "61f1c8b1-986a-485a-9d96-4e4285b6425a",
          "c043e276-fece-4cb9-a848-a0b16e6099b6",
          "e5d110e1-63ad-49ce-b9b7-22bbd7ef8a88",
          "7accb08d-acdb-46bc-bf7f-b9f678193115",
          "a52cd04b-41d6-40db-b050-00ef3a143f7e",
          "207fcf5e-c422-4333-9ec2-5dab38d240c7",
          "5ddd3f83-28a8-4b7f-9aec-203a3c2efbe5",
          "ccd4dd70-c0e4-42cf-870e-33d1013b201a",
          "e12314fe-f16a-4d85-95b4-e712ede450f6",
          "695461e3-283c-4b5b-9325-6b2588b67fd8",
          "8481be1e-0993-487d-8d73-b0eb72b304ee",
          "d7200791-4f1c-418f-8744-91b793486d9f"
        ],
        "slide_ids": [
          "bcbcc947-cab1-4400-aebc-1d9e251a3ce8",
          "cae8d0b9-3605-40af-bf99-7c23df8110a9"
        ],
        "submitter_sample_ids": [
          "TCGA-F2-6879-10A",
          "TCGA-F2-6879-01A"
        ]
      },
      {
        "sample_ids": [
          "3a66b5bd-7037-463c-9f8d-2ba3de9d5571",
          "84f603d6-9f71-48fb-b2e3-190424407452"
        ],
        "portion_ids": [
          "fe90de9f-8ee3-4d55-834f-a90538958cb7",
          "7a0042fd-07f0-4894-adb0-03cebce8aa02"
        ],
        "submitter_portion_ids": [
          "TCGA-VQ-A922-01A-11",
          "TCGA-VQ-A922-10A-01"
        ],
        "created_datetime": "2016-05-02T16:26:23.121974-05:00",
        "submitter_aliquot_ids": [
          "TCGA-VQ-A922-10A-01D-A412-01",
          "TCGA-VQ-A922-01A-11D-A40Z-01",
          "TCGA-VQ-A922-10A-01D-A413-08",
          "TCGA-VQ-A922-01A-01D-YYYY-23",
          "TCGA-VQ-A922-01A-11R-A414-31",
          "TCGA-VQ-A922-01A-11D-A410-08",
          "TCGA-VQ-A922-01A-11R-A415-13",
          "TCGA-VQ-A922-01A-11D-A411-05"
        ],
        "updated_datetime": "2016-05-02T16:26:23.121974-05:00",
        "submitter_analyte_ids": [
          "TCGA-VQ-A922-01A-11R",
          "TCGA-VQ-A922-10A-01D",
          "TCGA-VQ-A922-01A-11D"
        ],
        "analyte_ids": [
          "15bec495-04c7-412b-ad69-26b1f9274ccf",
          "26a24673-04a1-4837-b888-702b0578aef2",
          "2c0ecd67-b9ff-4e60-8d2f-7744c79a13aa"
        ],
        "submitter_id": "TCGA-VQ-A922",
        "case_id": "8bd783a3-d6c9-4c87-a2a1-09f903b9c7ca",
        "state": null,
        "aliquot_ids": [
          "58a121b4-265c-44ae-b6a9-79d087ee8b34",
          "76fbba49-0123-4524-89aa-a1818c5507cb",
          "0b0805bb-edaa-400f-ae9f-effed3dbb605",
          "3370d626-d572-4d13-9cd3-1823a5df3d34",
          "60934993-a9df-4389-b64d-da6844ef22df",
          "243f24ba-bb0f-44e0-bcb1-69a97b395981",
          "6cae9f2a-1c6c-4645-98b6-20719aec1413",
          "44d020d1-c516-4a15-94e8-bcf0cb9c2683"
        ],
        "slide_ids": [
          "0ff02899-57f8-419e-8872-c6ede53f4d3c"
        ],
        "submitter_sample_ids": [
          "TCGA-VQ-A922-10A",
          "TCGA-VQ-A922-01A"
        ]
      },
      {
        "sample_ids": [
          "5bb5bd60-cf47-413b-88fa-f14977e24035",
          "82fcf670-1646-4a28-9578-f7e5b2f426e5",
          "3b87fed0-cfbd-4ee3-b71d-ab595853e836"
        ],
        "portion_ids": [
          "18bf160e-702a-464a-9920-f115024b5484",
          "10a9c093-009d-4bc0-a344-2afd3f0f9b9f",
          "8ebd06e1-5eda-47ec-8888-61965ecf005e"
        ],
        "submitter_portion_ids": [
          "TCGA-HU-8243-11A-01",
          "TCGA-HU-8243-01A-11",
          "TCGA-HU-8243-10A-01"
        ],
        "created_datetime": "2016-05-02T16:17:09.754748-05:00",
        "submitter_aliquot_ids": [
          "TCGA-HU-8243-01A-01D-YYYY-23",
          "TCGA-HU-8243-01A-11D-2340-08",
          "TCGA-HU-8243-01A-11D-2338-01",
          "TCGA-HU-8243-01A-11D-2342-05",
          "TCGA-HU-8243-11A-01D-2338-01",
          "TCGA-HU-8243-11A-01D-2340-08",
          "TCGA-HU-8243-10A-01D-2339-01",
          "TCGA-HU-8243-01A-11R-2343-13",
          "TCGA-HU-8243-10A-01D-2341-08"
        ],
        "updated_datetime": "2016-05-02T16:17:09.754748-05:00",
        "submitter_analyte_ids": [
          "TCGA-HU-8243-11A-01D",
          "TCGA-HU-8243-10A-01D",
          "TCGA-HU-8243-01A-11R",
          "TCGA-HU-8243-01A-11D"
        ],
        "analyte_ids": [
          "89c9094d-5cf6-4c7d-ad24-41b7ad9427cc",
          "2c413e60-0122-426b-afb3-ae94810e2513",
          "57d41760-0fed-49d2-8606-48231cb244ea",
          "37ed51fd-b540-408e-8bd6-4447ae4aa84a"
        ],
        "submitter_id": "TCGA-HU-8243",
        "case_id": "77a8eab6-f6a1-4739-9031-75ead40d68cb",
        "state": null,
        "aliquot_ids": [
          "ace3edd6-14a9-42cc-84f3-6127237f2913",
          "a711abd1-f1c2-4e42-8b66-79b4514ac1c4",
          "6af7ba34-58f7-4472-8c7e-89fc91ad5ac1",
          "558ff67a-a584-46f8-9089-8f4a08015294",
          "71c0a224-5953-4b59-a49c-b7aa1e959f1e",
          "a460c222-bcac-4959-961f-4dbd73e1ce13",
          "6e5789d7-4988-457a-86eb-e618c7ab06eb",
          "ff31f56b-398c-45ee-b122-f10027774527",
          "9635cfd4-3d26-4fc6-846c-fd74d5b60098"
        ],
        "slide_ids": [
          "60b7c6b8-594a-40c3-9341-a0902e4e6938",
          "e55e00a0-2048-404a-b83a-f34106468694"
        ],
        "submitter_sample_ids": [
          "TCGA-HU-8243-10A",
          "TCGA-HU-8243-01A",
          "TCGA-HU-8243-11A"
        ]
      },
      {
        "sample_ids": [
          "2f5cc9c9-31a9-5eb3-952a-b21e7cef50ca",
          "4f3f4fc8-4465-5230-83ec-c0ef6aceb2ea"
        ],
        "updated_datetime": "2016-05-25T19:12:45.610324-05:00",
        "submitter_aliquot_ids": [
          "TARGET-30-PAUXFZ-01A-01D",
          "TARGET-30-PAUXFZ-10A-01D"
        ],
        "submitter_id": "TARGET-30-PAUXFZ",
        "case_id": "a7ccef7c-14c0-5232-b647-58b4a54fb343",
        "aliquot_ids": [
          "9e1e30a8-7607-5b7e-b33c-9a6c5828d5fb",
          "c56898f9-c394-516a-bdbb-bf32a5af9d3f"
        ],
        "submitter_sample_ids": [
          "TARGET-30-PAUXFZ-01A",
          "TARGET-30-PAUXFZ-10A"
        ]
      },
      {
        "sample_ids": [
          "c1bcb8d1-e13d-4af4-93f4-02d5f7f616a2",
          "52fcf737-cdcc-43ea-b33c-4018039b42dd"
        ],
        "portion_ids": [
          "e0e97a05-656a-468e-8418-0d08c38e76ab",
          "3e2a0eab-7d89-4f3c-9c0e-8942e53d3c45"
        ],
        "submitter_portion_ids": [
          "TCGA-KK-A8I9-01A-11",
          "TCGA-KK-A8I9-11A-11"
        ],
        "created_datetime": null,
        "submitter_aliquot_ids": [
          "TCGA-KK-A8I9-11A-11D-A361-01",
          "TCGA-KK-A8I9-11A-11D-A362-08",
          "TCGA-KK-A8I9-11A-11W-A446-08",
          "TCGA-KK-A8I9-01A-11R-A36G-07",
          "TCGA-KK-A8I9-11A-11D-A40C-01",
          "TCGA-KK-A8I9-01A-11D-A363-01",
          "TCGA-KK-A8I9-01A-11W-A447-08",
          "TCGA-KK-A8I9-01A-11D-A365-05",
          "TCGA-KK-A8I9-01A-11D-A364-08",
          "TCGA-KK-A8I9-01A-11R-A36B-13"
        ],
        "updated_datetime": "2016-05-02T15:57:29.451686-05:00",
        "submitter_analyte_ids": [
          "TCGA-KK-A8I9-11A-11W",
          "TCGA-KK-A8I9-01A-11R",
          "TCGA-KK-A8I9-11A-11D",
          "TCGA-KK-A8I9-01A-11W",
          "TCGA-KK-A8I9-01A-11D"
        ],
        "analyte_ids": [
          "ddec19cb-5e4c-4151-8b6d-741044abff1e",
          "96c5b539-8eb7-4156-81d0-7b7fecd68900",
          "ced38a45-7610-49d4-8bf9-d53a1fc2d489",
          "476f5deb-1b3f-4a35-8a31-f27763ba8d8a",
          "c284f2af-1e9b-40cc-8936-b61cfd251d62"
        ],
        "submitter_id": "TCGA-KK-A8I9",
        "case_id": "261c3d74-706e-4751-bd15-8f3c1a402ff0",
        "state": null,
        "aliquot_ids": [
          "4f76de2d-e07a-402b-9818-7f04d3704a43",
          "96802a73-b1db-47d7-8f5f-4504f3ece5ad",
          "f376fc45-370a-4d96-833b-9a1322e32a42",
          "d3e88dd3-66d7-40d4-978a-4ddab868373a",
          "06f1d087-75c9-4da8-8339-80aff3bfaa12",
          "50b1e243-b45a-42a1-8692-b7ae5d51250f",
          "0f1c00d3-f3dc-4d2b-bd8a-ecc31e4f4089",
          "986a3ed6-ba56-4025-a2bd-9909648e703a",
          "bebc84b6-9179-420b-8207-858b999e8c0c",
          "239d5e7e-5fb5-4df3-ae6b-a5a06ee296ae"
        ],
        "slide_ids": [
          "1e174ca5-9298-41b6-a705-728f111a3e7b",
          "a3e31324-9e06-4799-85b4-4f6236848009"
        ],
        "submitter_sample_ids": [
          "TCGA-KK-A8I9-11A",
          "TCGA-KK-A8I9-01A"
        ]
      },
      {
        "sample_ids": [
          "d43f727a-96d6-40b8-86ae-7a3e0aa46853",
          "b8329a6d-a87b-47f4-ad00-9e979e62647b"
        ],
        "portion_ids": [
          "8960ddcc-0950-4d6e-a557-8727b652c93b",
          "e36bfd07-c911-4a98-8424-e58e5e9aaa68"
        ],
        "submitter_portion_ids": [
          "TCGA-QR-A70H-10A-01",
          "TCGA-QR-A70H-01A-12"
        ],
        "created_datetime": null,
        "submitter_aliquot_ids": [
          "TCGA-QR-A70H-01A-12R-A35K-07",
          "TCGA-QR-A70H-01A-12R-A35M-13",
          "TCGA-QR-A70H-01A-12D-A35E-05",
          "TCGA-QR-A70H-10A-01D-A35A-01",
          "TCGA-QR-A70H-01A-12D-A35C-01",
          "TCGA-QR-A70H-01A-12W-A43Z-08",
          "TCGA-QR-A70H-10A-01D-A35B-08",
          "TCGA-QR-A70H-10A-01W-A441-08",
          "TCGA-QR-A70H-01A-12D-A35D-08"
        ],
        "updated_datetime": "2016-05-02T15:37:31.996088-05:00",
        "submitter_analyte_ids": [
          "TCGA-QR-A70H-10A-01D",
          "TCGA-QR-A70H-10A-01W",
          "TCGA-QR-A70H-01A-12D",
          "TCGA-QR-A70H-01A-12W",
          "TCGA-QR-A70H-01A-12R"
        ],
        "analyte_ids": [
          "c4a41555-dd45-4e10-a3be-50d49a1121a3",
          "957e01f6-eb3f-446e-9f45-b50c66337e2d",
          "1acde950-2e0c-4586-852b-b4ac4e1ea4a4",
          "67c033c0-9fe8-4004-967e-c605e1890f4d",
          "b0873010-5d60-4691-b700-e172950f1d7c"
        ],
        "submitter_id": "TCGA-QR-A70H",
        "case_id": "13b41b15-a785-4ab7-b864-ffff6d35dd45",
        "state": null,
        "aliquot_ids": [
          "d9120f00-7f10-49d5-ae84-6177e9424c7c",
          "31c6fa50-200a-46c1-a546-61b52592fd8f",
          "ab50f38c-2e7d-4d75-a216-27aeaa4d9305",
          "382d5e31-6c66-4df3-a695-6b8c29cfc681",
          "51d1fb14-c918-4439-b816-ef6cd3253c64",
          "f586d8d5-d0c6-4979-aaa7-10217a88fa4c",
          "2f9a60eb-602e-44bb-bc57-87e20d946f76",
          "fbafc85e-deff-46cd-a40f-479b9dc92a60",
          "cacbc8a6-0eb0-4277-931f-d0075c9b1de9"
        ],
        "slide_ids": [
          "2310e34c-0ea5-4876-9f87-bad0b7a44513"
        ],
        "submitter_sample_ids": [
          "TCGA-QR-A70H-01A",
          "TCGA-QR-A70H-10A"
        ]
      },
      {
        "sample_ids": [
          "19dee039-9c98-4d4a-8baf-eea1b6dda8eb",
          "fdf1e501-f34f-450c-9a5c-611157079a86"
        ],
        "portion_ids": [
          "10b6ccb4-3637-4769-8988-417c0306eaef",
          "92f8cd48-451d-4ed6-8e60-b15aa93d2c09",
          "d0d55efa-c91d-45de-92bf-cf6f0d263b21"
        ],
        "submitter_portion_ids": [
          "TCGA-BJ-A18Z-01A-21",
          "TCGA-BJ-A18Z-01A-11-A21L-20",
          "TCGA-BJ-A18Z-10A-01"
        ],
        "created_datetime": null,
        "submitter_aliquot_ids": [
          "TCGA-BJ-A18Z-01A-21D-A13U-02",
          "TCGA-BJ-A18Z-10A-01D-A13V-01",
          "TCGA-BJ-A18Z-01A-21R-A13Y-07",
          "TCGA-BJ-A18Z-01A-21W-A14T-08",
          "TCGA-BJ-A18Z-01A-21D-A13Z-05",
          "TCGA-BJ-A18Z-01A-21D-A37T-08",
          "TCGA-BJ-A18Z-10A-01D-A13W-08",
          "TCGA-BJ-A18Z-01A-21R-A13X-13",
          "TCGA-BJ-A18Z-01A-21D-A13W-08",
          "TCGA-BJ-A18Z-10A-01D-A13U-02",
          "TCGA-BJ-A18Z-10A-01W-A14T-08",
          "TCGA-BJ-A18Z-01A-21D-A13V-01"
        ],
        "updated_datetime": "2016-05-02T16:18:19.199189-05:00",
        "submitter_analyte_ids": [
          "TCGA-BJ-A18Z-01A-21W",
          "TCGA-BJ-A18Z-01A-21D",
          "TCGA-BJ-A18Z-01A-21R",
          "TCGA-BJ-A18Z-10A-01D",
          "TCGA-BJ-A18Z-10A-01W"
        ],
        "analyte_ids": [
          "119ebfa1-75b2-4f24-816a-4e9a5061f6b5",
          "f86759fd-ecc5-4f42-b5fe-b9f079d23968",
          "39691042-bd28-40ed-b66b-26414ecf1ba0",
          "76ea5056-d7fa-49fb-94bf-11171ca7c100",
          "71a822c9-b510-4a4c-8c30-18b8083acc2d"
        ],
        "submitter_id": "TCGA-BJ-A18Z",
        "case_id": "0d497faf-2c1c-4173-a5fe-770cca73323c",
        "state": null,
        "aliquot_ids": [
          "fa580596-e70f-4ed0-85a2-6fb594ca679a",
          "776cb4b1-8efd-4ea2-b53f-9dff7dd94b10",
          "85a7922f-0327-437c-bdf5-1bb67a1e932f",
          "6d532180-0175-4610-8bfa-cca3a7c3697a",
          "b5977e73-49d8-4e99-9e97-993cc44dad17",
          "918793fa-b35e-4745-ac75-4d1c868089f8",
          "ba9479a1-929f-4e4e-8bf5-e23cb280dfcf",
          "e9776ff5-69b9-4669-ab33-e4bb030461ec",
          "8ba98907-ab03-4c9e-a900-e31aa16ff810",
          "35e18649-183e-4223-b2f6-d812bdd9becd",
          "4aa17671-4420-4989-a6dd-379250f4aeda",
          "815c53c3-8add-4612-b93c-3ed4bfa530aa"
        ],
        "slide_ids": [
          "7c5b5c77-9fbc-4b48-81f5-48b5ede7c436"
        ],
        "submitter_sample_ids": [
          "TCGA-BJ-A18Z-01A",
          "TCGA-BJ-A18Z-10A"
        ]
      }
    ],
    "pagination": {
      "count": 10,
      "sort": "",
      "from": 0,
      "page": 1,
      "total": 6340,
      "pages": 634,
      "size": 10
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
submitter_id
TCGA-RC-A6M6
TCGA-B6-A0RV
TCGA-MB-A5Y8
TCGA-BQ-5876
TCGA-Z6-A9VB
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
				<submitter_id>TCGA-MQ-A4LV</submitter_id>
			</item>
			<item>
				<submitter_id>TCGA-N9-A4Q1</submitter_id>
			</item>
			<item>
				<submitter_id>TCGA-78-7154</submitter_id>
			</item>
			<item>
				<submitter_id>TCGA-S7-A7WX</submitter_id>
			</item>
			<item>
				<submitter_id>TCGA-XF-AAML</submitter_id>
			</item>
		</hits>
		<pagination>
			<count>5</count>
			<sort/>
			<from>0</from>
			<pages>2811</pages>
			<total>14052</total>
			<page>1</page>
			<size>5</size>
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
{"data": {"hits": [{"id": "f7af65fc-97e3-52ce-aa2c-b707650e747b", "submitter_id": "TARGET-00-NAAEMA"}, {"id": "513d0a2a-3c94-5a36-97a4-24c3656fc66e", "submitter_id": "TARGET-00-NAAEMB"}, {"id": "b5f20676-727b-50b0-9b5a-582cd8572d6d", "submitter_id": "TARGET-00-NAAEMC"}, {"id": "0c0b183f-0d4a-5a9d-9888-0617cebcc462", "submitter_id": "TARGET-20-PABGKN"}, {"id": "0f5ed7a7-226d-57bc-a4ce-8a6b18560c55", "submitter_id": "TARGET-20-PABHET"}], "pagination": {"count": 5, "sort": "submitter_id:asc", "from": 0, "page": 1, "total": 14551, "pages": 2911, "size": 5}}, "warnings": {}}
```
```Request2
curl  'https://api.gdc.cancer.gov/cases?fields=submitter_id&sort=submitter_id:asc&size=5&pretty=true'
```
```Response2
{
  "data": {
    "hits": [
      {
        "id": "f7af65fc-97e3-52ce-aa2c-b707650e747b",
        "submitter_id": "TARGET-00-NAAEMA"
      },
      {
        "id": "513d0a2a-3c94-5a36-97a4-24c3656fc66e",
        "submitter_id": "TARGET-00-NAAEMB"
      },
      {
        "id": "b5f20676-727b-50b0-9b5a-582cd8572d6d",
        "submitter_id": "TARGET-00-NAAEMC"
      },
      {
        "id": "0c0b183f-0d4a-5a9d-9888-0617cebcc462",
        "submitter_id": "TARGET-20-PABGKN"
      },
      {
        "id": "0f5ed7a7-226d-57bc-a4ce-8a6b18560c55",
        "submitter_id": "TARGET-20-PABHET"
      }
    ],
    "pagination": {
      "count": 5,
      "sort": "submitter_id:asc",
      "from": 0,
      "page": 1,
      "total": 14551,
      "pages": 2911,
      "size": 5
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
        "file_name": "NARKY_p_TCGAb69_SNP_N_GenomeWideSNP_6_H03_697832.grch38.seg.txt",
        "cases": [
          {
            "submitter_id": "TCGA-BP-4989"
          }
        ],
        "file_id": "3bd4d5dc-563a-481c-87a6-ec0017d0d58a",
        "file_size": 54200
      },
      {
        "file_name": "652ecf99-1af9-41fc-b0a5-d3e5c07a7b5d.FPKM.txt.gz",
        "cases": [
          {
            "submitter_id": "TCGA-60-2709"
          }
        ],
        "file_id": "b3286166-01f9-4149-81b5-a2ea5f27c50e",
        "file_size": 530665
      },
      {
        "file_name": "CUSKS_p_TCGAb47_SNP_1N_GenomeWideSNP_6_D05_628212.nocnv_grch38.seg.txt",
        "cases": [
          {
            "submitter_id": "TCGA-A8-A07Z"
          }
        ],
        "file_id": "282cc9d1-c5e9-49ff-b27b-e00c1e5529c6",
        "file_size": 15806
      },
      {
        "file_name": "REEDY_p_TCGAb65_SNP_N_GenomeWideSNP_6_F01_697686.nocnv_grch38.seg.txt",
        "cases": [
          {
            "submitter_id": "TCGA-CJ-4871"
          }
        ],
        "file_id": "fe44a644-eefc-42c5-aac7-a216bc1e88e1",
        "file_size": 6179
      },
      {
        "file_name": "84df7a8fee9fedb5e8e22849ec66d294_gdc_realn.bam",
        "cases": [
          {
            "submitter_id": "TCGA-A2-A0CO"
          }
        ],
        "file_id": "acd0ec73-c1fe-463e-912c-84e8416510e5",
        "file_size": 15545555724
      },
      {
        "file_name": "ed8c4bb6-891a-4cf2-80ba-42c5594760d0.vcf",
        "cases": [
          {
            "submitter_id": "TCGA-BQ-7059"
          }
        ],
        "file_id": "ed8c4bb6-891a-4cf2-80ba-42c5594760d0",
        "file_size": 264694
      },
      {
        "file_name": "nationwidechildrens.org_clinical.TCGA-IG-A6QS.xml",
        "cases": [
          {
            "submitter_id": "TCGA-IG-A6QS"
          }
        ],
        "file_id": "fe8cf009-f033-4536-95c7-836adcba5bf3",
        "file_size": 36996
      },
      {
        "file_name": "05f6f9f7-6fb7-4c95-b79c-fdfaba16539d.vep.reheader.vcf.gz",
        "cases": [
          {
            "submitter_id": "TCGA-DK-A3IV"
          }
        ],
        "file_id": "05f6f9f7-6fb7-4c95-b79c-fdfaba16539d",
        "file_size": 415044
      },
      {
        "file_name": "C484.TCGA-12-5301-01A-01D-1486-08.7_gdc_realn.bam",
        "cases": [
          {
            "submitter_id": "TCGA-12-5301"
          }
        ],
        "file_id": "3b0293c2-4a26-428c-b097-9489f23a2a2d",
        "file_size": 23661175335
      },
      {
        "file_name": "75a36e71-400d-46a5-93b0-7813cf0595ea.FPKM.txt.gz",
        "cases": [
          {
            "submitter_id": "TCGA-BF-A5EO"
          }
        ],
        "file_id": "28f763c7-8064-4151-ae0e-31e70cd9bfe8",
        "file_size": 488422
      }
    ],
    "pagination": {
      "count": 10,
      "sort": "",
      "from": 0,
      "page": 1,
      "total": 216435,
      "pages": 21644,
      "size": 10
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
```
{
  "data": {
    "data_type": "Aligned Reads",
    "updated_datetime": "2016-09-18T04:25:13.163601-05:00",
    "created_datetime": "2016-05-26T18:55:53.506549-05:00",
    "file_name": "000aa811c15656604161e8f0e3a0aae4_gdc_realn.bam",
    "md5sum": "200475f5f6e42520204e5f6aadfe954f",
    "data_format": "BAM",
    "acl": [
      "phs000178"
    ],
    "access": "controlled",
    "platform": "Illumina",
    "state": "submitted",
    "file_id": "ac2ddebd-5e5e-4aea-a430-5a87c6d9c878",
    "data_category": "Raw Sequencing Data",
    "file_size": 12667634731,
    "cases": [
      {
        "samples": [
          {
            "sample_type_id": "11",
            "updated_datetime": "2016-09-08T11:00:45.021005-05:00",
            "time_between_excision_and_freezing": null,
            "oct_embedded": "false",
            "tumor_code_id": null,
            "submitter_id": "TCGA-QQ-A5VA-11A",
            "intermediate_dimension": null,
            "sample_id": "b4e7558d-898e-4d68-a897-381edde0bbcc",
            "is_ffpe": false,
            "pathology_report_uuid": null,
            "created_datetime": null,
            "tumor_descriptor": null,
            "sample_type": "Solid Tissue Normal",
            "state": null,
            "current_weight": null,
            "composition": null,
            "time_between_clamping_and_freezing": null,
            "shortest_dimension": null,
            "tumor_code": null,
            "tissue_type": null,
            "days_to_sample_procurement": null,
            "freezing_method": null,
            "preservation_method": null,
            "days_to_collection": 5980,
            "initial_weight": 810.0,
            "longest_dimension": null
          }
        ]
      }
    ],
    "submitter_id": "32872121-d38a-4128-b96a-698a6f18f29d",
    "type": "aligned_reads",
    "file_state": "processed",
    "experimental_strategy": "WXS"
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
        "file_name": "unc.edu.276a1e00-cf3a-4463-a97b-d544381219ea.2363081.rsem.isoforms.normalized_results"
      },
      {
        "file_name": "nationwidechildrens.org_clinical.TCGA-EY-A5W2.xml"
      }
    ],
    "pagination": {
      "count": 2,
      "sort": "",
      "from": 0,
      "pages": 300936,
      "total": 601872,
      "page": 1,
      "size": 2
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
        "file_name": "OCULI_p_TCGA_159_160_SNP_N_GenomeWideSNP_6_E09_831242.grch38.seg.txt",
        "id": "1d959137-d8e6-4336-b357-8ab9c88eeca8"
      },
      {
        "file_name": "jhu-usc.edu_SKCM.HumanMethylation450.3.lvl-3.TCGA-EE-A3JI-06A-11D-A21B-05.gdc_hg38.txt",
        "id": "9c02ec95-4aa3-4112-8823-c0fa87f71773"
      },
      {
        "file_name": "jhu-usc.edu_LAML.HumanMethylation450.2.lvl-3.TCGA-AB-3002-03A-01D-0742-05.gdc_hg38.txt",
        "id": "731c3560-bcef-4ebf-bfbc-7320399a5bcb"
      },
      {
        "file_name": "CUSKS_p_TCGAb47_SNP_1N_GenomeWideSNP_6_B03_628222.grch38.seg.txt",
        "id": "a6f73a3e-faf8-49d9-9b68-77781bd302df"
      },
      {
        "file_name": "5496e9f1-a383-4874-95bb-f4d1b33f4594.vcf",
        "id": "5496e9f1-a383-4874-95bb-f4d1b33f4594"
      }
    ],
    "pagination": {
      "count": 5,
      "sort": "",
      "from": 101,
      "page": 21,
      "total": 274724,
      "pages": 54945,
      "size": 5
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
        "id": "f7af65fc-97e3-52ce-aa2c-b707650e747b",
        "submitter_id": "TARGET-00-NAAEMA"
      },
      {
        "id": "513d0a2a-3c94-5a36-97a4-24c3656fc66e",
        "submitter_id": "TARGET-00-NAAEMB"
      },
      {
        "id": "b5f20676-727b-50b0-9b5a-582cd8572d6d",
        "submitter_id": "TARGET-00-NAAEMC"
      },
      {
        "id": "0c0b183f-0d4a-5a9d-9888-0617cebcc462",
        "submitter_id": "TARGET-20-PABGKN"
      },
      {
        "id": "0f5ed7a7-226d-57bc-a4ce-8a6b18560c55",
        "submitter_id": "TARGET-20-PABHET"
      },
      {
        "id": "b2a560a4-5e52-5d78-90ef-d680fbaf44d0",
        "submitter_id": "TARGET-20-PABHKY"
      },
      {
        "id": "1e5c8323-383d-51a0-9199-1b9504b29c7e",
        "submitter_id": "TARGET-20-PABLDZ"
      },
      {
        "id": "c550a267-30bd-5bf3-9699-61341559e0d5",
        "submitter_id": "TARGET-20-PACDZR"
      },
      {
        "id": "0fe29a81-74fc-5158-ae13-0437bc272805",
        "submitter_id": "TARGET-20-PACEGD"
      },
      {
        "id": "dd2b23ec-46f4-56b2-9429-6015c6dc730f",
        "submitter_id": "TARGET-20-PADDXZ"
      }
    ],
    "pagination": {
      "count": 10,
      "sort": "submitter_id:asc",
      "from": 0,
      "page": 1,
      "total": 14551,
      "pages": 1456,
      "size": 10
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
    "pagination": {
      "count": 0,
      "sort": "program.name:asc",
      "from": 0,
      "page": 1,
      "total": 39,
      "pages": 39,
      "size": 0
    },
    "hits": [],
    "aggregations": {
      "program.name": {
        "buckets": [
          {
            "key": "TCGA",
            "doc_count": 33
          },
          {
            "key": "TARGET",
            "doc_count": 6
          }
        ]
      }
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
{
  "data": {
    "query": {
      "hits": [
        {
          "disease_type": [
            "Esophageal Carcinoma"
          ],
          "id": "UHJvamVjdDpUQ0dBLUVTQ0E=",
          "name": "Esophageal Carcinoma",
          "primary_site": [
            "Esophagus"
          ],
          "project_id": "TCGA-ESCA"
        },
        {
          "disease_type": [
            "Head and Neck Squamous Cell Carcinoma"
          ],
          "id": "UHJvamVjdDpUQ0dBLUhOU0M=",
          "name": "Head and Neck Squamous Cell Carcinoma",
          "primary_site": [
            "Head and Neck"
          ],
          "project_id": "TCGA-HNSC"
        },
        {
          "disease_type": [
            "Liver Hepatocellular Carcinoma"
          ],
          "id": "UHJvamVjdDpUQ0dBLUxJSEM=",
          "name": "Liver Hepatocellular Carcinoma",
          "primary_site": [
            "Liver"
          ],
          "project_id": "TCGA-LIHC"
        },
        {
          "disease_type": [
            "Colon Adenocarcinoma"
          ],
          "id": "UHJvamVjdDpUQ0dBLUNPQUQ=",
          "name": "Colon Adenocarcinoma",
          "primary_site": [
            "Colorectal"
          ],
          "project_id": "TCGA-COAD"
        },
        {
          "disease_type": [
            "Adrenocortical Carcinoma"
          ],
          "id": "UHJvamVjdDpUQ0dBLUFDQw==",
          "name": "Adrenocortical Carcinoma",
          "primary_site": [
            "Adrenal Gland"
          ],
          "project_id": "TCGA-ACC"
        }
      ]
    }
  }
}
```

This endpoint can be used to quickly retrieve information about a file.  For example, if a user wanted to know the UUID for `nationwidechildrens.org_biospecimen.TCGA-EL-A4K1.xml`, the following query could be used to quickly retrieve it programmatically:

```Shell
curl "https://api.gdc.cancer.gov/v0/all?query=nationwidechildrens.org_biospecimen.TCGA-EL-A4K1.xml&size=5"
```
```Response
{
  "data": {
    "query": {
      "hits": [
        {
          "file_id": "2a7a354b-e497-4ae6-8a85-a170951596c1",
          "file_name": "nationwidechildrens.org_biospecimen.TCGA-EL-A4K1.xml",
          "id": "RmlsZToyYTdhMzU0Yi1lNDk3LTRhZTYtOGE4NS1hMTcwOTUxNTk2YzE=",
          "submitter_id": null
        }
      ]
    }
  }
}
```

## Additional Examples

More examples of API functionality described in this section are provided in [Additional Examples](Additional_Examples.md).
