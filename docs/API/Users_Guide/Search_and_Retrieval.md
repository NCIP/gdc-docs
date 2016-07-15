# Search and Retrieval

## About Search and Retrieval API Queries

API queries that search and retrieve information stored in the GDC are constructed using [endpoints](#query-endpoints), [parameters](#query-parameters), and [filtering operators](#filtering-operators).

Queries can be executed using HTTP GET or HTTP POST. For simplicity, most examples in this section use GET, but for queries larger than a certain size POST is the only method that works, due to URL length limits. [This example](#example-http-post-request) explains how to construct POST API queries to search and retrieve information from the GDC.

**Note:** Requests to search and retrieve information stored in the GDC Legacy Archive must be directed to `legacy/` endpoints. See [Getting Started](Getting_Started.md#gdc-legacy-archive) for details.

**Note:** Queries described in this section work only on datasets that have been released to the GDC Data Portal. Data that is in the process of being submitted to GDC and is only available on the GDC Submission Portal cannot be queried using these methods. See [Submission](Submission.md) for information on how data submitters can query their unreleased data using GraphQL.

## Query Endpoints
The following search and retrieval endpoints are available on the GDC API:

| Endpoints | Description |
| --- | --- |
| files | Information about files stored in the GDC |
| cases | Information related to cases, or sample donors. |
| projects | Information about projects |
| annotations | Information about annotations to GDC data |
| \_mapping | Information about elements that can be used to query other endpoints |

### Project Endpoint
The GDC Project Endpoint `https://gdc-api.nci.nih.gov/projects` provides overall access to all the data served by GDC organized by Project such project(study) name, program, disease, primary site and state.

#### Example
This example is a query for projects contained in the GDC. It uses the [from](#from), [size](#size), [sort](#sort), and [pretty](#pretty) parameters, and returns the first two projects sorted by project id.

```shell
curl 'https://gdc-api.nci.nih.gov/projects?from=1&size=2&sort=project.project_id:asc&pretty=true'
```
``` Output

	{
	  "data": {
	    "hits": [
	      {
	        "state": "legacy",
	        "project_id": "TCGA-ACC",
	        "primary_site": "Adrenal Gland",
	        "disease_type": "Adrenocortical Carcinoma",
	        "name": "Adrenocortical Carcinoma"
	      },
	      {
	        "dbgap_accession_number": "phs000464",
	        "disease_type": "Acute Lymphoblastic Leukemia",
	        "state": "legacy",
	        "primary_site": "Blood",
	        "project_id": "TARGET-ALL-P2",
	        "name": "Acute Lymphoblastic Leukemia - Phase II"
	      }
	    ],
	    "pagination": {
	      "count": 2,
	      "sort": "project.project_id:asc",
	      "from": 1,
	      "pages": 22,
	      "total": 44,
	      "page": 1,
	      "size": 2
	    }
	  },
	  "warnings": {}
	}
```

#### Retrieval of project metadata using project_id

The `project` endpoint supports a simple query format that retrieves the metadata of a single project using its `project_id`:

```shell
curl 'https://gdc-api.nci.nih.gov/projects/TARGET-NBL?expand=summary,summary.experimental_strategies,summary.data_categories&pretty=true'
```
```Response
{
  "data": {
    "dbgap_accession_number": "phs000467",
    "name": "Neuroblastoma",
    "summary": {
      "data_categories": [
        {
          "case_count": 151,
          "file_count": 471,
          "data_category": "Transcriptome Profiling"
        },
        {
          "case_count": 216,
          "file_count": 1728,
          "data_category": "Simple Nucleotide Variation"
        },
        {
          "case_count": 1120,
          "file_count": 1,
          "data_category": "Clinical"
        },
        {
          "case_count": 270,
          "file_count": 599,
          "data_category": "Raw Sequencing Data"
        }
      ],
      "case_count": 1120,
      "file_count": 2799,
      "experimental_strategies": [
        {
          "case_count": 221,
          "file_count": 2170,
          "experimental_strategy": "WXS"
        },
        {
          "case_count": 151,
          "file_count": 628,
          "experimental_strategy": "RNA-Seq"
        }
      ],
      "file_size": 8157089415961
    },
    "released": true,
    "state": "legacy",
    "primary_site": "Nervous System",
    "project_id": "TARGET-NBL",
    "disease_type": "Neuroblastoma"
  },
  "warnings": {}
}
```

### Files Endpoint

The GDC Files Endpoint `https://gdc-api.nci.nih.gov/files` enables search and retrieval of information relating to files stored in the GDC, including file properties such as `file_name`, `md5sum`, `data_format`, and others.

#### Example

This example is a query for files contained in the GDC. It uses the [from](#from), [size](#size), [sort](#sort), and [pretty](#pretty) parameters, and returns only the first two files, sorted by file size, from smallest to largest.

```shell
curl 'https://gdc-api.nci.nih.gov/files?from=1&size=2&sort=file_size:asc&pretty=true'
```
``` Output
	{
	  "data": {
	    "hits": [
	      {
	        "origin": "migrated",
	        "data_type": "Raw microarray data",
	        "platform": "HG-U133_Plus_2",
	        "file_name": "TCGA-AB-2842-03A-01R-0757-21.CEL.README",
	        "md5sum": "56f9a6d58b450bf7e9f6431a86220b9d",
	        "data_format": "CEL",
	        "acl": "open",
	        "access": "open",
	        "uploaded_datetime": 1425340539,
	        "state": "live",
	        "data_subtype": "Raw intensities",
	        "file_id": "ca13321c-02aa-4141-bdb6-84d31e3c5711",
	        "file_size": 43,
	        "experimental_strategy": "Gene expression array"
	      },
	      {
	        "origin": "migrated",
	        "data_type": "Raw microarray data",
	        "platform": "HG-U133_Plus_2",
	        "file_name": "TCGA-AB-2809-03A-01R-0757-21.CEL.README",
	        "md5sum": "56f9a6d58b450bf7e9f6431a86220b9d",
	        "data_format": "CEL",
	        "acl": "open",
	        "access": "open",
	        "uploaded_datetime": 1425340539,
	        "state": "live",
	        "data_subtype": "Raw intensities",
	        "file_id": "299d500b-49e2-4c62-9111-c0691592dce1",
	        "file_size": 43,
	        "experimental_strategy": "Gene expression array"
	      }
	    ],
	    "pagination": {
	      "count": 2,
	      "sort": "file_size:asc",
	      "from": 1,
	      "pages": 300936,
	      "total": 601872,
	      "page": 1,
	      "size": 2
	    }
	  },
	  "warnings": {}
	}
```

#### Retrieval of file metadata using individual UUIDs:

The `files` endpoint supports a simple query format that retrieves the metadata of a single file using its UUID:

```Shell
curl 'https://gdc-api.nci.nih.gov/files/000225ad-497b-4a8c-967e-a72159c9b3c9?pretty=true'
```
```
{
  "data": {
    "data_type": "Raw Simple Somatic Mutation",
    "updated_datetime": "2016-06-04T23:42:25.428738-05:00",
    "created_datetime": "2016-06-03T19:04:32.950673-05:00",
    "file_name": "000225ad-497b-4a8c-967e-a72159c9b3c9.snp.Somatic.hc.vcf.gz",
    "md5sum": "bbe8a7157acbfc9133e47898650b5437",
    "data_format": "VCF",
    "acl": [
      "phs000178"
    ],
    "access": "controlled",
    "state": "submitted",
    "file_id": "000225ad-497b-4a8c-967e-a72159c9b3c9",
    "data_category": "Simple Nucleotide Variation",
    "file_size": 19690,
    "submitter_id": "TCGA-VR-A8ET-01A-11D-A403-09_TCGA-VR-A8ET-10B-01D-A403-09_varscan",
    "type": "simple_somatic_mutation",
    "file_state": "processed",
    "experimental_strategy": "WXS"
  },
  "warnings": {}
}
```

#### files/ids Endpoint

The `files/ids` endpoint corresponds to the "Quick Search" functionality of the GDC Data Portal. The API response includes all files for which the query matches the beginning (or entirety) of any of the following fields:

	project.project_id
	project.name
	project.disease_type.analyzed
	project.primary_site.analyzed
	case.aliquot_ids
	case.submitter_aliquot_ids
	case.analyte_ids
	case.submitter_analyte_ids
	case.case_id.raw
	case.submitter_id.raw
	case.portion_ids
	case.submitter_portion_ids
	case.sample_ids
	case.slide_ids
	case.submitter_slide_ids
	case.submitter_sample_ids
	file.file_id.raw
	file.file_name.raw
	file.submitter_id

Requests to this endpoint must be in the format `files/ids?query=`, for example:

```shell
curl 'https://gdc-api.nci.nih.gov/files/ids?query=nationwidechildrens.org_clinical.TCGA-EM&pretty=true'
```
``` Output
{
  "data": {
    "pagination": {
      "count": 5,
      "sort": "",
      "from": 1,
      "page": 1,
      "total": 81,
      "pages": 17,
      "size": 5
    },
    "hits": [
      {
        "_type": "file",
        "file_name": "nationwidechildrens.org_clinical.TCGA-EM-A3FQ.xml",
        "file_id": "efac6904-ac9f-4a44-bf9c-f7d9a822c127",
        "_score": 4.644438,
        "cases": [
          {
            "case_id": "fef9c64f-5959-4da0-aaa2-66b56fc7b4c3",
            "submitter_id": "TCGA-EM-A3FQ"
          }
        ],
        "_id": "efac6904-ac9f-4a44-bf9c-f7d9a822c127"
      },
      {
        "_type": "file",
        "file_name": "nationwidechildrens.org_clinical.TCGA-EM-A4FN.xml",
        "file_id": "07add35d-66f0-4384-bb2c-9d86661f4073",
        "_score": 4.644438,
        "cases": [
          {
            "case_id": "f854ce67-c586-4424-a674-2dd67ad0ed7f",
            "submitter_id": "TCGA-EM-A4FN"
          }
        ],
        "_id": "07add35d-66f0-4384-bb2c-9d86661f4073"
      },
      {
        "_type": "file",
        "file_name": "nationwidechildrens.org_clinical.TCGA-EM-A3AN.xml",
        "file_id": "889f222e-09d1-477c-a7b9-a514b65f322b",
        "_score": 4.644438,
        "cases": [
          {
            "case_id": "6491d025-c061-4180-bfd4-d7c6e6e55f66",
            "submitter_id": "TCGA-EM-A3AN"
          }
        ],
        "_id": "889f222e-09d1-477c-a7b9-a514b65f322b"
      },
      {
        "_type": "file",
        "file_name": "nationwidechildrens.org_clinical.TCGA-EM-A4FO.xml",
        "file_id": "37db38d4-64ca-4cd4-9adc-d504b812997b",
        "_score": 4.644438,
        "cases": [
          {
            "case_id": "e0ab95ef-5c96-4d00-8950-04e26e3b4672",
            "submitter_id": "TCGA-EM-A4FO"
          }
        ],
        "_id": "37db38d4-64ca-4cd4-9adc-d504b812997b"
      },
      {
        "_type": "file",
        "file_name": "nationwidechildrens.org_clinical.TCGA-EM-A3FN.xml",
        "file_id": "e2d13acf-3121-4c7b-b2ec-28db49eff699",
        "_score": 4.644438,
        "cases": [
          {
            "case_id": "59cd0969-a798-4d19-95ed-a311f39d2f38",
            "submitter_id": "TCGA-EM-A3FN"
          }
        ],
        "_id": "e2d13acf-3121-4c7b-b2ec-28db49eff699"
      }
    ],
    "_shards": {
      "successful": 5,
      "failed": 0,
      "total": 5
    },
    "took": 9,
    "timed_out": false
  },
  "warnings": {}
}
```

### Cases Endpoint

The GDC Cases Endpoint `https://gdc-api.nci.nih.gov/cases` enables search and retrieval of information related to a specific case.




#### Example

This example is a query for files contained in GDC. It returns case where submitter id is `TCGA-BH-A0EA`, using the [pretty](#pretty) and [filters](#filters) parameters and the following [filtering operators](#filtering-operators):

	{"op":"and","content":[{"op":"in","content":{"field":"submitter_id","value":["TCGA-BH-A0EA"]}}]}

Command:

```shell
curl 'https://gdc-api.nci.nih.gov/cases?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22submitter_id%22%2C%22value%22%3A%5B%22TCGA-BH-A0EA%22%5D%7D%7D%5D%7D%0A%0A&pretty=true'
```
``` Output
{
"data": {
	"hits": [
		{
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
		}
	],
	"pagination": {
		"count": 1,
		"sort": "",
		"from": 1,
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
curl 'https://gdc-api.nci.nih.gov/cases/1f601832-eee3-48fb-acf5-80c4a454f26e?pretty=true&expand=diagnoses'
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

### Annotations Endpoint

The GDC Annotation Endpoint `https://gdc-api.nci.nih.gov/annotations` enables search and retrieval of annotations stored in the GDC.


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
curl 'https://gdc-api.nci.nih.gov/annotations?filters=%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22entity_id%22%2C%22value%22%3A%5B%22e0d36cc0-652c-4224-bb10-09d15c7bd8f1%22%2C%2225ebc29a-7598-4ae4-ba7f-618d448882cc%22%2C%22fe660d7c-2746-4b50-ab93-b2ed99960553%22%5D%7D%7D&pretty=true'
```
``` Output
{
  "data": {
    "hits": [
      {
        "status": "Approved",
        "category": "Item flagged DNU",
        "entity_id": "fe660d7c-2746-4b50-ab93-b2ed99960553",
        "classification": "CenterNotification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2015-09-28T13:39:13-04:00",
        "annotation_id": "5ddadefe-8b57-5ce2-b8b2-918d63d99a59",
        "notes": "The aliquot failed Broad pipeline QC and not all files are suitable for use. Consult the SDRF file to determine which files are usable.",
        "entity_type": "aliquot",
        "submitter_id": "29087",
        "case_id": "41b59716-116f-4942-8b63-409870a87e26",
        "entity_submitter_id": "TCGA-DK-A3IM-10A-01D-A20B-01",
        "case_submitter_id": "TCGA-DK-A3IM"
      },
      {
        "status": "Approved",
        "category": "Item is noncanonical",
        "entity_id": "25ebc29a-7598-4ae4-ba7f-618d448882cc",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2012-07-12T15:00:15-04:00",
        "annotation_id": "d6500f94-618f-5334-a810-ade76b887ec9",
        "notes": "No Matching Normal",
        "entity_type": "sample",
        "submitter_id": "8009",
        "case_id": "bd114e05-5a97-41e2-a0d5-5d39a1e9d461",
        "entity_submitter_id": "TCGA-08-0514-01A",
        "case_submitter_id": "TCGA-08-0514"
      },
      {
        "status": "Approved",
        "category": "Prior malignancy",
        "entity_id": "e0d36cc0-652c-4224-bb10-09d15c7bd8f1",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2013-03-12T10:05:14-04:00",
        "annotation_id": "33336cdf-2cf0-5af2-bb52-fecd3427f180",
        "notes": "Patient had a prior lymphoma. Unknown radiation or systemic chemotherapy.",
        "entity_type": "case",
        "submitter_id": "15630",
        "case_id": "e0d36cc0-652c-4224-bb10-09d15c7bd8f1",
        "entity_submitter_id": "TCGA-FS-A1ZF",
        "case_submitter_id": "TCGA-FS-A1ZF"
      }
    ],
    "pagination": {
      "count": 3,
      "sort": "",
      "from": 1,
      "page": 1,
      "total": 3,
      "pages": 1,
      "size": 10
    }
  },
  "warnings": {}
}
```

#### Example

This example is a query for any annotations that are associated with the following *cases*, **directly or via a child entity**:

* the case with UUID 513c5f34-dc6e-4caa-81cc-907fd6a825b1
* the case with UUID 942c0088-c9a0-428c-a879-e16f8c5bfdb8

The query uses the [filters](#filters) parameter to specify entity UUIDs. Code samples below include the bare and percent-encoded filter JSON.

```Filter-JSON
{
   "op":"in",
   "content":{
      "field":"annotation.case_id",
      "value":[
         "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
         "942c0088-c9a0-428c-a879-e16f8c5bfdb8"
      ]
   }
}
```
```Filter-JSON-percent-encoded
%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22annotation.case_id%22%2C%22value%22%3A%5B%22513c5f34-dc6e-4caa-81cc-907fd6a825b1%22%2C%22942c0088-c9a0-428c-a879-e16f8c5bfdb8%22%5D%7D%7D
```
```shell
curl 'https://gdc-api.nci.nih.gov/annotations?filters=%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22annotation.case_id%22%2C%22value%22%3A%5B%22513c5f34-dc6e-4caa-81cc-907fd6a825b1%22%2C%22942c0088-c9a0-428c-a879-e16f8c5bfdb8%22%5D%7D%7D&pretty=true&size=30'
```
``` Output
{
  "data": {
    "hits": [
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "39c2e6c5-379e-4ffa-b02c-e1db298b34f7",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:25-04:00",
        "annotation_id": "8bed6748-11dd-5767-9a9c-577712f9b616",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "sample",
        "submitter_id": "22104",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-10A",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "e2bca154-9ff3-42a6-a5bc-0daa14080c92",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:32-04:00",
        "annotation_id": "a6997495-8544-58b1-81aa-7e628e46af7e",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "aliquot",
        "submitter_id": "22118",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-01A-11D-2395-08",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "History of unacceptable prior treatment related to a prior/other malignancy",
        "entity_id": "942c0088-c9a0-428c-a879-e16f8c5bfdb8",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2012-11-10T05:30:37-05:00",
        "annotation_id": "6d57211a-402b-52c9-b857-665164c63339",
        "notes": "Pt with synchronous B-cell lymphoma treated with cyclophosphamide, Adriamycin, vincristine, prednisone, and Rituxin prior to procurment of TCGA tumor.",
        "entity_type": "case",
        "submitter_id": "12063",
        "case_id": "942c0088-c9a0-428c-a879-e16f8c5bfdb8",
        "entity_submitter_id": "TCGA-CJ-4642",
        "case_submitter_id": "TCGA-CJ-4642"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "265c0257-54b6-4b79-8c9a-d56ca8bbdb48",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:41-04:00",
        "annotation_id": "0fe698fa-169a-51e2-b41b-6452d0b5cefc",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "aliquot",
        "submitter_id": "22122",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-01A-11R-A28V-07",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "a4d9f761-ece9-4a6b-8818-ecf4a8f7d380",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:28-04:00",
        "annotation_id": "6dd575b6-78a7-55ca-9dbe-654079eeca0f",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "analyte",
        "submitter_id": "22110",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-10A-01D",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "a5370483-6bdd-4e2a-9b4b-ec1eb381052b",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:25-04:00",
        "annotation_id": "3d268354-9193-50a9-a8b5-79f00d36d3c6",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "sample",
        "submitter_id": "22103",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-01A",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "2aa6f51a-e0fa-4bdb-830c-80b79c59cfc9",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:32-04:00",
        "annotation_id": "d4f40e94-f701-58d1-aa0c-b3d687844c7b",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "aliquot",
        "submitter_id": "22117",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-01A-11D-2391-01",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "d9fffb10-ae12-43ad-983b-1e2336f915a3",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:43-04:00",
        "annotation_id": "ab91b0d6-28fd-5600-8933-3457c0939687",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "aliquot",
        "submitter_id": "22126",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-10A-01D-2395-08",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "7402c87f-c9b7-451f-9b8d-a8979bf0d98b",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:30-04:00",
        "annotation_id": "b16b9f8f-adfd-5ab8-b10e-16b0ed277dbc",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "analyte",
        "submitter_id": "22114",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-11A-01R",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "Synchronous malignancy",
        "entity_id": "942c0088-c9a0-428c-a879-e16f8c5bfdb8",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2012-11-10T05:30:35-05:00",
        "annotation_id": "b02b2a31-7c5a-5ff0-a2ea-6ab54b721a86",
        "notes": "Pt with synchronous B-cell lymphoma treated with cyclophosphamide, Adriamycin, vincristine, prednisone, and Rituxin prior to procurment of TCGA tumor.",
        "entity_type": "case",
        "submitter_id": "12062",
        "case_id": "942c0088-c9a0-428c-a879-e16f8c5bfdb8",
        "entity_submitter_id": "TCGA-CJ-4642",
        "case_submitter_id": "TCGA-CJ-4642"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "5a676d7d-0b8d-48fb-b898-0eb023a871e1",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:27-04:00",
        "annotation_id": "96d82ca1-5e27-5048-82d5-ee5886437580",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "analyte",
        "submitter_id": "22107",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-01A-11H",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "2327bab1-e352-4c43-ab40-6162e8632c26",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:26-04:00",
        "annotation_id": "95873790-c215-5185-babe-0228e32d5689",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "analyte",
        "submitter_id": "22106",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-01A-11D",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "258a054e-b0ec-4215-81c9-8e1c715abc9b",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:30-04:00",
        "annotation_id": "f0a7ee98-09aa-5c46-b88a-66a09176c408",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "analyte",
        "submitter_id": "22113",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-11A-01H",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "7c9e5386-1540-40f2-96ed-c06d752598df",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:45-04:00",
        "annotation_id": "f34eb057-6f81-5d37-858c-9ed5cc8a2052",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "aliquot",
        "submitter_id": "22129",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-11A-01D-2391-01",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "77eccef5-254d-4ddc-a70f-576436278e0a",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:46-04:00",
        "annotation_id": "9bf047b4-bcfc-5bfc-bc42-584a1efd09b4",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "aliquot",
        "submitter_id": "22131",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-11A-01H-2402-13",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T08:54:02-04:00",
        "annotation_id": "53ccabc5-0c89-5e55-a267-130dd021d9ff",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "case",
        "submitter_id": "22102",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "9fe04d07-57ba-434f-8254-184fe00ab23e",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:46-04:00",
        "annotation_id": "100caf14-4fd6-51ef-8f3c-8e7ca62d0233",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "aliquot",
        "submitter_id": "22132",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-11A-01R-A28V-07",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "26368154-1ed7-4471-87b9-d50d5fe8046b",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:26-04:00",
        "annotation_id": "75e180f1-4216-591f-8bd4-426e923d0778",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "sample",
        "submitter_id": "22105",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-11A",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "Prior malignancy",
        "entity_id": "942c0088-c9a0-428c-a879-e16f8c5bfdb8",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2011-02-05T13:15:23-05:00",
        "annotation_id": "59f7f353-2119-58f2-a340-34bc39c3d7ef",
        "notes": "[intgen.org]: Prior Malignancy",
        "entity_type": "case",
        "submitter_id": "1272",
        "case_id": "942c0088-c9a0-428c-a879-e16f8c5bfdb8",
        "entity_submitter_id": "TCGA-CJ-4642",
        "case_submitter_id": "TCGA-CJ-4642"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "51ded92f-852f-4445-96b2-dd9db1544a9d",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:27-04:00",
        "annotation_id": "8ac6de9f-6b86-5421-a4be-e1fe8a4a8a38",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "analyte",
        "submitter_id": "22108",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-01A-11R",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "Molecular analysis outside specification",
        "entity_id": "942c0088-c9a0-428c-a879-e16f8c5bfdb8",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2010-10-29T00:00:00-04:00",
        "annotation_id": "e0fcd5e3-02f0-52b4-b949-95ed4de324f6",
        "notes": "Molecular results off spec",
        "entity_type": "case",
        "submitter_id": "831",
        "case_id": "942c0088-c9a0-428c-a879-e16f8c5bfdb8",
        "entity_submitter_id": "TCGA-CJ-4642",
        "case_submitter_id": "TCGA-CJ-4642"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "3cffe6b0-3aef-413a-a9a9-4a2fcefff6ad",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:45-04:00",
        "annotation_id": "dfad3f3e-9d62-51e0-a433-77f89b1c24ca",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "aliquot",
        "submitter_id": "22130",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-11A-01D-2395-08",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "6f688bcb-0534-4a50-a39d-a3a149165694",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:33-04:00",
        "annotation_id": "03716d9e-c294-570b-89ab-53af108850f3",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "aliquot",
        "submitter_id": "22120",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-01A-11H-2402-13",
        "case_submitter_id": "TCGA-56-8623"
      },
      {
        "status": "Approved",
        "category": "BCR Notification",
        "entity_id": "f1090757-a85d-4702-a16f-0ee35284c45d",
        "classification": "Notification",
        "updated_datetime": "2016-05-01T15:00:21.638875-05:00",
        "created_datetime": "2014-09-05T09:02:29-04:00",
        "annotation_id": "84236594-5c53-504d-bd27-1e8415cb768a",
        "notes": "Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues. The sequencing and characterization centers have reported that molecular signature of the tumor DNA and RNA cluster with the normal controls. Their analysis suggests that the purity of the tumor is very low. Therefore, the data from this patient should be used with caution.",
        "entity_type": "analyte",
        "submitter_id": "22112",
        "case_id": "513c5f34-dc6e-4caa-81cc-907fd6a825b1",
        "entity_submitter_id": "TCGA-56-8623-11A-01D",
        "case_submitter_id": "TCGA-56-8623"
      }
    ],
    "pagination": {
      "count": 24,
      "sort": "",
      "from": 1,
      "page": 1,
      "total": 24,
      "pages": 1,
      "size": 30
    }
  },
  "warnings": {}
}
```



### \_mapping Endpoint

Each Search and Retrieval endpoint is equipped with a ```_mapping``` endpoint that provides information about elements that can be used to query the Search and Retrieval endpoint.

The API response to a `_mapping` query is a list of objects, containing fields available through the API. The high-level structure of the response is as follows:

	"\_mapping": {}
	, "defaults": []
	, "expand": []
	, "fields": []
	, "multi": []
	, "nested": []

#### Example

```shell
curl 'https://gdc-api.nci.nih.gov/projects/_mapping'
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


## Query Parameters

The following query parameters can be used with all methods and resources in the GDC API. The use of any particular parameter is optional except where noted.

Parameter | Default | Description
--------- | ------- | -----------
format | JSON | Specifies the API response format: JSON, XML, or TSV
pretty | false | Returns response with indentations and line breaks in a human-readable format
fields | null | Query option to specify which fields to include in the response
size | 10 | Specifies the number of results to return
from   | 1 | Specifies the first record to return from the set resulting of a query
sort | null | Specifies sorting algorithm for the results in the API response
filters| null | Query option filters specify criteria for the returned response
facets | null | Provides a list of number of files available given current filters facet


### Format

Specifies the format of the API response. JSON is default, and `TSV` and `XML` options are available.

#### Examples

```shell1
curl  'https://gdc-api.nci.nih.gov/cases?fields=submitter_id&size=5&format=TSV'
```
```python1
import requests

cases_endpt = 'https://gdc-api.nci.nih.gov/cases'
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
curl  'https://gdc-api.nci.nih.gov/cases?fields=submitter_id&size=5&format=XML&pretty=true'
```
```python2
import requests

cases_endpt = 'https://gdc-api.nci.nih.gov/cases'
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
			<from>1</from>
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
curl  'https://gdc-api.nci.nih.gov/cases?fields=submitter_id&sort=submitter_id:asc&size=5'
```
```Response1
{"data": {"hits": [{"submitter_id": "TARGET-20-PABGKN"}, {"submitter_id": "TARGET-20-PABHET"}, {"submitter_id": "TARGET-20-PABHKY"}, {"submitter_id": "TARGET-20-PABLDZ"}, {"submitter_id": "TARGET-20-PACDZR"}], "pagination": {"count": 5, "sort": "submitter_id.raw:asc", "from": 1, "pages": 2811, "total": 14052, "page": 1, "size": 5}}, "warnings": {}}
```
```Request2
curl  'https://gdc-api.nci.nih.gov/cases?fields=submitter_id&sort=submitter_id:asc&size=5&pretty=true'
```
```Response2
{
  "data": {
    "hits": [
      {
        "submitter_id": "TARGET-20-PABGKN"
      },
      {
        "submitter_id": "TARGET-20-PABHET"
      },
      {
        "submitter_id": "TARGET-20-PABHKY"
      },
      {
        "submitter_id": "TARGET-20-PABLDZ"
      },
      {
        "submitter_id": "TARGET-20-PACDZR"
      }
    ],
    "pagination": {
      "count": 5,
      "sort": "submitter_id.raw:asc",
      "from": 1,
      "pages": 2811,
      "total": 14052,
      "page": 1,
      "size": 5
    }
  },
  "warnings": {}
}
```

### Fields

This query parameter specifies which fields are to be included in the API response. A listing of available fields for each endpoint is provided in [Appendix A](Appendix_A_Available_Fields.md).

#### Example

The following example requests case submitter ID, file UUID, file name and file size from the `files` endpoint.

```shell
curl 'https://gdc-api.nci.nih.gov/files?fields=cases.submitter_id,file_id,file_name,file_size&pretty=true'
```
```python
import requests
import json

files_endpt = 'https://gdc-api.nci.nih.gov/files'
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
      "from": 1,
      "page": 1,
      "total": 216435,
      "pages": 21644,
      "size": 10
    }
  },
  "warnings": {}
}
```

### Size and From

GDC API provides a pagination feature that limits the number of results returned by the API. It is implemented using `size` and `from` query parameters.

The `size` query parameter specifies the maximum number of results to return. Default `size` is 10. If the number of query results is greater than `size`, only some of the results will be returned.

The `from` query parameter specifies the first record to return out of the set of results. For example, if there are 20 cases returned from the `cases` endpoint, then setting `from` to `11` will return results 11 to 20. The `from` parameter can be used in conjunction with the `size` parameter to return a specific subset of results.


#### Example


``` Shell1
curl 'https://gdc-api.nci.nih.gov/files?fields=file_name&from=0&size=2&pretty=true'
```
``` Python1
import requests
import json

files_endpt = 'https://gdc-api.nci.nih.gov/files'
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
      "from": 1,
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
curl 'https://gdc-api.nci.nih.gov/files?fields=file_name&from=101&size=5&pretty=true'
```
``` Python2
import requests
import json

files_endpt = 'https://gdc-api.nci.nih.gov/files'
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
        "file_name": "unc.edu.956590e9-4962-497b-a59f-81ee0a1c0caf.1379536.junction_quantification.txt"
      },
      {
        "file_name": "MATZO_p_TCGAb40_SNP_1N_GenomeWideSNP_6_G09_667760.ismpolish.data.txt"
      },
      {
        "file_name": "GIRTH_p_TCGA_b108_137_SNP_N_GenomeWideSNP_6_D06_787864.hg18.seg.txt"
      },
      {
        "file_name": "PLENA_p_TCGAb63and64_SNP_N_GenomeWideSNP_6_B12_697382.CEL"
      },
      {
        "file_name": "TCGA-HU-8604-01A-11R-2402-13.isoform.quantification.txt"
      }
    ],
    "pagination": {
      "count": 5,
      "sort": "",
      "from": 100,
      "pages": 109553,
      "total": 547761,
      "page": 21,
      "size": 5
    }
  },
  "warnings": {}
}
```

### Sort

The `sort` query parameter sorts the results by a specific field, and with the sort direction specified using the `:asc` (ascending) or `:desc` (descending) prefix, e.g. `sort=field:desc`. A list of all valid _field_ names that can be used as facets is available in [Appendix A](Appendix_A_Available_Fields.md).

#### Example

Sort cases by `submitter_id` in ascending order:

``` shell
curl  'https://gdc-api.nci.nih.gov/cases?fields=submitter_id&sort=submitter_id:asc&pretty=true'
```
``` python
import requests
import json

cases_endpt = 'https://gdc-api.nci.nih.gov/cases'
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
        "submitter_id": "TARGET-20-PABGKN"
      },
      {
        "submitter_id": "TARGET-20-PABHET"
      },
      {
        "submitter_id": "TARGET-20-PABHKY"
      },
      {
        "submitter_id": "TARGET-20-PABLDZ"
      },
      {
        "submitter_id": "TARGET-20-PACDZR"
      },
      {
        "submitter_id": "TARGET-20-PACEGD"
      },
      {
        "submitter_id": "TARGET-20-PADDXZ"
      },
      {
        "submitter_id": "TARGET-20-PADYIR"
      },
      {
        "submitter_id": "TARGET-20-PADZCG"
      },
      {
        "submitter_id": "TARGET-20-PADZKD"
      }
    ],
    "pagination": {
      "count": 10,
      "sort": "submitter_id.raw:asc",
      "from": 1,
      "pages": 1406,
      "total": 14052,
      "page": 1,
      "size": 10
    }
  },
  "warnings": {}
}
```


### Filters

The `filters` parameter enables passing of complex search queries to the GDC API. The parameter carries a percent-encoded JSON object that contains a query.


#### Filtering Operators

Operators allow users to define query conditions.

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


Users can get a list of available values for a specific field in the filter by making a call to the appropriate API endpoint using the `facets` parameter, e.g. `https://gdc-api.nci.nih.gov/v0/cases?facets=demographic.gender&size=0&pretty=true`

#### Nested Operations

Filters support complex nested operations as well as simple queries on a single field. There are different types of operations available for many uses. For more examples see [Additional Examples](Additional_Examples.md).

#### Example: HTTP GET Request

This example requests `male` cases using HTTP GET.

The JSON object to be passed to the GDC API looks like:

	{"op": "=",
		  "content": {
			  "field": "cases.clinical.gender",
			  "value": ["male"]
		  }
	}

URL-encoding the above JSON object using [Percent-(URL)-encoding tool](http://text-rescue.com/string-escape/percent-url-encoding-tool.html) results in the following string:

	%7b%22op%22%3a+%22%3d%22%2c%0d%0a++++++%22content%22%3a+%7b%0d%0a++++++++++%22field%22%3a+%22cases.clinical.gender%22%2c%0d%0a++++++++++%22value%22%3a+%5b%22male%22%5d%0d%0a++++++%7d%0d%0a%7d

The above string can now be passed to the GDC API using the `filters` parameter:

```shell
 curl  'https://gdc-api.nci.nih.gov/cases?filters=%7b%22op%22%3a+%22%3d%22%2c%0d%0a++++++%22content%22%3a+%7b%0d%0a++++++++++%22field%22%3a+%22cases.clinical.gender%22%2c%0d%0a++++++++++%22value%22%3a+%5b%22male%22%5d%0d%0a++++++%7d%0d%0a%7d&pretty=true'
```
```python
import requests
import json
cases_endpt = 'https://gdc-api.nci.nih.gov/cases'
filt = {"op":"=",
        "content":{
            "field": "cases.clinical.gender",
            "value": ["male"]
        }
}
params = {'filters':json.dumps(filt), 'sort':'clinical.gender:asc'}
# requests URL-encodes automatically
response = requests.get(cases_endpt, params = params)
print json.dumps(response.json(), indent=2)
```
``` Output
{
  "data": {
    "hits": [
      {
        "sample_ids": "4cea843f-80ae-4284-b70f-8989717eee7c",
        "portion_ids": "eda1ff5e-a9e4-472a-8231-1a39ac2d2157",
        "submitter_portion_ids": "TCGA-KL-8344-01A-11",
        "submitter_aliquot_ids": "TCGA-KL-8344-01A-11W-2331-10",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-KL-8344-01A-11R",
        "analyte_ids": "f18fe926-0fb3-465b-98a5-40734ad550f9",
        "submitter_id": "TCGA-KL-8344",
        "case_id": "a5af1391-89d3-4f41-930e-7e8272afec98",
        "slide_ids": "58b31314-ca6e-450f-9fb2-f43e9fb27fc8",
        "submitter_sample_ids": "TCGA-KL-8344-01A",
        "project_id": "TCGA-KICH",
        "aliquot_ids": "3820fe10-c8a9-4250-b88d-07a511cf2751",
        "submitter_slide_ids": "TCGA-KL-8344-01A-01-TS1"
      },
      {
        "sample_ids": "70e1edd0-29a3-493a-8e03-107a5b00a4d7",
        "portion_ids": "a0e13920-f019-44b1-83c6-e8a2e508b315",
        "submitter_portion_ids": "TCGA-FA-A4XK-10A-01",
        "submitter_aliquot_ids": "TCGA-FA-A4XK-10A-01D-A31H-26",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-FA-A4XK-10A-01D",
        "analyte_ids": "521d5cb7-a976-4805-9855-4fc63378ee49",
        "submitter_id": "TCGA-FA-A4XK",
        "case_id": "a5b188f0-a6d3-4d4a-b04f-36d47ec05338",
        "slide_ids": "2b0fc256-4d70-4eb6-97e0-152ab7f6039e",
        "submitter_sample_ids": "TCGA-FA-A4XK-10A",
        "project_id": "TCGA-DLBC",
        "aliquot_ids": "8e020cb7-dab9-42b8-8ffe-c01fe58e0151",
        "submitter_slide_ids": "TCGA-FA-A4XK-01A-01-TSA"
      },
      {
        "sample_ids": "815ffaae-31a8-4064-b51a-17123a21e9e6",
        "portion_ids": "8c066e1f-1ed8-4ccb-8348-2a991ebe896b",
        "submitter_portion_ids": "TCGA-EF-5830-01A-21-1934-20",
        "submitter_aliquot_ids": "TCGA-EF-5830-01A-01D-1657-10",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-EF-5830-01A-01R",
        "analyte_ids": "bc9064c9-5af2-44ca-ab25-c01ffe03142c",
        "submitter_id": "TCGA-EF-5830",
        "case_id": "0b6b1937-4024-4f2c-aeca-46387277755f",
        "slide_ids": "b237171d-d2c7-413a-8c56-6ec374de9149",
        "submitter_sample_ids": "TCGA-EF-5830-01A",
        "project_id": "TCGA-READ",
        "aliquot_ids": "61fe838a-d218-407a-922a-7cb9b4fc8aa8",
        "submitter_slide_ids": "TCGA-EF-5830-01A-01-BS1"
      },
      {
        "sample_ids": "077bda39-7c8c-4b20-9036-53487f512299",
        "portion_ids": "adf04f1f-ad12-431e-9b7e-0afa49c3400c",
        "submitter_portion_ids": "TCGA-A3-A6NN-10A-01",
        "submitter_aliquot_ids": "TCGA-A3-A6NN-10A-01D-A33K-10",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-A3-A6NN-10A-01D",
        "analyte_ids": "f8a21ea5-10af-4a0b-809b-52c6053e842e",
        "submitter_id": "TCGA-A3-A6NN",
        "case_id": "060b2104-a015-400e-86c6-6febfb02bbd3",
        "slide_ids": "6a83abea-7ace-4df3-83a6-cf4f0dd05850",
        "submitter_sample_ids": "TCGA-A3-A6NN-10A",
        "project_id": "TCGA-KIRC",
        "aliquot_ids": "56cad354-34c7-417f-b806-391caa999cb8",
        "submitter_slide_ids": "TCGA-A3-A6NN-01A-01-TS1"
      },
      {
        "sample_ids": "45c9117d-f512-46d5-b330-9c9b12627c31",
        "portion_ids": "b16079d1-d585-48da-acb2-5c398e78400c",
        "submitter_portion_ids": "TCGA-2Z-A9JD-10A-01",
        "submitter_aliquot_ids": "TCGA-2Z-A9JD-10A-01D-A42L-01",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-2Z-A9JD-10A-01D",
        "analyte_ids": "4b4d1bd0-a299-4093-8da8-8c008dfb6ff6",
        "submitter_id": "TCGA-2Z-A9JD",
        "case_id": "cff68090-09df-492b-874c-0caeb29f9361",
        "slide_ids": "4acf7a4c-ce0e-40d3-841e-84b1ea3e55aa",
        "submitter_sample_ids": "TCGA-2Z-A9JD-10A",
        "project_id": "TCGA-KIRP",
        "aliquot_ids": "1a98795b-767e-45ac-9e0f-f50ab97a76df",
        "submitter_slide_ids": "TCGA-2Z-A9JD-01A-01-TS1"
      },
      {
        "sample_ids": "a3e7242d-47a7-44a3-9cf0-357d398a9735",
        "portion_ids": "1a7c879c-1486-4e86-883f-d25ba8c699f7",
        "submitter_portion_ids": "TCGA-EB-A5VU-01A-21",
        "submitter_aliquot_ids": "TCGA-EB-A5VU-01A-21R-A32K-13",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-EB-A5VU-01A-21R",
        "analyte_ids": "230b97aa-76d2-479a-be09-5910c0a10f92",
        "submitter_id": "TCGA-EB-A5VU",
        "case_id": "b6ee9305-c935-4692-81f3-956956d5b4c4",
        "slide_ids": "ed27865d-f0b4-46e8-ae93-a85ae6c6b7c3",
        "submitter_sample_ids": "TCGA-EB-A5VU-01A",
        "project_id": "TCGA-SKCM",
        "aliquot_ids": "d4d3ac8a-9d06-489f-a644-db77770e6dfb",
        "submitter_slide_ids": "TCGA-EB-A5VU-01A-02-TSB"
      },
      {
        "sample_ids": "30f43a40-040f-4c42-a84e-0ba7ec9c8945",
        "portion_ids": "28feac7b-f5ef-4a20-9b46-ca42e51fac15",
        "submitter_portion_ids": "TCGA-KR-A7K8-01A-11",
        "submitter_aliquot_ids": "TCGA-KR-A7K8-01A-11D-A33K-10",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-KR-A7K8-01A-11D",
        "analyte_ids": "32ab0de2-a816-44be-b9e9-f81f83e07dda",
        "submitter_id": "TCGA-KR-A7K8",
        "case_id": "4b0e544c-4090-4ebb-b3e5-2623aa91356c",
        "slide_ids": "c1708f35-0cce-4ac5-9af1-48682b756b24",
        "submitter_sample_ids": "TCGA-KR-A7K8-01A",
        "project_id": "TCGA-LIHC",
        "aliquot_ids": "a4fc5ba8-74ba-4047-8857-5d6d6cbc978b",
        "submitter_slide_ids": "TCGA-KR-A7K8-01A-01-TS1"
      },
      {
        "sample_ids": "3eae6b68-d333-4b36-9e70-d2bed19c1dfa",
        "portion_ids": "25222281-b9bd-49b2-a66e-a1ab930ef2d0",
        "submitter_portion_ids": "TCGA-EJ-5510-01A-01",
        "submitter_aliquot_ids": "TCGA-EJ-5510-01A-01R-1579-13",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-EJ-5510-01A-01R",
        "analyte_ids": "b252b7eb-c317-417d-99c3-c21f19a5512a",
        "submitter_id": "TCGA-EJ-5510",
        "case_id": "4a3d9efe-f24e-4135-8a70-e202651e7a81",
        "slide_ids": "a01a6b12-1019-4b95-8a1e-d2ad1848c9af",
        "submitter_sample_ids": "TCGA-EJ-5510-01A",
        "project_id": "TCGA-PRAD",
        "aliquot_ids": "7c0db1d7-66a9-4e71-8c97-005d9a15cb76",
        "submitter_slide_ids": "TCGA-EJ-5510-01A-01-TS1"
      },
      {
        "sample_ids": "4b13bcc0-d9ff-4616-b6cf-543a86881828",
        "portion_ids": "fc4e2d10-6b8e-447d-9757-0d125bd78b02",
        "submitter_portion_ids": "TCGA-AA-3979-01A-01",
        "submitter_aliquot_ids": "TCGA-AA-3979-01A-01D-1637-02",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-AA-3979-01A-01X",
        "analyte_ids": "d77c659b-5d0e-4ee7-ab70-ab938a271438",
        "submitter_id": "TCGA-AA-3979",
        "case_id": "af17788e-a562-43aa-bd36-fd3e031f483c",
        "slide_ids": "93a2c329-a9fe-44de-9115-82b433434b54",
        "submitter_sample_ids": "TCGA-AA-3979-01A",
        "project_id": "TCGA-COAD",
        "aliquot_ids": "5e52933c-840a-4feb-adbe-9d0cfe69e939",
        "submitter_slide_ids": "TCGA-AA-3979-01A-01-TS1"
      },
      {
        "sample_ids": "ef1ab6b8-af5e-4d24-a035-9e4e636e6be7",
        "portion_ids": "253a641f-1814-4b56-8f48-491ad508f770",
        "submitter_portion_ids": "TCGA-BQ-5877-11A-01",
        "submitter_aliquot_ids": "TCGA-BQ-5877-11A-01R-1591-13",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-BQ-5877-11A-01W",
        "analyte_ids": "4797d59e-1d9b-48b8-9e91-779d374f941f",
        "submitter_id": "TCGA-BQ-5877",
        "case_id": "53d40ba3-c33b-4139-967d-5110884237a3",
        "slide_ids": "e769ffb2-da9b-4f8a-abb3-4c84d9c71bd8",
        "submitter_sample_ids": "TCGA-BQ-5877-11A",
        "project_id": "TCGA-KIRP",
        "aliquot_ids": "97f70a79-624f-495c-adb5-723399d2249b",
        "submitter_slide_ids": "TCGA-BQ-5877-11A-01-TS1"
      }
    ],
    "pagination": {
      "count": 10,
      "sort": "",
      "from": 1,
      "page": 1,
      "total": 6303,
      "pages": 631,
      "size": 10
    }
  },
  "warnings": {}
}
```

More in depth examples of various filter types supported in GDC are available in the [Appendix A](Appendix_A_Available_Fields.md)

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
    "fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id",
    "size":"100"
}
```
```Shell
curl --request POST --header "Content-Type: application/json" --data @Payload.txt 'https://gdc-api.nci.nih.gov/files' > File_metadata.txt
```
```File_metadata_txt
cases_0_submitter_id	cases_0_case_id	data_type	cases_0_samples_0_sample_type	cases_0_samples_0_tissue_type	file_name	cases_0_samples_0_submitter_id	cases_0_samples_0_tumor_descriptor	file_id	data_category	cases_0_samples_0_sample_id
TCGA-G7-6790	e7a1cbe2-793c-4747-8412-8be794f2382b	Aligned Reads	Blood Derived Normal		C489.TCGA-G7-6790-10A-01D-1962-08.2_gdc_realn.bam	TCGA-G7-6790-10A		004e2a2c-1acc-4873-9379-ef1aa12283b6	Raw Sequencing Data	4be83d0f-8b09-4e9e-8318-358371d34332
TCGA-B0-5117	ae55b2d3-62a1-419e-9f9a-5ddfac356db4	Aligned Reads	Solid Tissue Normal		C345.TCGA-B0-5117-11A-01D-1421-08.5_gdc_realn.bam	TCGA-B0-5117-11A		003143c8-bbbf-46b9-a96f-f58530f4bb82	Raw Sequencing Data	b1116541-bece-4df3-b3dd-cec50aeb277b
TCGA-B0-5094	8aaa4e25-5c12-4ace-96dc-91aaa0c4457c	Aligned Reads	Solid Tissue Normal		C345.TCGA-B0-5094-11A-01D-1421-08.5_gdc_realn.bam	TCGA-B0-5094-11A		0001801b-54b0-4551-8d7a-d66fb59429bf	Raw Sequencing Data	7519d7a8-c3ee-417b-9cfc-111bc5ad0637
TCGA-B9-A69E	a4225cb2-7b4b-4122-b6b9-629c26e3ea56	Aligned Reads	Blood Derived Normal		TCGA-B9-A69E-10A-01D-A31X-10_Illumina_gdc_realn.bam	TCGA-B9-A69E-10A		0084a614-780b-42ec-b85f-7a1b83128cd3	Raw Sequencing Data	5d6d6cd4-6a7b-499d-936a-1be9bf74b07f
TCGA-EE-A2GU	24faa36a-268d-4a13-b3ae-eacd431a2bcc	Aligned Reads	Blood Derived Normal		C828.TCGA-EE-A2GU-10A-01D-A198-08.2_gdc_realn.bam	TCGA-EE-A2GU-10A		00c74a8b-10aa-40cc-991e-3365ea1f3fce	Raw Sequencing Data	cc4a5ed8-376a-4842-a25d-ffb07d8e1ca0
TCGA-CE-A484	e62a728d-390f-428a-bea1-fc8c9814fb11	Aligned Reads	Blood Derived Normal		C499.TCGA-CE-A484-10A-01D-A23U-08.3_gdc_realn.bam	TCGA-CE-A484-10A		00df5a50-bce3-4edf-a078-641e54800dcb	Raw Sequencing Data	27a8008e-044a-4966-b518-cc6905e292ca
TCGA-DA-A1IB	8fc9cc74-f388-49f0-b957-debb62638634	Aligned Reads	Blood Derived Normal		C828.TCGA-DA-A1IB-10A-01D-A198-08.2_gdc_realn.bam	TCGA-DA-A1IB-10A		00ab2b5a-b59e-4ec9-b297-76f74ff1d3fb	Raw Sequencing Data	432952c5-6505-4220-a581-f65270a45281
TCGA-AX-A2HG	7a2cf5ce-8317-4fff-946e-b9937afab815	Aligned Reads	Blood Derived Normal		6c2a8ea343da8d6cc0fd2043492f16df_gdc_realn.bam	TCGA-AX-A2HG-10A		00c5f14e-a398-4076-95d1-25f320ee3a37	Raw Sequencing Data	ef4b80ec-b453-48ec-8ad8-ccac83e1e4db
TCGA-EC-A24G	b5c1e511-baf2-45b3-9919-110e8941e3c2	Aligned Reads	Blood Derived Normal		671333b193812fc2bd2744053b383459_gdc_realn.bam	TCGA-EC-A24G-10A		00a5e471-a79f-4d56-8a4c-4847ac037400	Raw Sequencing Data	61cf2e54-1b8d-40a0-9c73-a7449cbd570a
TCGA-B5-A0K0	29c8f468-5ac1-4d6c-8376-e36e6d246926	Aligned Reads	Blood Derived Normal		TCGA-B5-A0K0-10A-01W-A062-09_IlluminaGA-DNASeq_exome_gdc_realn.bam	TCGA-B5-A0K0-10A		007ce9b5-3268-441e-9ffd-b40d1127a319	Raw Sequencing Data	1df69e2e-f392-465f-8e61-4671ba2fcd35
TCGA-C8-A27B	f0d8a1fe-e313-44f1-99cc-b965cbeeff0e	Aligned Reads	Blood Derived Normal		3c99d98ea8eb6acbf819e67fc77623d9_gdc_realn.bam	TCGA-C8-A27B-10A		006b8839-31e5-4697-b912-8e3f4124dd15	Raw Sequencing Data	31139082-7978-45aa-9d8f-ac4789ac5cec
TCGA-E9-A295	fec0da58-1047-44d2-b6d1-c18cceed43dc	Aligned Reads	Blood Derived Normal		fd4421a6bbf3efd4e3d5c17fdd610314_gdc_realn.bam	TCGA-E9-A295-10A		002c67f2-ff52-4246-9d65-a3f69df6789e	Raw Sequencing Data	e74183e1-f0b4-412a-8dac-a62d404add78
TCGA-EB-A44O	c787c4da-c564-44f1-89eb-dd9da107acb1	Aligned Reads	Blood Derived Normal		C828.TCGA-EB-A44O-10A-01D-A25O-08.3_gdc_realn.bam	TCGA-EB-A44O-10A		006ce9a8-cf38-462e-bb99-7f08499244ab	Raw Sequencing Data	5b738547-1825-4684-81bd-864bf2eb43ef
TCGA-A2-A3XX	53886143-c1c6-40e9-88e6-e4e5e0271fc8	Aligned Reads	Blood Derived Normal		b40998d4778f18ed80d6dd8bff0eb761_gdc_realn.bam	TCGA-A2-A3XX-10A		0043d981-3c6b-463f-b512-ab1d076d3e62	Raw Sequencing Data	c6eb6218-ad71-40a6-88b7-a4f1a015b816
TCGA-EB-A3XB	a9255dcb-b236-4777-ac43-555e3a5386c3	Aligned Reads	Blood Derived Normal		C828.TCGA-EB-A3XB-10B-01D-A23B-08.1_gdc_realn.bam	TCGA-EB-A3XB-10B		005239a8-2e63-4ff1-9cd4-714f81837a61	Raw Sequencing Data	0e1d4c7c-204d-4765-b090-68ed4cd83835
```


The GDC API also supports requests with `Content-Type: application/x-www-form-urlencoded` (curl default), which require POST payloads in the following format:

	filters=%7B%0A%20%20%20%20%22op%22%3A%22in%22%2C%0A%20%20%20%20%22content%22%3A%7B%0A%20%20%20%20%20%20%20%20%22field%22%3A%22files.file_id%22%2C%0A%20%20%20%20%20%20%20%20%22value%22%3A%5B%0A%20%20%20%20%20%20%20%20%20%20%20%20%220001801b-54b0-4551-8d7a-d66fb59429bf%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22002c67f2-ff52-4246-9d65-a3f69df6789e%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22003143c8-bbbf-46b9-a96f-f58530f4bb82%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%220043d981-3c6b-463f-b512-ab1d076d3e62%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22004e2a2c-1acc-4873-9379-ef1aa12283b6%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22005239a8-2e63-4ff1-9cd4-714f81837a61%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22006b8839-31e5-4697-b912-8e3f4124dd15%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22006ce9a8-cf38-462e-bb99-7f08499244ab%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22007ce9b5-3268-441e-9ffd-b40d1127a319%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%220084a614-780b-42ec-b85f-7a1b83128cd3%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%2200a5e471-a79f-4d56-8a4c-4847ac037400%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%2200ab2b5a-b59e-4ec9-b297-76f74ff1d3fb%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%2200c5f14e-a398-4076-95d1-25f320ee3a37%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%2200c74a8b-10aa-40cc-991e-3365ea1f3fce%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%2200df5a50-bce3-4edf-a078-641e54800dcb%22%0A%20%20%20%20%20%20%20%20%5D%0A%20%20%20%20%7D%0A%7D&fields=file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id&format=tsv&size=100


### Facets
The `facets` query parameter provides aggregated data based on a search query. The primary intended use of this parameter is for displaying aggregate information in the GDC Data Portal. For example, to get a count of projects in each program, `facets=program.name` can be passed to the `projects` endpoint.

#### Example

```shell
curl  'https://gdc-api.nci.nih.gov/projects?facets=program.name&from=1&size=0&sort=program.name:asc&pretty=true'
```
```python
import requests
import json

projects_endpt = 'https://gdc-api.nci.nih.gov/projects'
params = {'facets':'program.name',
          'from':1, 'size':0,
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
		  "from": 1,
		  "pages": 46,
		  "total": 46,
		  "page": 1,
		  "size": 0
		},
		"hits": [],
		"aggregations": {
		  "program.name": {
			"buckets": [
			  {
				"key": "TCGA",
				"doc_count": 37
			  },
			  {
				"key": "TARGET",
				"doc_count": 9
			  }
			]
		  }
		}
	  },
	  "warnings": {}
	}
```
