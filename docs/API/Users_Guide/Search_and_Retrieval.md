# Search and Retrieval

## Constructing Queries

API queries for selecting, filtering, and sorting GDC data and metadata are constructed using [endpoints](#query-endpoints), [parameters](#query-parameters), and [filtering operators](#filtering-operators). These queries work only on datasets that have been released to the GDC Data Portal.

Data that is in the process of being submitted to GDC and is only available on the GDC Submission Portal cannot be queried using these methods. See [Submission](Submission.md) for information on how data submitters can query their unreleased data using GraphQL.

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
The GDC Project Endpoint `https://gdc-api.nci.nih.gov/projects` provides overall access to all the data served by GDC organized by Project such project(study) name, program,disease, primary site and state.

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

### Cases Endpoint

The GDC Cases Endpoint `https://gdc-api.nci.nih.gov/cases` enables search and retrieval of information related to a specific case, or sample donor.


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
	        "sample_ids": "7f791228-dd77-4ab0-8227-d784a4c7fea1",
	        "portion_ids": "8629bf5a-cdaf-4f6a-90bb-27dd4a7565c5",
	        "submitter_portion_ids": "TCGA-BH-A0EA-01A-21-A13C-20",
	        "submitter_aliquot_ids": "TCGA-BH-A0EA-01A-11D-A10Y-09",
	        "days_to_index": 0,
	        "submitter_analyte_ids": "TCGA-BH-A0EA-01A-11D",
	        "analyte_ids": "66ed0f86-5ca5-4dec-ba76-7ee4dcf31831",
	        "submitter_id": "TCGA-BH-A0EA",
	        "case_id": "1f601832-eee3-48fb-acf5-80c4a454f26e",
	        "slide_ids": "90154ea1-6b76-4445-870e-d531d6fa1239",
	        "submitter_sample_ids": "TCGA-BH-A0EA-01A",
	        "aliquot_ids": "561b8777-801a-49ed-a306-e7dafeb044b6",
	        "submitter_slide_ids": "TCGA-BH-A0EA-01A-01-TSA"
	      }
	    ],
	    "pagination": {
	      "count": 1,
	      "sort": "",
	      "from": 1,
	      "pages": 1,
	      "total": 1,
	      "page": 1,
	      "size": 10
	    }
	  },
	  "warnings": {}
	}
```

### Annotations Endpoint

The GDC Annotation Endpoint `https://gdc-api.nci.nih.gov/annotations` enables search and retrieval of annotations stored in the GDC.


#### Example

This example is a query for Annotations contained in the GDC. It uses the [from](#from), [size](#size), and [pretty](#pretty) parameters, and returns the first two annotations.

```shell
curl 'https://gdc-api.nci.nih.gov/annotations?from=1&size=2&pretty=true'
```
``` Output

	{
	  "data": {
	    "hits": [
	      {
	        "category": "Item flagged DNU",
	        "status": "Approved",
	        "entity_id": "2b61b856-b988-43ca-8dc5-9f97600118ec",
	        "classification": "CenterNotification",
	        "entity_type": "aliquot",
	        "created_datetime": 1294525038,
	        "annotation_id": "7d01080f-e82d-5e58-98a6-910c041ee2b3",
	        "notes": "SDRF in broad.mit.edu_READ.Genome_Wide_SNP_6.mage-tab.1.1003.0 flagged aliquot to be excluded for analysis based on file 'SCENA_p_TCGAb29and30_SNP_N_GenomeWideSNP_6_C04_569122.ismpolish.data.txt'.",
	        "creator": "DCC",
	        "submitter_id": "1099",
	        "case_id": "e7503a51-6647-4cc2-80dd-645d0df4db43",
	        "entity_submitter_id": "TCGA-AG-A008-10A-01D-A003-01"
	      },
	      {
	        "category": "Item flagged DNU",
	        "status": "Approved",
	        "entity_id": "d1f35d46-c6c9-4cff-ad95-e86d88b38b51",
	        "classification": "CenterNotification",
	        "entity_type": "aliquot",
	        "created_datetime": 1414794925,
	        "annotation_id": "c6a9e076-bb56-5dd9-89e7-c340594fa8f7",
	        "notes": "SDRF in broad.mit.edu_COAD.Genome_Wide_SNP_6.mage-tab.1.2010.0 flagged aliquot to be excluded for analysis based on file 'SNORT_p_TCGA_b89_SNP_N_GenomeWideSNP_6_E05_777376.birdseed.data.txt'.",
	        "creator": "DCC",
	        "submitter_id": "23507",
	        "case_id": "57b0f89f-1b75-453e-922c-01cd4d44ca49",
	        "entity_submitter_id": "TCGA-CK-5914-10A-01D-1649-01"
	      }
	    ],
	    "pagination": {
	      "count": 2,
	      "sort": "",
	      "from": 1,
	      "pages": 12296,
	      "total": 24592,
	      "page": 1,
	      "size": 2
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
facets | false | Provides a list document counts for each included facet
fields | false | Query option to specify which fields to include in the response
filters| false | Query option filters specify criteria for the returned response
from   | false | Specifies the first record to return from the set resulting of a query
size | false | determines the number of results to return
sort | false | specifies a field to sort the returned results by sort order: + use asc for ascending order + use des for descending order
pretty | false | Returns response with indentations and line breaks in a human-readable format
format | false | Returns response in XML or TSV format, JSON is default

### Facets
The `facets` query parameter provides aggregated data based on a search query. In the simplest case, a terms facet can return facet counts for various facet values for a specific field.

#### Example

To get a count of projects in each program, `facets=program.name` can be passed to the `projects` endpoint.

A list of all valid _field_ names that can be used as facets is available in [Appendix A](Appendix_A_Available_Fields.md).

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
``` Output

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


### Fields

This query option specifies which fields to include in the response. Using fields can help improve performance. A complete listing of all valid fields for each endpoint type are available in [Appendix A](Appendix_A_Available_Fields.md).

#### Example

To get back only the file names for each file, `fields=file_name` can be passed to the `files` endpoint.

```shell
curl 'https://gdc-api.nci.nih.gov/files?fields=file_name&pretty=true'
```
```python
import requests
import json

files_endpt = 'https://gdc-api.nci.nih.gov/files'
params = {'fields':'file_name'}
response = requests.get(files_endpt, params = params)
print json.dumps(response.json(), indent=2)
```
``` Output
	{
	  "data": {
		"hits": [
		  {
			"file_name": "Collagen_VI-R-V_GBL1112757.tif"
		  },
		  {
			"file_name": "DUNGS_p_TCGA_b84_115_SNP_N_GenomeWideSNP_6_C09_771624.birdseed.data.txt"
		  },
		  {
			"file_name": "unc.edu.ba350477-3c49-4258-8409-f39447687497.2145690.junction_quantification.txt"
		  },
		  {
			"file_name": "SWEDE_p_TCGAb322_23_24_25_26NSP_GenomeWideSNP_6_C06_1364960.hg18.seg.txt"
		  },
		  {
			"file_name": "unc.edu.a696be4d-8576-43a3-8137-30778dff9b89.2445604.junction_quantification.txt"
		  },
		  {
			"file_name": "MSK_252152921979_S01_CGH-v4_10_27Aug08__GCN_V3_A1__CBS_out.txt"
		  },
		  {
			"file_name": "9721366028_R05C01_Red.idat"
		  },
		  {
			"file_name": "28bdce76404c4a5a9a766c6e93f390ed.bam"
		  },
		  {
			"file_name": "PKC-delta_pS664-R-V_GBL9016761.tif"
		  },
		  {
			"file_name": "jhu-usc.edu_LAML.HumanMethylation450.2.lvl-2.TCGA-AB-3011-03A-01D-0742-05.txt"
		  }
		],
		"pagination": {
		  "count": 10,
		  "sort": "",
		  "from": 1,
		  "pages": 54777,
		  "total": 547761,
		  "page": 1,
		  "size": 10
		}
	  },
	  "warnings": {}
	}
```

### Filters

Using the query option `filters` lets users specify criteria for the returned response. The `filters` syntax is a JSON object that contains the query filters that are translatable to Elastic Search JSON-based queries used by the GDC middleware. Users can get a list of available values for a specific field in the filter by making a call to the appropriate API endpoint using the `facets` parameter.


#### Filtering Operators

Operators allow users to define query conditions. These can be used to restrict facet values and then to connect these in logical statements.

Operators can relate an operation to one field (Single Field Operators, e.g. A = B) or multiple fields (e.g. [and (a,b,c,d)].

Operators (**op** in the examples in Section 6.2) can take different values depending of the context and type of data.

| Type | Possible Values |
| --- | --- |
| Single field | =, != , <, <=, =, >, >=, in, is, not, range, exclude |
| Multiple fields | and, or |

When using multiple fields, operator content requires nested data containing additional operators.

#### Example

To get a list of available values for the `clinical.gender` field, users can query the cases endpoint using the `facets` parameter:

```shell
curl 'https://gdc-api.nci.nih.gov/cases?facets=clinical.gender&from=1&size=0&sort=clinical.gender:asc&pretty=true'
```
```python
import requests
import json

cases_endpt = 'https://gdc-api.nci.nih.gov/cases'
params = {'facets':'clinical.gender',
          'from':1, 'size':0,
          'sort':'clinical.gender:asc'}
response = requests.get(cases_endpt, params = params)
print json.dumps(response.json(), indent=2)
```
``` Output
{
	"data": {
	"pagination": {
		"count": 0,
		"sort": "clinical.gender:asc",
		"from": 1,
		"pages": 14052,
		"total": 14052,
		"page": 1,
		"size": 0
	},
	"hits": [],
	"aggregations": {
		"clinical.gender": {
		"buckets": [
			{
			"key": "female",
			"doc_count": 6598
			},
			{
			"key": "male",
			"doc_count": 6303
			},
			{
			"key": "_missing",
			"doc_count": 1151
			}
		]
		}
	}
	},
	"warnings": {}
}
```

#### Nested Operations

Filters support complex nested operations as well as simple queries on a single field. There are different types of operations available for many uses. For more examples see [Additional Examples](Additional_Examples.md).

It is possible to obtain multiple values from multiple fields in one single query, e.g. `(facets=field1,field2)`.

#### Example

This example returns `male` cases.

The JSON object to be passed to the `filter` parameter looks like:

	{"op": "=",
		  "content": {
			  "field": "cases.clinical.gender",
			  "value": ["male"]
		  }
	}

URL-encoding the above JSON object (see [Tools](Getting_Started.md#tools-for-communicating-with-the-gdc-api)) results in the following string:

	%7b%22op%22%3a+%22%3d%22%2c%0d%0a++++++%22content%22%3a+%7b%0d%0a++++++++++%22field%22%3a+%22cases.clinical.gender%22%2c%0d%0a++++++++++%22value%22%3a+%5b%22male%22%5d%0d%0a++++++%7d%0d%0a%7d

The above string can now be passed to the `filters` parameter in an API call:

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


### From

The GDC API uses pagination, the `from` query parameter specifies the first record to return out of the set of results. For example, if there are 20 cases returned from the `case` endpoint, `from` can be set to 10 and results 10-20 will be returned. The `from` parameter can be used in conjunction with the `size` parameter to return a specific subset of results. For more information see [Filters](#filters) examples.


#### Example

Getting 5 file results starting from the 100th result from a set of 500.

``` Shell
curl 'https://gdc-api.nci.nih.gov/files?fields=file_name&from=101&size=5&pretty=true'
```
``` Python
import requests
import json

files_endpt = 'https://gdc-api.nci.nih.gov/files'
params = {'fields':'file_name',
          'from':101, 'size':5}
response = requests.get(files_endpt, params = params)
print json.dumps(response.json(), indent=2)

```
``` Output
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

### Size

The `size` query parameter specifies the number of results to return. When `size` is not specified the default is 10.

#### Example

This example returns the names of the first two files:

``` Shell
curl 'https://gdc-api.nci.nih.gov/files?fields=file_name&from=0&size=2&pretty=true'
```
``` Python
import requests
import json

files_endpt = 'https://gdc-api.nci.nih.gov/files'
params = {'fields':'file_name',
          'from':0, 'size':2}
response = requests.get(files_endpt, params = params)
print json.dumps(response.json(), indent=2)

```
``` Output
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

### Sort

The `sort` query parameter sorts the results by a specific field, and with the sort direction specified using the `:asc` (ascending) or `:dsc` (descending) prefix, e.g. `sort=field:desc`.

#### Example

Sort cases by submitter_id in ascending order:

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

### Pretty

Returns response with indentations and line breaks in a human-readable format.

#### Example: pretty=false

```Shell
curl  'https://gdc-api.nci.nih.gov/cases?fields=submitter_id&sort=submitter_id:asc&size=5'
```
```Output
{"data": {"hits": [{"submitter_id": "TARGET-20-PABGKN"}, {"submitter_id": "TARGET-20-PABHET"}, {"submitter_id": "TARGET-20-PABHKY"}, {"submitter_id": "TARGET-20-PABLDZ"}, {"submitter_id": "TARGET-20-PACDZR"}], "pagination": {"count": 5, "sort": "submitter_id.raw:asc", "from": 1, "pages": 2811, "total": 14052, "page": 1, "size": 5}}, "warnings": {}}
```

####  Example: pretty=true

Same query as above with pretty=true returns a more human-readable format

```Shell
curl  'https://gdc-api.nci.nih.gov/cases?fields=submitter_id&sort=submitter_id:asc&size=5&pretty=true'
```
```Output
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

### Format

Specifies the format of the API response. JSON is default, and `TSV` and `XML` options are available.

#### Example: TSV

A query with `format=TSV` returns results in a tab separated values (TSV) format.

```shell
curl  'https://gdc-api.nci.nih.gov/cases?fields=submitter_id&size=5&format=TSV'
```
```python
import requests

cases_endpt = 'https://gdc-api.nci.nih.gov/cases'
params = {'fields':'submitter_id',
          'format':'TSV'}
response = requests.get(cases_endpt, params = params)
print response.content
```
```Output
submitter_id
TCGA-RC-A6M6
TCGA-B6-A0RV
TCGA-MB-A5Y8
TCGA-BQ-5876
TCGA-Z6-A9VB
```

#### Example: XML

A query with `format=XML` returns results in XML format:

```shell
curl  'https://gdc-api.nci.nih.gov/cases?fields=submitter_id&size=5&format=XML&pretty=true'
```
```python
import requests

cases_endpt = 'https://gdc-api.nci.nih.gov/cases'
params = {'fields':'submitter_id',
          'format':'XML',
          'pretty':'true'}
response = requests.get(cases_endpt, params = params)
print response.content
```
```Output
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
