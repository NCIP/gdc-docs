## 8.1 Data Search and Retrieval

### 8.1.1 Endpoint Examples

This section contains additional examples for using endpoints.

#### Status Endpoint Example
<a name="status-endpoint-ex"></a>

#### Project Endpoint Example
<a name="project-endpoint-ex"></a>

This example is a query for Projects contained in GDC. It returns only the first five projects sorted by project name.

Query

```
$ curl 'https://gdc-api.nci.nih.gov/projects?from=1&size=5&sort=project.name:asc&pretty=true'
```
Response

```
{
  "data": {
    "hits": [
      {
        "state": "legacy",
        "project_id": "TARGET-AML",
        "primary_site": "Blood",
        "disease_type": "Acute Myeloid Leukemia",
        "name": "Acute Myeloid Leukemia"
      },
      {
        "state": "legacy",
        "project_id": "TCGA-LAML",
        "primary_site": "Blood",
        "disease_type": "Acute Myeloid Leukemia",
        "name": "Acute Myeloid Leukemia"
      },
      {
        "state": "legacy",
        "project_id": "TARGET-AML-IF",
        "primary_site": "Blood",
        "disease_type": "Acute Myeloid Leukemia Induction Failure",
        "name": "Acute Myeloid Leukemia Induction Failure"
      },
      {
        "state": "legacy",
        "project_id": "TARGET-ALL-P2",
        "primary_site": "Blood",
        "disease_type": "Acute Lymphoblastic Leukemia",
        "name": "Acute Lymphoblastic Leukemia - Phase II"
      },
      {
        "state": "legacy",
        "project_id": "TARGET-ALL-P1",
        "primary_site": "Blood",
        "disease_type": "Acute Lymphoblastic Leukemia",
        "name": "Acute Lymphoblastic Leukemia - Phase I"
      }
    ],
    "pagination": {
      "count": 5,
      "sort": "project.name:asc",
      "from": 1,
      "pages": 10,
      "total": 46,
      "page": 1,
      "size": 5
    }
  },
  "warnings": {}
}
```
#### Files Endpoint Example
<a name="files-endpoint-ex"></a>

This example is a query for files contained in GDC. It returns only the first two files, sorted by file size, from smallest to largest.

Query

```
$ curl 'https://gdc-api.nci.nih.gov/files?from=1&size=2&sort=file_size:asc&pretty=true'
```

Response

```
{
  "data": {
    "hits": [
      {
        "state": "legacy",
        "project_id": "TARGET-AML",
        "primary_site": "Blood",
        "disease_type": "Acute Myeloid Leukemia",
        "name": "Acute Myeloid Leukemia"
      },
      {
        "state": "legacy",
        "project_id": "TCGA-LAML",
        "primary_site": "Blood",
        "disease_type": "Acute Myeloid Leukemia",
        "name": "Acute Myeloid Leukemia"
      },
      {
        "state": "legacy",
        "project_id": "TARGET-AML-IF",
        "primary_site": "Blood",
        "disease_type": "Acute Myeloid Leukemia Induction Failure",
        "name": "Acute Myeloid Leukemia Induction Failure"
      },
      {
        "state": "legacy",
        "project_id": "TARGET-ALL-P2",
        "primary_site": "Blood",
        "disease_type": "Acute Lymphoblastic Leukemia",
        "name": "Acute Lymphoblastic Leukemia - Phase II"
      },
      {
        "state": "legacy",
        "project_id": "TARGET-ALL-P1",
        "primary_site": "Blood",
        "disease_type": "Acute Lymphoblastic Leukemia",
        "name": "Acute Lymphoblastic Leukemia - Phase I"
      }
    ],
    "pagination": {
      "count": 5,
      "sort": "project.name:asc",
      "from": 1,
      "pages": 10,
      "total": 46,
      "page": 1,
      "size": 5
    }
  },
  "warnings": {}
}
```

#### Cases Endpoint Example
<a name="cases-endpoint-ex"></a>

This example is a query for cases contained in GDC. It returns only the first five files.

Query

```
$ curl 'https://gdc-api.nci.nih.gov/cases?from=1&size=5&pretty=true'
```

Response

```
{
  "data": {
    "hits": [
      {
        "sample_ids": "fae164e6-16ed-4547-9872-15d53c79bb45",
        "portion_ids": "0a5fa1fd-aa9b-49d1-8a32-3522271a56e8",
        "submitter_portion_ids": "TCGA-78-7535-10A-01",
        "submitter_aliquot_ids": "TCGA-78-7535-10A-01W-2107-08",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-78-7535-10A-01W",
        "analyte_ids": "14081a57-a8ee-497d-a944-3f24ef8efddb",
        "submitter_id": "TCGA-78-7535",
        "case_id": "46592b7b-6968-42a6-83af-0917c9f4a9a5",
        "submitter_sample_ids": "TCGA-78-7535-10A",
        "aliquot_ids": "22036caf-c6c9-4ad4-8a69-912b8e56aace"
      },
      {
        "sample_ids": "094cf919-3e36-4d9e-9d37-a00ae04736ee",
        "portion_ids": "1a723c9e-ac2e-40fd-b342-6f9fe7795681",
        "submitter_portion_ids": "TCGA-DJ-A2Q9-01A-11-A21M-20",
        "submitter_aliquot_ids": "TCGA-DJ-A2Q9-01A-21R-A18B-13",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-DJ-A2Q9-01A-21R",
        "analyte_ids": "f888f0c4-7f33-4a64-8975-316d88e214b3",
        "submitter_id": "TCGA-DJ-A2Q9",
        "case_id": "061fab24-727a-4551-a205-89eeb9f530ea",
        "submitter_sample_ids": "TCGA-DJ-A2Q9-01A",
        "aliquot_ids": "29c9b306-d3e9-4d09-bc54-21e46f92ad8a"
      },
      {
        "sample_ids": "8a59b137-8e5e-4484-a9a9-65a596a47ef8",
        "portion_ids": "9f590b62-a6ab-489d-92b9-6e7802812a15",
        "submitter_portion_ids": "TCGA-J4-A83I-01A-11",
        "submitter_aliquot_ids": "TCGA-J4-A83I-01A-11W-A447-08",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-J4-A83I-01A-11W",
        "analyte_ids": "f12a0133-c4c3-4240-a013-1ce159cc08f6",
        "submitter_id": "TCGA-J4-A83I",
        "case_id": "5dc7e186-7e01-4a54-8ae8-350dace2297b",
        "submitter_sample_ids": "TCGA-J4-A83I-01A",
        "aliquot_ids": "dffd09c5-965d-45ac-93ca-f4a547d78684"
      },
      {
        "sample_ids": "587bf402-b61b-444d-af40-6f67bb04c323",
        "portion_ids": "b2e05da1-d75d-4075-9704-2efbd7ed51f3",
        "submitter_portion_ids": "TCGA-BG-A0VW-10A-01",
        "submitter_aliquot_ids": "TCGA-BG-A0VW-10A-01D-A122-09",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-BG-A0VW-10A-01W",
        "analyte_ids": "2f3e28c3-60ec-433a-8f13-6466dd68c5ac",
        "submitter_id": "TCGA-BG-A0VW",
        "case_id": "62d71839-4fba-42e9-9929-8d937f0fe287",
        "submitter_sample_ids": "TCGA-BG-A0VW-10A",
        "aliquot_ids": "a9cd1596-89b6-46fd-9864-8a62f47d1f8b"
      },
      {
        "sample_ids": "26a0b6bc-f8aa-45f1-a215-c90e4e840607",
        "portion_ids": "0dfe9858-d9bf-4e75-84f5-cae466ece831",
        "submitter_portion_ids": "TCGA-P4-AAVL-11A-11",
        "submitter_aliquot_ids": "TCGA-P4-AAVL-11A-11D-A42M-10",
        "days_to_index": 0,
        "submitter_analyte_ids": "TCGA-P4-AAVL-11A-11D",
        "analyte_ids": "f310ac58-97b8-4bba-9dab-ae7768185375",
        "submitter_id": "TCGA-P4-AAVL",
        "case_id": "39847790-c951-4f9d-b23c-88c7f44d20a0",
        "submitter_sample_ids": "TCGA-P4-AAVL-11A",
        "aliquot_ids": "fb226698-844d-4a24-86c2-29571549b9bb"
      }
    ],
    "pagination": {
      "count": 5,
      "sort": "",
      "from": 1,
      "pages": 2822,
      "total": 14108,
      "page": 1,
      "size": 5
    }
  },
  "warnings": {}
}
```
#### Annotations Endpoint Example
<a name="annotations-endpoint-ex"></a>

This example is a query for Annotations contained in GDC. It returns only the first two annotations.

Query

```
$ curl 'https://gdc-api.nci.nih.gov/annotations?from=1&size=2&pretty=true'
```

Response

```
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

#### Data Endpoint Examples
<a name="data-endpoint-ex"></a>

#### Submission Endpoint Examples
<a name="submission-endpoint-ex"></a>

### 8.1.2 Filters Examples

This section contains additional examples for using endpoints.

####Basic Syntax

The basic syntax for the query filter is as below:

Query

```
$ curl 'https://gdc-api.nci.nih.gov/projects?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22primary_site%22%2C%22value%22%3A%5B%22Blood%22%5D%7D%7D%5D%7D&pretty=true'
```

**Note:** the above filter string is encoded to be compliant with curl and GET, the non-encode version would be:

```
{"op":"and","content":[{"op":"in","content":{"field":"primary_site","value":["Blood"]}}]}
```

Response

```
{
  "data": {
    "hits": [
      {
        "state": "legacy",
        "project_id": "TARGET-AML",
        "primary_site": "Blood",
        "disease_type": "Acute Myeloid Leukemia",
        "name": "Acute Myeloid Leukemia"
      },
      {
        "state": "legacy",
        "project_id": "TCGA-LAML",
        "primary_site": "Blood",
        "disease_type": "Acute Myeloid Leukemia",
        "name": "Acute Myeloid Leukemia"
      },
      {
        "state": "legacy",
        "project_id": "TARGET-AML-IF",
        "primary_site": "Blood",
        "disease_type": "Acute Myeloid Leukemia Induction Failure",
        "name": "Acute Myeloid Leukemia Induction Failure"
      },
      {
        "state": "legacy",
        "project_id": "TARGET-ALL-P2",
        "primary_site": "Blood",
        "disease_type": "Acute Lymphoblastic Leukemia",
        "name": "Acute Lymphoblastic Leukemia - Phase II"
      },
      {
        "state": "legacy",
        "project_id": "TCGA-LCML",
        "primary_site": "Blood",
        "disease_type": "Chronic Myelogenous Leukemia",
        "name": "Chronic Myelogenous Leukemia"
      },
      {
        "state": "legacy",
        "project_id": "TARGET-ALL-P1",
        "primary_site": "Blood",
        "disease_type": "Acute Lymphoblastic Leukemia",
        "name": "Acute Lymphoblastic Leukemia - Phase I"
      }
    ],
    "pagination": {
      "count": 6,
      "sort": "",
      "from": 1,
      "pages": 1,
      "total": 6,
      "page": 1,
      "size": 10
    }
  },
  "warnings": {}
}
```

**Note:** URL encoding is used in the above example to handle the special JSON characters. After URL decoding it becomes more apparent that this query finds projects where the primary site is Blood.

#### Example: Filter cases keeping only 'male'

Value-based filter

```
{
   "op" : "=" ,
   "content" : {
       "field" : "cases.clinical.gender" ,
       "value" : [ "male" ]
   }
}
```
'Query returning males only with and fields Case_id and gender'

```
$ curl 'https://gdc-api.nci.nih.gov/casesfilters=%7B%0A%20%20%20%22op%22%20%3A%20%22%3D%22%20%2C%0%20%20%20%22content%22%20%3A%20%7B%0A%20%20%20%20%20%20%20%22field%22%20%3A%20%22cases.clinical.gender%22%20%2C%0A%20%20%20%20%20%20%20%22value%22%20%3A%20%5B%20%22male%22%20%5D%0A%20%20%20%7D%0A%7D&fields=clinical.gender,case_id&pretty=true'
```
#### Example: Filter using a range

Range Filter

```
{
    "op": "and",
    "content": [
        {
            "op": ">=",
            "content": {
                "field": "cases.clinical.age_at_diagnosis",
                "value": [
                    14600
                ]
            }
        },
        {
            "op": "<=",
            "content": {
                "field": "cases.clinical.age_at_diagnosis",
                "value": [
                    25550
                ]
            }
        }
    ]
}
```

Query returning case_id with age at diagnosis between 40 and 70 years

```
$ curl 'https://gdc-api.nci.nih.gov/cases?filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22%3E%3D%22,%22content%22:%7B%22field%22:%22cases.clinical.age_at_diagnosis%22,%22value%22:%5B14600%5D%7D%7D,%7B%22op%22:%22%3C%3D%22,%22content%22:%7B%22field%22:%22cases.clinical.age_at_diagnosis%22,%22value%22:%5B25550%5D%7D%7D%5D%7D&fields=clinical.age_at_diagnosis,case_id&pretty=true'
```

**Note:** age_at_diagnosis is expressed in days.

#### Example: Filter all project with primary_site being Kidney or Brain and program.name being TCGA (Multiple fields operation)

Multiple fields operation Filter

```
{
     "op" : "and" ,
     "content" : [{
             "op" : "in" ,
             "content" : {
                 "field" : "primary_site" ,
                 "value" : [
                     "Kidney" ,
                     "Brain"
                 ]
             }
         }, {
             "op" : "in" ,
             "content" : {
                 "field" : "program.name" ,
                 "value" : [
                     "TCGA"
                 ]
             }
         }]
}
```

Query returning case_id with age at diagnosis between 40 and 70 years

```
$ curl 'https://gdc-api.nci.nih.gov/projects?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22primary_site%22%2C%22value%22%3A%5B%22Kidney%22%2C%22Brain%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%5D%7D&pretty=true'
```

## 8.2 File Download

The example below will download two files (2442bd46-fca4-4a03-9d58-e30929da6ed8 and bb8182f8-8cc1-4504-b900-de0ab3ace714).

``` bash
curl  -O https://gdc-api.nci.nih.gov/data/2442bd46-fca4-4a03-9d58-e30929da6ed8,bb8182f8-8cc1-4504-b900-de0ab3ace714
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 10647    0 10647    0     0  14746      0 --:--:-- --:--:-- --:--:-- 14746

```
