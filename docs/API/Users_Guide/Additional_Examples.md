# Additional Examples

## Data Search and Retrieval

### Endpoint Examples

This section contains additional examples for using endpoints.

#### Project Endpoint Example

This example is a query for Projects contained in GDC. It returns only the first five projects sorted by project name.

```Query
curl 'https://api.gdc.cancer.gov/projects?from=0&size=5&sort=project.name:asc&pretty=true'
```
```Response
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
      "from": 0,
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

This example is a query for files contained in GDC. It returns only the first two files, sorted by file size, from smallest to largest.

``` Query
curl 'https://api.gdc.cancer.gov/files?from=0&size=2&sort=file_size:asc&pretty=true'
```
```Response
{
  "data": {
    "hits": [
      {
        "data_type": "Raw Simple Somatic Mutation",
        "updated_datetime": "2017-03-04T16:45:40.925270-06:00",
        "file_name": "9f78a291-2d50-472c-8f56-5f8fbd09ab2a.snp.Somatic.hc.vcf.gz",
        "submitter_id": "TCGA-13-0757-01A-01W-0371-08_TCGA-13-0757-10A-01W-0371-08_varscan",
        "file_id": "9f78a291-2d50-472c-8f56-5f8fbd09ab2a",
        "file_size": 1120,
        "id": "9f78a291-2d50-472c-8f56-5f8fbd09ab2a",
        "created_datetime": "2016-05-04T14:50:54.560567-05:00",
        "md5sum": "13c1ceb3519615e2c67128b350365fbf",
        "data_format": "VCF",
        "acl": [
          "phs000178"
        ],
        "access": "controlled",
        "state": "live",
        "data_category": "Simple Nucleotide Variation",
        "type": "simple_somatic_mutation",
        "file_state": "submitted",
        "experimental_strategy": "WXS"
      },
      {
        "data_type": "Raw Simple Somatic Mutation",
        "updated_datetime": "2017-03-04T16:45:40.925270-06:00",
        "file_name": "7780009b-abb6-460b-903d-accdac626c2e.snp.Somatic.hc.vcf.gz",
        "submitter_id": "TCGA-HC-8261-01A-11D-2260-08_TCGA-HC-8261-10A-01D-2260-08_varscan",
        "file_id": "7780009b-abb6-460b-903d-accdac626c2e",
        "file_size": 1237,
        "id": "7780009b-abb6-460b-903d-accdac626c2e",
        "created_datetime": "2016-05-08T13:54:38.369393-05:00",
        "md5sum": "fd9bb46c8022b96af730c48dc00e2c41",
        "data_format": "VCF",
        "acl": [
          "phs000178"
        ],
        "access": "controlled",
        "state": "live",
        "data_category": "Simple Nucleotide Variation",
        "type": "simple_somatic_mutation",
        "file_state": "submitted",
        "experimental_strategy": "WXS"
      }
    ],
    "pagination": {
      "count": 2,
      "sort": "file_size:asc",
      "from": 0,
      "page": 1,
      "total": 274724,
      "pages": 137362,
      "size": 2
    }
  },
  "warnings": {}
}
```

#### Cases Endpoint Example

This example is a query for cases contained in GDC. It returns only the first five files.

```Query
curl 'https://api.gdc.cancer.gov/cases?from=0&size=5&pretty=true'
```
```Response
{
  "data": {
    "hits": [
      {
        "updated_datetime": "2017-03-09T10:01:14.834935-06:00",
        "submitter_analyte_ids": [
          "TCGA-ER-A193-06A-12D",
          "TCGA-ER-A193-06A-12R",
          "TCGA-ER-A193-06A-12W",
          "TCGA-ER-A193-10A-01W",
          "TCGA-ER-A193-10A-01D"
        ],
        "analyte_ids": [
          "62e14ca4-95f5-4af3-848f-83f7273c3b70",
          "6178b8aa-6afb-4951-bc92-bf9bfc57b9c7",
          "e16b701c-7809-4fb5-a9e0-4ff71e5d1d84",
          "5bfa8c9f-6797-4b2b-9122-854f8ab3bbba",
          "9b73d64e-c973-45b6-be31-a486fb8d1708"
        ],
        "submitter_id": "TCGA-ER-A193",
        "case_id": "8ab09143-daf6-40a9-85d3-0fe9de7b3e06",
        "id": "8ab09143-daf6-40a9-85d3-0fe9de7b3e06",
        "disease_type": "Skin Cutaneous Melanoma",
        "sample_ids": [
          "378b3d8a-adbb-4912-a0bf-6b74a282113e",
          "7a384d44-8b05-4197-9921-7d020ada2437"
        ],
        "portion_ids": [
          "6680bbf2-9cf1-4f93-9ec3-04318cffb5ba",
          "690d3b12-a61d-42fd-af2a-5a7a9a3e5de8",
          "824d724e-6836-423e-a751-fee3260ef4d2"
        ],
        "submitter_portion_ids": [
          "TCGA-ER-A193-06A-21-A20N-20",
          "TCGA-ER-A193-10A-01",
          "TCGA-ER-A193-06A-12"
        ],
        "created_datetime": null,
        "slide_ids": [
          "d2751354-a8b7-4f7a-a4f1-d062de5ceb14"
        ],
        "state": "live",
        "aliquot_ids": [
          "dc9f9544-6c76-4b45-b5c3-dd2fecd5acfe",
          "390b3574-ba23-4ecb-acf8-f5ad8a958bd2",
          "33f43961-b32d-46fc-ba11-264f1101e78d",
          "cd17367c-3270-42ae-8ac5-941a3453ea33",
          "b17269a2-79aa-459e-9c3d-589b7efe6fd9",
          "28a7d729-7555-4545-924b-3dec49b54230",
          "13256e77-0b0b-49e3-9959-3b6730d68732",
          "87ca642a-dd4c-47ea-b81f-2d3402f2157a",
          "8a1bfe0e-c97a-41c4-815f-cf5bb5cfc69f",
          "5e1e9c82-99fd-49de-9dfb-a349d4d8ac94",
          "67f00459-e423-4900-be23-9283b0478620",
          "d939c477-a01f-4d54-bcfb-c9fdd957f2ec"
        ],
        "primary_site": "Skin",
        "submitter_aliquot_ids": [
          "TCGA-ER-A193-06A-12D-A18Y-02",
          "TCGA-ER-A193-10A-01D-A193-01",
          "TCGA-ER-A193-10A-01D-A190-02",
          "TCGA-ER-A193-06A-12D-A197-08",
          "TCGA-ER-A193-06A-12R-A18S-07",
          "TCGA-ER-A193-06A-12W-A20H-08",
          "TCGA-ER-A193-10A-01D-A199-08",
          "TCGA-ER-A193-10A-01D-A38R-08",
          "TCGA-ER-A193-10A-01W-A20J-08",
          "TCGA-ER-A193-06A-12R-A18V-13",
          "TCGA-ER-A193-06A-12D-A19C-05",
          "TCGA-ER-A193-06A-12D-A191-01"
        ],
        "submitter_sample_ids": [
          "TCGA-ER-A193-10A",
          "TCGA-ER-A193-06A"
        ],
        "submitter_slide_ids": [
          "TCGA-ER-A193-06A-01-TSA"
        ]
      },
      {
        "updated_datetime": "2017-03-04T16:39:19.244769-06:00",
        "submitter_analyte_ids": [
          "TCGA-VR-AA4G-10A-01W",
          "TCGA-VR-AA4G-01A-11R",
          "TCGA-VR-AA4G-10A-01D",
          "TCGA-VR-AA4G-01A-11D",
          "TCGA-VR-AA4G-01A-11W"
        ],
        "analyte_ids": [
          "152d7d7a-c746-4b58-8c3f-4252454c7b7c",
          "9090d556-bd2e-4851-8a0c-46e22cc61408",
          "7118f4c3-b635-4428-8240-8db85281f2d9",
          "1d8223ff-685a-4427-a3d1-f53887f2a19d",
          "60dfb30a-bea0-426d-b11d-d5813ba39cfc"
        ],
        "submitter_id": "TCGA-VR-AA4G",
        "case_id": "df5bd25c-d70b-4126-89cb-6c838044ae3b",
        "id": "df5bd25c-d70b-4126-89cb-6c838044ae3b",
        "disease_type": "Esophageal Carcinoma",
        "sample_ids": [
          "21456849-38a9-4190-9ece-ed69b3c24fda",
          "6ee6d239-2af6-41cd-bc32-c5cdaf7742b0"
        ],
        "portion_ids": [
          "484b40d5-d77c-4e6f-9e80-1ef27ffbc8a5",
          "fdc56e67-52ab-44fd-823a-5a3124876ff7"
        ],
        "submitter_portion_ids": [
          "TCGA-VR-AA4G-10A-01",
          "TCGA-VR-AA4G-01A-11"
        ],
        "created_datetime": null,
        "slide_ids": [
          "e950eba2-7d6e-4ffd-a2d5-e0eb6486848a"
        ],
        "state": "live",
        "aliquot_ids": [
          "db6beed3-a5a2-469f-8dc8-00d838c1f37f",
          "f5db4d36-034b-429b-a7be-26a872b702ee",
          "16421a96-b843-4f7e-9f7c-64d2fb5b2a25",
          "5d938cb5-7064-40bc-877d-57faa94c3333",
          "d231404d-ece5-43c0-a8a3-e9f294ceb777",
          "8c77dc3e-2ea3-4626-88f5-e74f242bedf3",
          "993624d4-1c28-41a5-a0b6-094a0e442c36",
          "105a18c9-df7e-4573-b1a2-6a987e57d553",
          "af81c3bb-3b9e-41cb-b85a-b55c6437d05b",
          "38938066-5fd9-415c-b00e-65efff14085e",
          "20139afe-ad04-4571-b779-0c4a51e74ada"
        ],
        "primary_site": "Esophagus",
        "submitter_aliquot_ids": [
          "TCGA-VR-AA4G-10A-01W-A44M-09",
          "TCGA-VR-AA4G-01A-11D-A37B-01",
          "TCGA-VR-AA4G-01A-11D-A37D-05",
          "TCGA-VR-AA4G-10A-01D-A37F-09",
          "TCGA-VR-AA4G-01A-11D-A37R-26",
          "TCGA-VR-AA4G-01A-11R-A37J-13",
          "TCGA-VR-AA4G-01A-11R-A37I-31",
          "TCGA-VR-AA4G-01A-11D-A37C-09",
          "TCGA-VR-AA4G-10A-01D-A37R-26",
          "TCGA-VR-AA4G-10A-01D-A37E-01",
          "TCGA-VR-AA4G-01A-11W-A44L-09"
        ],
        "submitter_sample_ids": [
          "TCGA-VR-AA4G-01A",
          "TCGA-VR-AA4G-10A"
        ],
        "submitter_slide_ids": [
          "TCGA-VR-AA4G-01A-01-TS1"
        ]
      },
      {
        "updated_datetime": "2017-03-04T16:39:19.244769-06:00",
        "submitter_analyte_ids": [
          "TCGA-D1-A174-01A-11D",
          "TCGA-D1-A174-01A-11W",
          "TCGA-D1-A174-10A-01D",
          "TCGA-D1-A174-10A-01W",
          "TCGA-D1-A174-01A-11R"
        ],
        "analyte_ids": [
          "96203028-f824-4a90-9758-22340285062c",
          "f4878e33-b773-43b5-83a5-9fd8e539e668",
          "8627ccd0-0575-4d03-b589-ca45642d523d",
          "1183f7c6-992d-4084-946e-adce7c52f9cc",
          "5343f6a8-8ac2-4446-ace5-a27d21e76844"
        ],
        "submitter_id": "TCGA-D1-A174",
        "case_id": "fc7315b0-9f48-4206-b197-2268c0518eb4",
        "id": "fc7315b0-9f48-4206-b197-2268c0518eb4",
        "disease_type": "Uterine Corpus Endometrial Carcinoma",
        "sample_ids": [
          "df9a1f44-9b3f-48b2-96af-54aaabdfd243",
          "ad5a9cb6-b3f9-4651-b6d1-13c78010bd88"
        ],
        "portion_ids": [
          "79dd516c-bae3-4f6e-b4cb-901de030acb7",
          "6e55e6d9-902f-439b-b6f1-ca296c123fd3"
        ],
        "submitter_portion_ids": [
          "TCGA-D1-A174-01A-11",
          "TCGA-D1-A174-10A-01"
        ],
        "created_datetime": null,
        "slide_ids": [
          "7602727e-b46d-40fc-bd03-5ccf631041f8"
        ],
        "state": "live",
        "aliquot_ids": [
          "5c15542b-cd63-44b5-b278-e211410fb0aa",
          "d661cfb9-248a-49e6-b0db-865ca257e8dc",
          "83bd3bdb-9bd3-46fa-888c-f6f5efec530f",
          "c46551c9-c0d0-4140-8d0a-946b53e504e2",
          "96b511df-3a69-4168-908c-662060b4f976",
          "0182d4e1-f835-46b5-a8f0-53decf5868de",
          "e9563a06-0b86-4986-976e-43d4040f1d61",
          "6bb2de6e-5b85-4e97-a930-1f2c6bf663a1",
          "f6ee5558-a1b6-4b11-8f48-c17186fff39a",
          "67f6f0d9-6581-4946-a9c7-a6629da86888",
          "39e9a948-054a-4b50-b108-7d7aee686363",
          "ddb4ca26-655d-4bdc-a00d-7caf26cadafe"
        ],
        "primary_site": "Uterus",
        "submitter_aliquot_ids": [
          "TCGA-D1-A174-01A-11D-A12F-02",
          "TCGA-D1-A174-01A-01D-YYYY-23",
          "TCGA-D1-A174-01A-11W-A139-09",
          "TCGA-D1-A174-10A-01W-A139-09",
          "TCGA-D1-A174-01A-11D-A12K-05",
          "TCGA-D1-A174-10A-01D-A12F-02",
          "TCGA-D1-A174-10A-01D-A12G-01",
          "TCGA-D1-A174-01A-11R-A12I-07",
          "TCGA-D1-A174-01A-11D-A12J-09",
          "TCGA-D1-A174-10A-01D-A12J-09",
          "TCGA-D1-A174-01A-11R-A12H-13",
          "TCGA-D1-A174-01A-11D-A12G-01"
        ],
        "submitter_sample_ids": [
          "TCGA-D1-A174-01A",
          "TCGA-D1-A174-10A"
        ],
        "submitter_slide_ids": [
          "TCGA-D1-A174-01A-01-TS1"
        ]
      },
      {
        "updated_datetime": "2017-03-04T16:39:19.244769-06:00",
        "submitter_analyte_ids": [
          "TCGA-XM-A8RL-10A-01D",
          "TCGA-XM-A8RL-01A-11R",
          "TCGA-XM-A8RL-01A-11D"
        ],
        "analyte_ids": [
          "2c483e72-92b0-425d-ac1b-b75a169cf531",
          "57f88d4f-8b3a-4349-88b0-3d2e58a95ed9",
          "499bfbe1-639c-479c-abaa-42cbb11c0568"
        ],
        "submitter_id": "TCGA-XM-A8RL",
        "case_id": "dd240b82-b1d6-4c0f-aa3e-6fcfe1364ec1",
        "id": "dd240b82-b1d6-4c0f-aa3e-6fcfe1364ec1",
        "disease_type": "Thymoma",
        "sample_ids": [
          "cb091cc1-7bbe-43a4-8460-01215af3aa21",
          "cabc9729-c1e1-4f08-9959-985dcb7a00d5"
        ],
        "portion_ids": [
          "e8ea57c9-729e-46ea-b1da-2db7a00b02bc",
          "8e2edb92-753f-4cb0-a5b8-8c45dbefaf36",
          "650fa4f2-9fa2-4d3a-8b63-ff4a9bd8c33e"
        ],
        "submitter_portion_ids": [
          "TCGA-XM-A8RL-01A-21-A45R-20",
          "TCGA-XM-A8RL-10A-01",
          "TCGA-XM-A8RL-01A-11"
        ],
        "created_datetime": null,
        "slide_ids": [
          "08cedd34-aafd-4b47-891f-cf66ee1f627b"
        ],
        "state": "live",
        "aliquot_ids": [
          "df9d8553-8d5b-4c65-8b28-74030a8f8e76",
          "47b7f634-b36f-49e9-a4dc-d8f5508fdc0a",
          "e692ebed-9721-40db-8986-fcaba07d68f1",
          "189ee080-95d1-4ccb-8618-955605c7bd55",
          "83af7ff3-45be-4378-a8b5-5dff3584e95d",
          "42ebb1f0-e236-48ae-847f-69a153969903",
          "e8a4938f-6b93-4ad1-9324-31c97dd1d477"
        ],
        "primary_site": "Thymus",
        "submitter_aliquot_ids": [
          "TCGA-XM-A8RL-10A-01D-A426-09",
          "TCGA-XM-A8RL-01A-11D-A423-09",
          "TCGA-XM-A8RL-01A-11D-A422-01",
          "TCGA-XM-A8RL-01A-11R-A42C-07",
          "TCGA-XM-A8RL-10A-01D-A425-01",
          "TCGA-XM-A8RL-01A-11R-A42W-13",
          "TCGA-XM-A8RL-01A-11D-A424-05"
        ],
        "submitter_sample_ids": [
          "TCGA-XM-A8RL-10A",
          "TCGA-XM-A8RL-01A"
        ],
        "submitter_slide_ids": [
          "TCGA-XM-A8RL-01A-01-TSA"
        ]
      },
      {
        "updated_datetime": "2017-03-04T16:39:19.244769-06:00",
        "submitter_analyte_ids": [
          "TCGA-B0-5120-01A-01W",
          "TCGA-B0-5120-01A-01D",
          "TCGA-B0-5120-01A-01R",
          "TCGA-B0-5120-11A-01W",
          "TCGA-B0-5120-11A-01D"
        ],
        "analyte_ids": [
          "996336e6-fad7-4100-96ae-60adb5c276f1",
          "0eb7da02-0b90-4f6d-abd2-b048a9cb2995",
          "fa2861b9-67c1-486a-a1e0-95d8f8adf65b",
          "7e9f5639-a462-493e-98f8-1b7aeee383c7",
          "d51e9fd4-0c99-49ec-9de5-db3946b0bf43"
        ],
        "submitter_id": "TCGA-B0-5120",
        "case_id": "c5bf474c-6919-47b4-ba59-34ab20c087d5",
        "id": "c5bf474c-6919-47b4-ba59-34ab20c087d5",
        "disease_type": "Kidney Renal Clear Cell Carcinoma",
        "sample_ids": [
          "b50d3c6f-fdec-488b-ab26-a9b690fad34f",
          "f3148210-ecae-4314-b5f8-9bee2315a093"
        ],
        "portion_ids": [
          "b8fcbf00-4c5a-42c3-95e9-fb6e169a8da9",
          "34443e91-0210-4477-9511-53026ae62b38",
          "e466f011-79a1-4158-b796-f8e9dda32d68"
        ],
        "submitter_portion_ids": [
          "TCGA-B0-5120-01A-01",
          "TCGA-B0-5120-11A-01",
          "TCGA-B0-5120-01A-21-1740-20"
        ],
        "created_datetime": null,
        "slide_ids": [
          "e5a29e92-4125-4acb-a797-86822b4961a2",
          "78d873e0-037f-4aef-8725-7c651598b1f8",
          "43d8cec7-f5a0-45d5-a5f8-cc77d6b7b539"
        ],
        "state": "live",
        "aliquot_ids": [
          "b35280fe-dbfa-4e45-8f49-3d0489e68743",
          "a2e3a2f2-c32b-44a1-9b29-911145d700b8",
          "a064d108-e8b2-46fa-b277-0a7a89904a3a",
          "59be71a1-50e3-4565-852a-173afc8a6851",
          "136dff0e-b181-49c9-8305-b3289625ea2e",
          "8fbb983b-53ad-44a9-976a-7945628eaa51",
          "cecf40f8-7301-4db9-b276-a14317d4dd59",
          "fac8b066-bf2c-4f08-b42b-251035596a28",
          "fa55c92f-54e8-436b-b8c4-04cb68a24e93",
          "007e3098-aaf9-4ee7-9ae1-f94b131a5ae0",
          "6ce58fbc-6742-4ade-84b0-cd025266e030",
          "9668e15e-a3fa-4ead-ad42-322c5700e0db",
          "c1167003-0730-41d5-bdd5-1cbf501c1463",
          "73aab074-cbd1-45f2-8266-9ef6f7c559bc"
        ],
        "primary_site": "Kidney",
        "submitter_aliquot_ids": [
          "TCGA-B0-5120-11A-01D-1416-02",
          "TCGA-B0-5120-11A-01D-2099-10",
          "TCGA-B0-5120-11A-01D-1418-05",
          "TCGA-B0-5120-01A-01W-1475-10",
          "TCGA-B0-5120-01A-01D-1421-08",
          "TCGA-B0-5120-01A-01D-1416-02",
          "TCGA-B0-5120-01A-01R-1419-13",
          "TCGA-B0-5120-01A-01R-1420-07",
          "TCGA-B0-5120-11A-01D-1421-08",
          "TCGA-B0-5120-01A-01D-1417-01",
          "TCGA-B0-5120-01A-01D-1418-05",
          "TCGA-B0-5120-11A-01W-1475-10",
          "TCGA-B0-5120-01A-01D-2099-10",
          "TCGA-B0-5120-11A-01D-1417-01"
        ],
        "submitter_sample_ids": [
          "TCGA-B0-5120-11A",
          "TCGA-B0-5120-01A"
        ],
        "submitter_slide_ids": [
          "TCGA-B0-5120-11A-01-TS1",
          "TCGA-B0-5120-01A-01-BS1",
          "TCGA-B0-5120-01A-01-TS1"
        ]
      }
    ],
    "pagination": {
      "count": 5,
      "sort": "",
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

#### Annotations Endpoint Example

This example is a query for annotations contained in the GDC. It returns only the first two annotations.

```Query
curl 'https://api.gdc.cancer.gov/annotations?from=0&size=2&pretty=true'
```
```Response
{
  "data": {
    "hits": [
      {
        "category": "History of unacceptable prior treatment related to a prior/other malignancy",
        "status": "Approved",
        "entity_id": "51c37449-6a2e-4c3d-a7cc-06f901e1224f",
        "classification": "Notification",
        "entity_type": "case",
        "created_datetime": "2014-06-16T00:00:00",
        "annotation_id": "3d086829-de62-5d08-b848-ce0724188ff0",
        "notes": "unknown treatment history",
        "updated_datetime": "2017-03-09T12:32:36.305475-06:00",
        "submitter_id": "20743",
        "state": "submitted",
        "case_id": "51c37449-6a2e-4c3d-a7cc-06f901e1224f",
        "case_submitter_id": "TCGA-AG-A014",
        "entity_submitter_id": "TCGA-AG-A014",
        "id": "3d086829-de62-5d08-b848-ce0724188ff0"
      },
      {
        "category": "Center QC failed",
        "status": "Approved",
        "entity_id": "733f0607-6c6b-4385-9868-fa6f155a9a2e",
        "classification": "CenterNotification",
        "entity_type": "aliquot",
        "created_datetime": "2012-07-20T00:00:00",
        "annotation_id": "5cf05f41-ce70-58a3-8ecb-6bfaf6264437",
        "notes": "RNA-seq:INSUFFICIENT INPUT MATERIAL,LOW SEQUENCE YIELD/DIVERSITY;LOW 5/3 COVERAGE RATIO",
        "updated_datetime": "2017-03-09T13:51:45.396638-06:00",
        "submitter_id": "8764",
        "state": "submitted",
        "case_id": "3e8a51bf-7e1f-4eab-af83-3c60d04db1bf",
        "case_submitter_id": "TCGA-13-0913",
        "entity_submitter_id": "TCGA-13-0913-02A-01R-1564-13",
        "id": "5cf05f41-ce70-58a3-8ecb-6bfaf6264437"
      }
    ],
    "pagination": {
      "count": 2,
      "sort": "",
      "from": 0,
      "page": 1,
      "total": 2361,
      "pages": 1181,
      "size": 2
    }
  },
  "warnings": {}
}
```

### Filters Examples

This section contains additional examples for using the `filters` parameter.

#### Example: Basic syntax

The following is an example of `filters` syntax, including the JSON object passed to the `filters` parameter, the corresponding API query, and the JSON object returned by the API. The example finds projects where the primary site is Blood.

```Filter
{
  "op": "and",
  "content": [
    {
      "op": "in",
      "content": {
        "field": "primary_site",
        "value": [
          "Blood"
        ]
      }
    }
  ]
}
```
```Query
curl 'https://api.gdc.cancer.gov/projects?filters=%7b%0d%0a++%22op%22%3a+%22and%22%2c%0d%0a++%22content%22%3a+%5b%0d%0a++++%7b%0d%0a++++++%22op%22%3a+%22in%22%2c%0d%0a++++++%22content%22%3a+%7b%0d%0a++++++++%22field%22%3a+%22primary_site%22%2c%0d%0a++++++++%22value%22%3a+%5b%0d%0a++++++++++%22Blood%22%0d%0a++++++++%5d%0d%0a++++++%7d%0d%0a++++%7d%0d%0a++%5d%0d%0a%7d&pretty=true'
```
```Response
{
  "data": {
    "hits": [
      {
        "dbgap_accession_number": "phs000465",
        "disease_type": [
          "Acute Myeloid Leukemia"
        ],
        "released": true,
        "state": "legacy",
        "primary_site": [
          "Blood"
        ],
        "project_id": "TARGET-AML",
        "id": "TARGET-AML",
        "name": "Acute Myeloid Leukemia"
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

#### Example: Filter cases keeping only 'male'

This is an example of a value-based filter:


```Filter
{
   "op" : "=" ,
   "content" : {
       "field" : "cases.demographic.gender" ,
       "value" : [ "male" ]
   }
}
```
```Query
curl 'https://api.gdc.cancer.gov/cases?filters=%7b%0d%0a+++%22op%22+%3a+%22%3d%22+%2c%0d%0a+++%22content%22+%3a+%7b%0d%0a+++++++%22field%22+%3a+%22cases.demographic.gender%22+%2c%0d%0a+++++++%22value%22+%3a+%5b+%22male%22+%5d%0d%0a+++%7d%0d%0a%7d%0d%0a&fields=demographic.gender,case_id&pretty=true'
```

#### Example: Filter using a range

This is an example of filtering for age at diagnosis. The request is for cases where the age at diagnosis is between 40 and 70 years. >**Note:** `age_at_diagnosis` is expressed in days.

```Filter
{
    "op": "and",
    "content": [
        {
            "op": ">=",
            "content": {
                "field": "cases.diagnoses.age_at_diagnosis",
                "value": [
                    14600
                ]
            }
        },
        {
            "op": "<=",
            "content": {
                "field": "cases.diagnoses.age_at_diagnosis",
                "value": [
                    25550
                ]
            }
        }
    ]
}
```
```Query
curl 'https://api.gdc.cancer.gov/cases?filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22%3E%3D%22,%22content%22:%7B%22field%22:%22cases.diagnoses.age_at_diagnosis%22,%22value%22:%5B14600%5D%7D%7D,%7B%22op%22:%22%3C%3D%22,%22content%22:%7B%22field%22:%22cases.diagnoses.age_at_diagnosis%22,%22value%22:%5B25550%5D%7D%7D%5D%7D&fields=diagnoses.age_at_diagnosis,case_id&pretty=true'
```


#### Example: Multiple fields

Filter projects for primary_site being Kidney or Brain and program.name being TCGA

```Filter
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
```Query
curl 'https://api.gdc.cancer.gov/projects?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22primary_site%22%2C%22value%22%3A%5B%22Kidney%22%2C%22Brain%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%5D%7D&pretty=true'
```
