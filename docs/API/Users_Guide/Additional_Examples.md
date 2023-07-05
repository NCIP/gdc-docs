# Additional Examples

## Data Search and Retrieval

### Endpoint Examples

This section contains additional examples for using endpoints.

#### Project Endpoint Example

This example is a query for Projects contained in GDC. It returns only the first five projects sorted by project name.

```Query
curl 'https://api.gdc.cancer.gov/projects?from=0&size=5&sort=name:asc&pretty=true'
```
```Response
{
  "data": {
    "hits": [
      {
        "id": "GENIE-DFCI",
        "primary_site": [
          "Other and unspecified parts of biliary tract",
          "Thymus",
          "Gallbladder",
          "Bladder",
          "Colon",
          "Small intestine",
          "Prostate gland",
          "Pancreas",
          "Hypopharynx",
          "Peripheral nerves and autonomic nervous system",
          "Rectum",
          "Other and ill-defined sites in lip, oral cavity and pharynx",
          "Adrenal gland",
          "Corpus uteri",
          "Nasal cavity and middle ear",
          "Meninges",
          "Penis",
          "Other and ill-defined digestive organs",
          "Vagina",
          "Stomach",
          "Ovary",
          "Heart, mediastinum, and pleura",
          "Kidney",
          "Anus and anal canal",
          "Nasopharynx",
          "Testis",
          "Brain",
          "Breast",
          "Bronchus and lung",
          "Liver and intrahepatic bile ducts",
          "Other endocrine glands and related structures",
          "Eye and adnexa",
          "Esophagus",
          "Hematopoietic and reticuloendothelial systems",
          "Skin",
          "Other and ill-defined sites",
          "Oropharynx",
          "Other and unspecified female genital organs",
          "Spinal cord, cranial nerves, and other parts of central nervous system",
          "Thyroid gland",
          "Other and unspecified major salivary glands",
          "Other and unspecified urinary organs",
          "Cervix uteri",
          "Larynx",
          "Uterus, NOS",
          "Unknown",
          "Bones, joints and articular cartilage of other and unspecified sites",
          "Retroperitoneum and peritoneum",
          "Connective, subcutaneous and other soft tissues"
        ],
        "dbgap_accession_number": null,
        "project_id": "GENIE-DFCI",
        "disease_type": [
          "Mesothelial Neoplasms",
          "Adenomas and Adenocarcinomas",
          "Blood Vessel Tumors",
          "Fibromatous Neoplasms",
          "Mucoepidermoid Neoplasms",
          "Thymic Epithelial Neoplasms",
          "Trophoblastic neoplasms",
          "Adnexal and Skin Appendage Neoplasms",
          "Miscellaneous Bone Tumors",
          "Immunoproliferative Diseases",
          "Neuroepitheliomatous Neoplasms",
          "Basal Cell Neoplasms",
          "Specialized Gonadal Neoplasms",
          "Mesonephromas",
          "Complex Epithelial Neoplasms",
          "Myeloid Leukemias",
          "Neoplasms of Histiocytes and Accessory Lymphoid Cells",
          "Lipomatous Neoplasms",
          "Nevi and Melanomas",
          "Ductal and Lobular Neoplasms",
          "Granular Cell Tumors and Alveolar Soft Part Sarcomas",
          "Acinar Cell Neoplasms",
          "Other Leukemias",
          "Soft Tissue Tumors and Sarcomas, NOS",
          "Paragangliomas and Glomus Tumors",
          "Hodgkin Lymphoma",
          "Mature B-Cell Lymphomas",
          "Myxomatous Neoplasms",
          "Lymphoid Leukemias",
          "Chronic Myeloproliferative Disorders",
          "Myomatous Neoplasms",
          "Lymphatic Vessel Tumors",
          "Epithelial Neoplasms, NOS",
          "Plasma Cell Tumors",
          "Gliomas",
          "Fibroepithelial Neoplasms",
          "Neoplasms, NOS",
          "Mature T- and NK-Cell Lymphomas",
          "Transitional Cell Papillomas and Carcinomas",
          "Meningiomas",
          "Synovial-like Neoplasms",
          "Malignant Lymphomas, NOS or Diffuse",
          "Myelodysplastic Syndromes",
          "Germ Cell Neoplasms",
          "Cystic, Mucinous and Serous Neoplasms",
          "Squamous Cell Neoplasms",
          "Giant Cell Tumors",
          "Osseous and Chondromatous Neoplasms",
          "Precursor Cell Lymphoblastic Lymphoma",
          "Miscellaneous Tumors",
          "Nerve Sheath Tumors",
          "Complex Mixed and Stromal Neoplasms"
        ],
        "name": "AACR Project GENIE - Contributed by Dana-Farber Cancer Institute",
        "releasable": true,
        "state": "open",
        "released": true
      },
      {
        "id": "GENIE-GRCC",
        "primary_site": [
          "Other and unspecified parts of biliary tract",
          "Thymus",
          "Gallbladder",
          "Bladder",
          "Colon",
          "Small intestine",
          "Prostate gland",
          "Pancreas",
          "Hypopharynx",
          "Peripheral nerves and autonomic nervous system",
          "Rectum",
          "Other and ill-defined sites in lip, oral cavity and pharynx",
          "Adrenal gland",
          "Corpus uteri",
          "Nasal cavity and middle ear",
          "Penis",
          "Other and ill-defined digestive organs",
          "Vagina",
          "Stomach",
          "Ovary",
          "Heart, mediastinum, and pleura",
          "Kidney",
          "Nasopharynx",
          "Anus and anal canal",
          "Testis",
          "Brain",
          "Breast",
          "Bronchus and lung",
          "Liver and intrahepatic bile ducts",
          "Esophagus",
          "Skin",
          "Other and ill-defined sites",
          "Oropharynx",
          "Connective, subcutaneous and other soft tissues",
          "Spinal cord, cranial nerves, and other parts of central nervous system",
          "Other and unspecified female genital organs",
          "Thyroid gland",
          "Other and unspecified major salivary glands",
          "Other and unspecified urinary organs",
          "Cervix uteri",
          "Uterus, NOS",
          "Unknown",
          "Bones, joints and articular cartilage of other and unspecified sites",
          "Retroperitoneum and peritoneum",
          "Larynx"
        ],
        "dbgap_accession_number": null,
        "project_id": "GENIE-GRCC",
        "disease_type": [
          "Mesothelial Neoplasms",
          "Adenomas and Adenocarcinomas",
          "Blood Vessel Tumors",
          "Fibromatous Neoplasms",
          "Mucoepidermoid Neoplasms",
          "Thymic Epithelial Neoplasms",
          "Adnexal and Skin Appendage Neoplasms",
          "Miscellaneous Bone Tumors",
          "Neuroepitheliomatous Neoplasms",
          "Basal Cell Neoplasms",
          "Specialized Gonadal Neoplasms",
          "Complex Epithelial Neoplasms",
          "Odontogenic Tumors",
          "Lipomatous Neoplasms",
          "Nevi and Melanomas",
          "Ductal and Lobular Neoplasms",
          "Granular Cell Tumors and Alveolar Soft Part Sarcomas",
          "Acinar Cell Neoplasms",
          "Soft Tissue Tumors and Sarcomas, NOS",
          "Paragangliomas and Glomus Tumors",
          "Myomatous Neoplasms",
          "Epithelial Neoplasms, NOS",
          "Gliomas",
          "Neoplasms, NOS",
          "Transitional Cell Papillomas and Carcinomas",
          "Synovial-like Neoplasms",
          "Germ Cell Neoplasms",
          "Squamous Cell Neoplasms",
          "Cystic, Mucinous and Serous Neoplasms",
          "Osseous and Chondromatous Neoplasms",
          "Nerve Sheath Tumors",
          "Complex Mixed and Stromal Neoplasms"
        ],
        "name": "AACR Project GENIE - Contributed by Institut Gustave Roussy",
        "releasable": true,
        "state": "open",
        "released": true
      },
      {
        "id": "GENIE-JHU",
        "primary_site": [
          "Thymus",
          "Bladder",
          "Colon",
          "Small intestine",
          "Prostate gland",
          "Pancreas",
          "Peripheral nerves and autonomic nervous system",
          "Rectum",
          "Corpus uteri",
          "Meninges",
          "Other and ill-defined digestive organs",
          "Vagina",
          "Stomach",
          "Ovary",
          "Heart, mediastinum, and pleura",
          "Anus and anal canal",
          "Brain",
          "Breast",
          "Bronchus and lung",
          "Liver and intrahepatic bile ducts",
          "Other endocrine glands and related structures",
          "Eye and adnexa",
          "Esophagus",
          "Hematopoietic and reticuloendothelial systems",
          "Skin",
          "Other and ill-defined sites",
          "Oropharynx",
          "Spinal cord, cranial nerves, and other parts of central nervous system",
          "Thyroid gland",
          "Uterus, NOS",
          "Unknown",
          "Bones, joints and articular cartilage of other and unspecified sites",
          "Connective, subcutaneous and other soft tissues"
        ],
        "dbgap_accession_number": null,
        "project_id": "GENIE-JHU",
        "disease_type": [
          "Mesothelial Neoplasms",
          "Adenomas and Adenocarcinomas",
          "Leukemias, NOS",
          "Thymic Epithelial Neoplasms",
          "Adnexal and Skin Appendage Neoplasms",
          "Neuroepitheliomatous Neoplasms",
          "Complex Epithelial Neoplasms",
          "Myeloid Leukemias",
          "Neoplasms of Histiocytes and Accessory Lymphoid Cells",
          "Lipomatous Neoplasms",
          "Nevi and Melanomas",
          "Ductal and Lobular Neoplasms",
          "Acinar Cell Neoplasms",
          "Other Leukemias",
          "Soft Tissue Tumors and Sarcomas, NOS",
          "Myxomatous Neoplasms",
          "Mature B-Cell Lymphomas",
          "Lymphoid Leukemias",
          "Chronic Myeloproliferative Disorders",
          "Myomatous Neoplasms",
          "Epithelial Neoplasms, NOS",
          "Gliomas",
          "Neoplasms, NOS",
          "Mast Cell Tumors",
          "Transitional Cell Papillomas and Carcinomas",
          "Meningiomas",
          "Synovial-like Neoplasms",
          "Malignant Lymphomas, NOS or Diffuse",
          "Myelodysplastic Syndromes",
          "Cystic, Mucinous and Serous Neoplasms",
          "Squamous Cell Neoplasms",
          "Osseous and Chondromatous Neoplasms",
          "Complex Mixed and Stromal Neoplasms"
        ],
        "name": "AACR Project GENIE - Contributed by Johns Hopkins Sidney Kimmel Comprehensive Cancer Center",
        "releasable": true,
        "state": "open",
        "released": true
      },
      {
        "id": "GENIE-MDA",
        "primary_site": [
          "Other and unspecified parts of biliary tract",
          "Thymus",
          "Gallbladder",
          "Bladder",
          "Colon",
          "Small intestine",
          "Prostate gland",
          "Pancreas",
          "Peripheral nerves and autonomic nervous system",
          "Adrenal gland",
          "Corpus uteri",
          "Meninges",
          "Other and ill-defined digestive organs",
          "Vagina",
          "Stomach",
          "Ovary",
          "Heart, mediastinum, and pleura",
          "Kidney",
          "Anus and anal canal",
          "Nasopharynx",
          "Testis",
          "Brain",
          "Breast",
          "Bronchus and lung",
          "Liver and intrahepatic bile ducts",
          "Other endocrine glands and related structures",
          "Eye and adnexa",
          "Esophagus",
          "Hematopoietic and reticuloendothelial systems",
          "Skin",
          "Other and ill-defined sites",
          "Other and unspecified female genital organs",
          "Spinal cord, cranial nerves, and other parts of central nervous system",
          "Thyroid gland",
          "Other and unspecified major salivary glands",
          "Other and unspecified urinary organs",
          "Cervix uteri",
          "Uterus, NOS",
          "Unknown",
          "Bones, joints and articular cartilage of other and unspecified sites",
          "Retroperitoneum and peritoneum",
          "Connective, subcutaneous and other soft tissues"
        ],
        "dbgap_accession_number": null,
        "project_id": "GENIE-MDA",
        "disease_type": [
          "Mesothelial Neoplasms",
          "Adenomas and Adenocarcinomas",
          "Leukemias, NOS",
          "Blood Vessel Tumors",
          "Fibromatous Neoplasms",
          "Mucoepidermoid Neoplasms",
          "Thymic Epithelial Neoplasms",
          "Trophoblastic neoplasms",
          "Miscellaneous Bone Tumors",
          "Neuroepitheliomatous Neoplasms",
          "Basal Cell Neoplasms",
          "Specialized Gonadal Neoplasms",
          "Complex Epithelial Neoplasms",
          "Lipomatous Neoplasms",
          "Nevi and Melanomas",
          "Ductal and Lobular Neoplasms",
          "Hodgkin Lymphoma",
          "Mature B-Cell Lymphomas",
          "Myomatous Neoplasms",
          "Epithelial Neoplasms, NOS",
          "Fibroepithelial Neoplasms",
          "Gliomas",
          "Neoplasms, NOS",
          "Transitional Cell Papillomas and Carcinomas",
          "Meningiomas",
          "Synovial-like Neoplasms",
          "Germ Cell Neoplasms",
          "Malignant Lymphomas, NOS or Diffuse",
          "Cystic, Mucinous and Serous Neoplasms",
          "Squamous Cell Neoplasms",
          "Osseous and Chondromatous Neoplasms",
          "Miscellaneous Tumors",
          "Nerve Sheath Tumors",
          "Complex Mixed and Stromal Neoplasms"
        ],
        "name": "AACR Project GENIE - Contributed by MD Anderson Cancer Center",
        "releasable": true,
        "state": "open",
        "released": true
      },
      {
        "id": "GENIE-MSK",
        "primary_site": [
          "Other and unspecified parts of biliary tract",
          "Thymus",
          "Gallbladder",
          "Bladder",
          "Colon",
          "Small intestine",
          "Prostate gland",
          "Pancreas",
          "Hypopharynx",
          "Peripheral nerves and autonomic nervous system",
          "Rectum",
          "Adrenal gland",
          "Other and ill-defined sites in lip, oral cavity and pharynx",
          "Corpus uteri",
          "Nasal cavity and middle ear",
          "Meninges",
          "Penis",
          "Other and ill-defined digestive organs",
          "Vagina",
          "Ovary",
          "Heart, mediastinum, and pleura",
          "Stomach",
          "Kidney",
          "Nasopharynx",
          "Anus and anal canal",
          "Testis",
          "Brain",
          "Breast",
          "Bronchus and lung",
          "Liver and intrahepatic bile ducts",
          "Other endocrine glands and related structures",
          "Eye and adnexa",
          "Esophagus",
          "Hematopoietic and reticuloendothelial systems",
          "Skin",
          "Other and ill-defined sites",
          "Other and unspecified female genital organs",
          "Oropharynx",
          "Spinal cord, cranial nerves, and other parts of central nervous system",
          "Thyroid gland",
          "Other and unspecified major salivary glands",
          "Other and unspecified urinary organs",
          "Cervix uteri",
          "Larynx",
          "Uterus, NOS",
          "Unknown",
          "Bones, joints and articular cartilage of other and unspecified sites",
          "Retroperitoneum and peritoneum",
          "Connective, subcutaneous and other soft tissues"
        ],
        "dbgap_accession_number": null,
        "project_id": "GENIE-MSK",
        "disease_type": [
          "Mesothelial Neoplasms",
          "Adenomas and Adenocarcinomas",
          "Blood Vessel Tumors",
          "Fibromatous Neoplasms",
          "Mucoepidermoid Neoplasms",
          "Thymic Epithelial Neoplasms",
          "Trophoblastic neoplasms",
          "Adnexal and Skin Appendage Neoplasms",
          "Miscellaneous Bone Tumors",
          "Immunoproliferative Diseases",
          "Neuroepitheliomatous Neoplasms",
          "Basal Cell Neoplasms",
          "Specialized Gonadal Neoplasms",
          "Mesonephromas",
          "Complex Epithelial Neoplasms",
          "Neoplasms of Histiocytes and Accessory Lymphoid Cells",
          "Myeloid Leukemias",
          "Odontogenic Tumors",
          "Lipomatous Neoplasms",
          "Nevi and Melanomas",
          "Ductal and Lobular Neoplasms",
          "Granular Cell Tumors and Alveolar Soft Part Sarcomas",
          "Acinar Cell Neoplasms",
          "Soft Tissue Tumors and Sarcomas, NOS",
          "Paragangliomas and Glomus Tumors",
          "Hodgkin Lymphoma",
          "Mature B-Cell Lymphomas",
          "Myxomatous Neoplasms",
          "Lymphoid Leukemias",
          "Myomatous Neoplasms",
          "Epithelial Neoplasms, NOS",
          "Fibroepithelial Neoplasms",
          "Gliomas",
          "Plasma Cell Tumors",
          "Neoplasms, NOS",
          "Mature T- and NK-Cell Lymphomas",
          "Mast Cell Tumors",
          "Transitional Cell Papillomas and Carcinomas",
          "Meningiomas",
          "Synovial-like Neoplasms",
          "Germ Cell Neoplasms",
          "Malignant Lymphomas, NOS or Diffuse",
          "Myelodysplastic Syndromes",
          "Squamous Cell Neoplasms",
          "Cystic, Mucinous and Serous Neoplasms",
          "Osseous and Chondromatous Neoplasms",
          "Miscellaneous Tumors",
          "Nerve Sheath Tumors",
          "Complex Mixed and Stromal Neoplasms"
        ],
        "name": "AACR Project GENIE - Contributed by Memorial Sloan Kettering Cancer Center",
        "releasable": true,
        "state": "open",
        "released": true
      }
    ],
    "pagination": {
      "count": 5,
      "total": 78,
      "size": 5,
      "from": 0,
      "sort": "None",
      "page": 1,
      "pages": 16
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
        "id": "0286c31b-a704-4d7d-99e3-0bc4e8975b8b",
        "lost_to_followup": "Yes",
        "days_to_lost_to_followup": null,
        "disease_type": "Ductal and Lobular Neoplasms",
        "analyte_ids": [
          "1e1213a7-4c2b-401f-acf2-38ecd6a084b7",
          "b86795d2-93d6-47ac-b197-ae7c4a3effcd",
          "92b26409-6096-4b23-96e3-cbd63fb22580",
          "684a5a28-dfdc-4554-837f-3fb03d2ed23f",
          "f4a22f4b-5279-42e7-bfd6-f14b23c89268"
        ],
        "submitter_id": "HCM-CSHL-0084-C25",
        "submitter_analyte_ids": [
          "HCM-CSHL-0084-C25-11A-11D",
          "HCM-CSHL-0084-C25-01A-11D",
          "HCM-CSHL-0084-C25-85A-01R",
          "HCM-CSHL-0084-C25-01A-11R",
          "HCM-CSHL-0084-C25-85B-01D"
        ],
        "aliquot_ids": [
          "a686fbd1-cfa8-424b-a4a9-3658515348e3",
          "da425644-5762-455a-9d09-992bb0cf1f7a",
          "c8af1fef-147f-4d2c-b41d-75930368b369",
          "407645cb-481e-4ded-8709-48eebed12d51",
          "2f755ce8-ed34-43b2-8a1c-75e55faef9bf"
        ],
        "submitter_aliquot_ids": [
          "HCM-CSHL-0084-C25-85B-01D-A78M-36",
          "HCM-CSHL-0084-C25-85A-01R-A78N-41",
          "HCM-CSHL-0084-C25-11A-11D-A78M-36",
          "HCM-CSHL-0084-C25-01A-11R-A78N-41",
          "HCM-CSHL-0084-C25-01A-11D-A78M-36"
        ],
        "created_datetime": "2019-09-19T08:58:49.953315-05:00",
        "diagnosis_ids": [
          "70da7d37-590b-41a5-97c7-17238c0c8652"
        ],
        "sample_ids": [
          "b253895e-9f1a-4aef-9353-be454fac8062",
          "06c1a6af-3c39-4454-9667-d3b1b3f68eaa",
          "1483d5e0-9037-4c92-aa5d-c655666af41e",
          "77b0c06f-7407-4af2-9cc0-a9160c105043"
        ],
        "submitter_sample_ids": [
          "HCM-CSHL-0084-C25-85A",
          "HCM-CSHL-0084-C25-01A",
          "HCM-CSHL-0084-C25-85B",
          "HCM-CSHL-0084-C25-11A"
        ],
        "primary_site": "Pancreas",
        "submitter_diagnosis_ids": [
          "HCM-CSHL-0084-C25_diagnosis"
        ],
        "updated_datetime": "2023-02-22T07:39:25.979291-06:00",
        "case_id": "0286c31b-a704-4d7d-99e3-0bc4e8975b8b",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "e4a3579c-2d21-4620-8d5a-2966cd84c4cd",
          "11c9a726-6aab-4259-9d78-dfebf9cd0249"
        ],
        "submitter_portion_ids": [
          "HCM-CSHL-0084-C25-01A-11",
          "HCM-CSHL-0084-C25-11A-11"
        ]
      },
      {
        "id": "02f6d684-b6b5-419a-b0e1-b74d0a384a30",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Gliomas",
        "analyte_ids": [
          "30773bae-0ca1-467c-becc-574e3c9e7b5d",
          "9ef57718-b3f0-48d7-868e-c48452c5464a",
          "6e6b33fb-ce09-4c2c-a372-d468bdc1a11c"
        ],
        "submitter_id": "HCM-BROD-0408-C71",
        "submitter_analyte_ids": [
          "HCM-BROD-0408-C71-85R-01D",
          "HCM-BROD-0408-C71-85R-01R",
          "HCM-BROD-0408-C71-10A-01D"
        ],
        "days_to_consent": null,
        "aliquot_ids": [
          "65c0b96f-c527-489c-bede-4e2c8080c513",
          "7b6f05de-cada-49c2-a7cc-df943fda0fdd",
          "29b97383-b859-4895-8f25-18a943a7319b"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0408-C71-85R-01D-A80U-36",
          "HCM-BROD-0408-C71-85R-01R-A80V-41",
          "HCM-BROD-0408-C71-10A-01D-A80U-36"
        ],
        "created_datetime": "2020-06-15T13:40:11.094740-05:00",
        "diagnosis_ids": [
          "9fc0fac7-b813-4c0f-995d-4045f0db272c"
        ],
        "sample_ids": [
          "d30766c0-22b9-4f63-8990-609b81bc2cfe",
          "29e0bdfe-1618-41dc-90da-03907199151d"
        ],
        "consent_type": null,
        "submitter_sample_ids": [
          "HCM-BROD-0408-C71-85R",
          "HCM-BROD-0408-C71-10A"
        ],
        "primary_site": "Brain",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0408-C71_diagnosis"
        ],
        "updated_datetime": "2023-02-22T07:39:25.979291-06:00",
        "case_id": "02f6d684-b6b5-419a-b0e1-b74d0a384a30",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "0ce874cd-4701-4746-822f-5f60587c6f11"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0408-C71-10A-01"
        ]
      },
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
        "id": "04cbceab-f945-482b-956b-840756a17a4a",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Gliomas",
        "analyte_ids": [
          "6852014a-7cc2-447f-b09f-ba7c53403ab9",
          "45446a95-ef00-4f45-8068-62e733f29129",
          "b256dfa3-f988-4595-a276-d324b83928ac",
          "643c8f31-9fd7-42a8-bcca-aa4528162718",
          "d40b1ebb-c895-4732-9abc-0c9f818037ac"
        ],
        "submitter_id": "HCM-BROD-0421-C71",
        "submitter_analyte_ids": [
          "HCM-BROD-0421-C71-10A-01D",
          "HCM-BROD-0421-C71-01B-11D",
          "HCM-BROD-0421-C71-85A-01R",
          "HCM-BROD-0421-C71-01B-11R",
          "HCM-BROD-0421-C71-85A-01D"
        ],
        "days_to_consent": null,
        "aliquot_ids": [
          "1689e951-b01b-4f81-8c00-83e3904058d0",
          "e3747c4d-7822-494c-9ed2-30a632ffa8a6",
          "ae4cc7ec-72ad-4d6c-b061-f67615d4ad2d",
          "2cf507a2-a338-4b22-a176-bdb91dde35de",
          "44270718-c802-4fb9-b784-fc8bc8d7a1aa"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0421-C71-01B-11D-A80U-36",
          "HCM-BROD-0421-C71-01B-11R-A80V-41",
          "HCM-BROD-0421-C71-10A-01D-A80U-36",
          "HCM-BROD-0421-C71-85A-01R-A80V-41",
          "HCM-BROD-0421-C71-85A-01D-A80U-36"
        ],
        "created_datetime": "2020-07-01T12:59:56.694185-05:00",
        "diagnosis_ids": [
          "252e0d9f-6cdb-4849-b981-e5a5ef0bb151"
        ],
        "sample_ids": [
          "009edc9d-3acf-4e70-b49b-5e83e54ae258",
          "4c8febaa-bbc4-4870-8960-328da9cac009",
          "6039b795-fe2d-4a71-b9b9-2e6caedb5cb5"
        ],
        "consent_type": null,
        "submitter_sample_ids": [
          "HCM-BROD-0421-C71-10A",
          "HCM-BROD-0421-C71-85A",
          "HCM-BROD-0421-C71-01B"
        ],
        "primary_site": "Brain",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0421-C71_diagnosis"
        ],
        "updated_datetime": "2023-02-22T07:39:25.979291-06:00",
        "case_id": "04cbceab-f945-482b-956b-840756a17a4a",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "20322504-bde0-4b32-8fda-f5d7dc5bcc47",
          "87467751-20cc-4889-8d94-4adda1d74d76"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0421-C71-01B-11",
          "HCM-BROD-0421-C71-10A-01"
        ]
      }
    ],
    "pagination": {
      "count": 5,
      "total": 86962,
      "size": 5,
      "from": 0,
      "sort": "",
      "page": 1,
      "pages": 17393
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
        "id": "99a26957-3e1b-4c9f-835f-96859c4475ef",
        "entity_submitter_id": "2432eaeb-d928-4e02-ab48-8349da12e909",
        "notes": "Real somatic mutations were mistakenly labeled as LOH (Loss of Heterozygosity) in certain SomaticSniper VCF files.",
        "classification": "Notification",
        "entity_id": "b191c772-48b8-4b41-a16d-432bb36517a1",
        "created_datetime": "2022-02-27T22:14:47.165866-06:00",
        "annotation_id": "99a26957-3e1b-4c9f-835f-96859c4475ef",
        "entity_type": "annotated_somatic_mutation",
        "updated_datetime": "2022-02-27T22:14:47.165866-06:00",
        "case_id": "03bfeb7c-cecf-4691-8263-33cdfe391ea9",
        "state": "released",
        "category": "General",
        "status": "Approved",
        "case_submitter_id": "HCM-BROD-0124-C25"
      },
      {
        "id": "581cd77a-632b-4ce0-ae4f-287814271a56",
        "entity_submitter_id": "ce49859e-afb5-4f8d-b687-e503c2b7ddc8",
        "notes": "Real somatic mutations were mistakenly labeled as LOH (Loss of Heterozygosity) in certain SomaticSniper VCF files.",
        "classification": "Notification",
        "entity_id": "9ba2dd7b-7f5d-4385-abc0-fc5af745ff05",
        "created_datetime": "2022-02-28T05:10:55.897983-06:00",
        "annotation_id": "581cd77a-632b-4ce0-ae4f-287814271a56",
        "entity_type": "annotated_somatic_mutation",
        "updated_datetime": "2023-02-22T07:39:25.979291-06:00",
        "case_id": "78b05ced-6b5e-44e0-b8e2-07ab9b69e440",
        "state": "released",
        "category": "General",
        "status": "Approved",
        "case_submitter_id": "HCM-BROD-0716-C43"
      }
    ],
    "pagination": {
      "count": 2,
      "total": 75111,
      "size": 2,
      "from": 0,
      "sort": "",
      "page": 1,
      "pages": 37556
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
          "Skin"
        ]
      }
    }
  ]
}
```
```Query
curl 'https://api.gdc.cancer.gov/cases?filters=%7B%0A%20%20%22op%22%3A%20%22and%22%2C%0A%20%20%22content%22%3A%20%5B%0A%20%20%20%20%7B%0A%20%20%20%20%20%20%22op%22%3A%20%22in%22%2C%0A%20%20%20%20%20%20%22content%22%3A%20%7B%0A%20%20%20%20%20%20%20%20%22field%22%3A%20%22primary_site%22%2C%0A%20%20%20%20%20%20%20%20%22value%22%3A%20%5B%0A%20%20%20%20%20%20%20%20%20%20%22Skin%22%0A%20%20%20%20%20%20%20%20%5D%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%7D%0A%20%20%5D%0A%7D&pretty=true'
```
```Response
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
        "id": "1c54d889-a752-4025-90ce-af896cdb95d4",
        "lost_to_followup": null,
        "slide_ids": [
          "21864471-d2bf-46d7-aa99-de91e2bd975f",
          "bb901d03-1075-4ad9-88ce-61f2fde2505d"
        ],
        "submitter_slide_ids": [
          "HCM-BROD-0643-C43-06A-01-S2-HE",
          "HCM-BROD-0643-C43-06A-01-S1-HE"
        ],
        "days_to_lost_to_followup": null,
        "disease_type": "Nevi and Melanomas",
        "analyte_ids": [
          "8c8b84e5-af2e-497e-beb6-addea8e64a4b",
          "9dda998f-416d-4878-88ee-ab2e6aa800eb",
          "4d8d9d56-957d-4d65-a3ba-73b266977f2e",
          "b818d10c-b6c0-499c-8997-438d646a3b6f",
          "d35a1c25-f02e-4cde-a48c-050528cf1923"
        ],
        "submitter_id": "HCM-BROD-0643-C43",
        "submitter_analyte_ids": [
          "HCM-BROD-0643-C43-85N-01R",
          "HCM-BROD-0643-C43-85N-01D",
          "HCM-BROD-0643-C43-06A-11D",
          "HCM-BROD-0643-C43-10A-01D",
          "HCM-BROD-0643-C43-06A-11R"
        ],
        "days_to_consent": null,
        "aliquot_ids": [
          "eda70710-057b-47dc-8a3b-b070b51e5cf7",
          "5e664a16-9ed4-4d47-876c-560bbd460211",
          "62565e48-482f-43a3-8345-f126a906fcf9",
          "d1567468-ec2c-4aae-ab8e-7487023fdae0",
          "0986f594-ef25-4f51-b802-ddfacd38ee34"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0643-C43-06A-11D-A80U-36",
          "HCM-BROD-0643-C43-10A-01D-A80U-36",
          "HCM-BROD-0643-C43-06A-11R-A80V-41",
          "HCM-BROD-0643-C43-85N-01D-A80U-36",
          "HCM-BROD-0643-C43-85N-01R-A80V-41"
        ],
        "created_datetime": "2020-07-01T13:17:19.998786-05:00",
        "diagnosis_ids": [
          "3ea0fe35-fd02-4e73-b294-965afe66a719",
          "379deec7-50a8-4219-9bd2-1737e7908d50"
        ],
        "sample_ids": [
          "b68abcbc-6310-4422-a2ae-8134f2fb9bf5",
          "6dbb5557-ba67-49d4-b4fa-b333246788f7",
          "90aa1d3f-7ef3-4561-ad1e-9153215e756f"
        ],
        "consent_type": null,
        "submitter_sample_ids": [
          "HCM-BROD-0643-C43-06A",
          "HCM-BROD-0643-C43-10A",
          "HCM-BROD-0643-C43-85N"
        ],
        "primary_site": "Skin",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0643-C43_diagnosis2",
          "HCM-BROD-0643-C43_diagnosis"
        ],
        "updated_datetime": "2021-03-11T14:33:32.013919-06:00",
        "case_id": "1c54d889-a752-4025-90ce-af896cdb95d4",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "53aa58c5-f1f1-47a7-b2f4-6ee81b772568",
          "f6aebf1d-448b-4fd2-bd2f-08ac8de2c666"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0643-C43-06A-11",
          "HCM-BROD-0643-C43-10A-01"
        ]
      },
      {
        "id": "1c906cee-df63-4348-9354-0a523ab88dfa",
        "disease_type": "Nevi and Melanomas",
        "analyte_ids": [
          "4b6c6d11-172a-42a6-8a30-59a2055ed2d2",
          "e340feee-6025-4f73-9f23-fade5a4169b7",
          "f12b690e-bd8d-47b4-bf51-a1b488730e55",
          "3a6ed58f-61c8-4516-8dc9-b27f68878e93",
          "415de04e-7c98-47d4-8a98-896d718c2831"
        ],
        "submitter_id": "HCM-BROD-0758-C43",
        "submitter_analyte_ids": [
          "HCM-BROD-0758-C43-06A-01D",
          "HCM-BROD-0758-C43-10A-01D",
          "HCM-BROD-0758-C43-85N-01R",
          "HCM-BROD-0758-C43-85N-01D",
          "HCM-BROD-0758-C43-06A-01R"
        ],
        "aliquot_ids": [
          "14ec708e-68de-4613-86fa-3214337f37af",
          "ca35c0b9-62a7-4f04-b930-918d4b4867eb",
          "e9b9a3e8-90f8-40b3-a0f9-7c07bfe30084",
          "bb21ec5d-7bc0-4329-9315-7ce37379e7f0",
          "be10bb2b-f25b-4636-8402-7a15b71e6e7b"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0758-C43-06A-01R-A85D-41",
          "HCM-BROD-0758-C43-85N-01R-A85D-41",
          "HCM-BROD-0758-C43-10A-01D-A85C-36",
          "HCM-BROD-0758-C43-85N-01D-A85C-36",
          "HCM-BROD-0758-C43-06A-01D-A85C-36"
        ],
        "created_datetime": "2021-06-03T11:55:23.082303-05:00",
        "diagnosis_ids": [
          "640cca22-1494-46df-82e5-8f9bffd253ce",
          "8a52bb8c-67a6-4130-b4e6-19e0d0ed4c52"
        ],
        "sample_ids": [
          "a61d1220-138c-482e-be8d-91f13a394d34",
          "933dc8cb-176f-4df8-8b72-b492e9a4fb9e",
          "6d16bdfb-35ee-4efa-a7ee-e95c6ad861af"
        ],
        "submitter_sample_ids": [
          "HCM-BROD-0758-C43-85N",
          "HCM-BROD-0758-C43-06A",
          "HCM-BROD-0758-C43-10A"
        ],
        "primary_site": "Skin",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0758-C43_diagnosis",
          "HCM-BROD-0758-C43_diagnosis2"
        ],
        "updated_datetime": "2023-02-22T07:39:25.979291-06:00",
        "case_id": "1c906cee-df63-4348-9354-0a523ab88dfa",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "6c5bf5f0-4dac-4397-89f9-4a06ce49883d",
          "da103752-c408-46b4-aed4-69d51bf80170"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0758-C43-06A-01",
          "HCM-BROD-0758-C43-10A-01"
        ]
      },
      {
        "id": "3b84c894-91d7-4c3e-8f17-4d6adc8bb468",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Nevi and Melanomas",
        "analyte_ids": [
          "38dbed8f-f9c9-4f98-a76f-b44031a1e83c",
          "6d410531-100f-46f5-8ad0-f05b0dc6eba5",
          "f14b1ef0-f2a8-4c03-a451-4b4b8cc076c2",
          "34550b5b-3cfa-49cc-92ae-d9eebfc6dac7"
        ],
        "submitter_id": "HCM-BROD-0339-C43",
        "submitter_analyte_ids": [
          "HCM-BROD-0339-C43-10A-01D",
          "HCM-BROD-0339-C43-06A-01D",
          "HCM-BROD-0339-C43-85M-01D",
          "HCM-BROD-0339-C43-85M-01R"
        ],
        "days_to_consent": null,
        "aliquot_ids": [
          "8b878538-7036-4f10-a4a6-619151156020",
          "cb9b9bbe-34fb-451a-923d-4b9f6103b3cd",
          "8adc893c-8bb2-4fdb-b716-40845a595e21",
          "0f698098-853a-481e-82a4-d9be5793b562"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0339-C43-10A-01D-A79C-36",
          "HCM-BROD-0339-C43-85M-01D-A79C-36",
          "HCM-BROD-0339-C43-85M-01R-A79D-41",
          "HCM-BROD-0339-C43-06A-01D-A79C-36"
        ],
        "created_datetime": "2020-01-08T16:03:51.205115-06:00",
        "diagnosis_ids": [
          "f2100f46-97ee-4a53-82a8-95c975e1b147",
          "1c9f6d6d-12e9-4e4d-b301-8e9ec59c12d8"
        ],
        "sample_ids": [
          "816e34cd-335c-48f4-b2d0-3eb6c4bbe39f",
          "b129099e-1ad5-4b5d-a47d-2dbf1f44a5c5",
          "b37c5fa5-2dc5-462c-a957-b09ebd55e2d8"
        ],
        "consent_type": null,
        "submitter_sample_ids": [
          "HCM-BROD-0339-C43-06A",
          "HCM-BROD-0339-C43-85M",
          "HCM-BROD-0339-C43-10A"
        ],
        "primary_site": "Skin",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0339-C43_diagnosis2",
          "HCM-BROD-0339-C43_diagnosis"
        ],
        "updated_datetime": "2021-01-06T22:55:10.531130-06:00",
        "case_id": "3b84c894-91d7-4c3e-8f17-4d6adc8bb468",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "3829676e-9cdf-4a2b-bb42-e29536e760ed",
          "255a1955-5d1f-4764-aba1-f3b5334424d3"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0339-C43-06A-01",
          "HCM-BROD-0339-C43-10A-01"
        ]
      },
      {
        "id": "47322ea3-6bbe-442b-a656-c48469cc99c1",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Nevi and Melanomas",
        "analyte_ids": [
          "f9a2e102-f495-4ec8-b038-25beaa892cc1",
          "107931a1-6be8-4b2c-8f7c-de22523b3045",
          "e7fd0047-27bb-4ffd-837f-6ff05d54baaf"
        ],
        "submitter_id": "HCM-BROD-0335-C43",
        "submitter_analyte_ids": [
          "HCM-BROD-0335-C43-10A-01D",
          "HCM-BROD-0335-C43-85M-01D",
          "HCM-BROD-0335-C43-85M-01R"
        ],
        "aliquot_ids": [
          "839a15d8-9c22-4f82-a89a-c6dcf8c1e0bd",
          "06d1a983-858f-454a-9b0d-448bf29670cc",
          "31c27aec-0049-4da6-b748-0dd5d93d76c5"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0335-C43-85M-01R-A78V-41",
          "HCM-BROD-0335-C43-10A-01D-A78T-36",
          "HCM-BROD-0335-C43-85M-01D-A78T-36"
        ],
        "created_datetime": "2019-10-14T10:42:44.170730-05:00",
        "diagnosis_ids": [
          "280cf92a-e646-417e-bb6b-8d34b6d13f43",
          "207053c5-a053-449c-94c4-a3e36f19b161"
        ],
        "sample_ids": [
          "2ccdf84e-0d53-4b91-9dcc-939951f722b6",
          "6cacdfa9-89e3-4984-926c-69e94ee7ca1a"
        ],
        "submitter_sample_ids": [
          "HCM-BROD-0335-C43-85M",
          "HCM-BROD-0335-C43-10A"
        ],
        "primary_site": "Skin",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0335-C43_diagnosis2",
          "HCM-BROD-0335-C43_diagnosis"
        ],
        "updated_datetime": "2021-07-12T12:25:55.528644-05:00",
        "case_id": "47322ea3-6bbe-442b-a656-c48469cc99c1",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "b848a38c-2e03-4518-b5e3-0f6ab0144719"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0335-C43-10A-01"
        ]
      },
      {
        "id": "4e828c35-5fc3-40af-a2ee-10b4281c77b8",
        "lost_to_followup": null,
        "slide_ids": [
          "b4917938-a3f4-4fb5-a0a8-f8d4fbdf0d3d",
          "e3c181f3-a364-43a4-946a-d4882d7ef68e"
        ],
        "submitter_slide_ids": [
          "HCM-BROD-0569-C43-06A-01-S2-HE",
          "HCM-BROD-0569-C43-06A-01-S1-HE"
        ],
        "days_to_lost_to_followup": null,
        "disease_type": "Nevi and Melanomas",
        "analyte_ids": [
          "e516dc5f-729f-46b3-8a44-a0db98f47eb8",
          "fee1b1e1-6604-4da1-950e-991f68f3dc05",
          "3d51694f-9814-483a-b0d0-25e81b453624",
          "b2b8b533-a8fb-4521-a09a-483cc833c08f",
          "a829df29-620c-406b-ac8c-ca8565ce9aa9"
        ],
        "submitter_id": "HCM-BROD-0569-C43",
        "submitter_analyte_ids": [
          "HCM-BROD-0569-C43-10A-01D",
          "HCM-BROD-0569-C43-85M-01R",
          "HCM-BROD-0569-C43-06A-11D",
          "HCM-BROD-0569-C43-85M-01D",
          "HCM-BROD-0569-C43-06A-11R"
        ],
        "days_to_consent": null,
        "aliquot_ids": [
          "68b0d1e1-1da8-4331-9178-89a736b07fff",
          "56bb2b64-77e5-4155-8a62-328eb67036e3",
          "e44bf205-ac59-4a29-bf50-1b9da39671a4",
          "8e4a9a37-f5b3-4720-86e9-8999126c2a6f",
          "8bd0ad30-1df5-49bf-847b-0aecbf846043"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0569-C43-06A-11D-A80U-36",
          "HCM-BROD-0569-C43-85M-01D-A80U-36",
          "HCM-BROD-0569-C43-06A-11R-A80V-41",
          "HCM-BROD-0569-C43-10A-01D-A80U-36",
          "HCM-BROD-0569-C43-85M-01R-A80V-41"
        ],
        "created_datetime": "2020-07-01T13:12:41.873804-05:00",
        "diagnosis_ids": [
          "5191a8e6-4ec6-4a84-af18-70d24c53d97d",
          "e795772c-a156-49cc-8578-4ba5901f374e"
        ],
        "sample_ids": [
          "dafc8b78-87eb-4210-bb30-6356aa644f8a",
          "68e1cc96-9b30-4724-9842-cc3ff548d569",
          "ad9ab2ba-a70b-4c57-9aa6-8965d0be6071"
        ],
        "consent_type": null,
        "submitter_sample_ids": [
          "HCM-BROD-0569-C43-10A",
          "HCM-BROD-0569-C43-06A",
          "HCM-BROD-0569-C43-85M"
        ],
        "primary_site": "Skin",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0569-C43_diagnosis2",
          "HCM-BROD-0569-C43_diagnosis"
        ],
        "updated_datetime": "2021-03-11T14:31:21.933263-06:00",
        "case_id": "4e828c35-5fc3-40af-a2ee-10b4281c77b8",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "da022483-8c21-4fbd-9c0b-2a64dea62a55",
          "a2285ada-e582-466f-b977-d0e00921ce41"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0569-C43-06A-11",
          "HCM-BROD-0569-C43-10A-01"
        ]
      },
      {
        "id": "6ded7d89-1235-4b52-90aa-52c57c669301",
        "lost_to_followup": null,
        "slide_ids": [
          "a47ad461-0c22-45f6-b1bd-b133d74a2efe",
          "e6feb344-ef68-4396-a492-0bc3a45fda12"
        ],
        "submitter_slide_ids": [
          "HCM-BROD-0557-C43-06A-01-S1-HE",
          "HCM-BROD-0557-C43-06A-01-S2-HE"
        ],
        "days_to_lost_to_followup": null,
        "disease_type": "Nevi and Melanomas",
        "analyte_ids": [
          "88a810b1-530b-4e27-b702-c47fb1687e77",
          "9fdad739-b2aa-443c-8c1c-3837864d90fa",
          "43e7a638-003e-4186-93c2-7804a28e1d52",
          "f30e086a-453c-4d9c-a81f-925c35300534",
          "b2870252-b41e-4053-ab3c-e4fd14299f09"
        ],
        "submitter_id": "HCM-BROD-0557-C43",
        "submitter_analyte_ids": [
          "HCM-BROD-0557-C43-06A-11R",
          "HCM-BROD-0557-C43-85M-01D",
          "HCM-BROD-0557-C43-06A-11D",
          "HCM-BROD-0557-C43-10A-01D",
          "HCM-BROD-0557-C43-85M-01R"
        ],
        "days_to_consent": null,
        "aliquot_ids": [
          "51c2fa6b-eefa-4510-9a9c-adcc06865865",
          "75fda0d3-eb62-4ddf-bc41-3ade57851870",
          "a8eb1adc-2ae7-423f-9c86-2a4ac6cf7c23",
          "b95fb706-4a05-41cb-bfd0-fa46afaef8d9",
          "da3182a4-2327-495a-88ec-2069757dab87"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0557-C43-85M-01R-A80V-41",
          "HCM-BROD-0557-C43-06A-11D-A80U-36",
          "HCM-BROD-0557-C43-06A-11R-A80V-41",
          "HCM-BROD-0557-C43-10A-01D-A80U-36",
          "HCM-BROD-0557-C43-85M-01D-A80U-36"
        ],
        "created_datetime": "2020-07-01T13:20:18.566669-05:00",
        "diagnosis_ids": [
          "4e02019a-674d-49a1-9d7c-af3daeaec684",
          "cce3cf4e-5bc1-475b-908e-113edf8aca5d"
        ],
        "sample_ids": [
          "61123977-af2f-4d0c-8fb4-39e35c76569a",
          "fd17be83-b1ce-4c36-867b-d3fb0c6c7142",
          "e22af195-1cda-484f-ab99-ad508dd8bd8d"
        ],
        "consent_type": null,
        "submitter_sample_ids": [
          "HCM-BROD-0557-C43-85M",
          "HCM-BROD-0557-C43-10A",
          "HCM-BROD-0557-C43-06A"
        ],
        "primary_site": "Skin",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0557-C43_diagnosis",
          "HCM-BROD-0557-C43_diagnosis2"
        ],
        "updated_datetime": "2021-03-11T14:30:11.765496-06:00",
        "case_id": "6ded7d89-1235-4b52-90aa-52c57c669301",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "34a0076b-12b5-4fbf-8568-05177e163ff0",
          "903a130b-bea0-4e3f-9a7b-31ba1aa3ca21"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0557-C43-10A-01",
          "HCM-BROD-0557-C43-06A-11"
        ]
      },
      {
        "id": "71dd56c3-0048-4ec9-89cb-53b5e966255a",
        "lost_to_followup": "Yes",
        "disease_type": "Nevi and Melanomas",
        "analyte_ids": [
          "7ede93ee-b336-4985-9a3f-a3ae69489e74",
          "1271a757-a30d-409d-bff8-5f8373dba1c7",
          "c9ab7c21-69cd-40a0-be3c-b1387db18b41"
        ],
        "submitter_id": "HCM-BROD-0702-C43",
        "submitter_analyte_ids": [
          "HCM-BROD-0702-C43-85M-01R",
          "HCM-BROD-0702-C43-85M-01D",
          "HCM-BROD-0702-C43-10A-01D"
        ],
        "aliquot_ids": [
          "724da8ab-341f-45f0-afe0-e1f41acb7def",
          "2c853383-6f59-48e4-a746-e7ceda4614d2",
          "b7cecd5a-7213-45c2-98bd-9efe487220a5"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0702-C43-10A-01D-A82H-36",
          "HCM-BROD-0702-C43-85M-01D-A82H-36",
          "HCM-BROD-0702-C43-85M-01R-A82J-41"
        ],
        "created_datetime": "2020-10-30T18:06:47.894645-05:00",
        "diagnosis_ids": [
          "9ce4d590-fc03-4b3f-9cb3-012ae0639d37",
          "228d20c0-aa60-46f5-938b-3b10fcdea7dd"
        ],
        "sample_ids": [
          "811ac38b-49ce-4151-8d2a-b18d63976993",
          "22d64afc-53fe-4f15-a785-012efd7a75e3"
        ],
        "submitter_sample_ids": [
          "HCM-BROD-0702-C43-10A",
          "HCM-BROD-0702-C43-85M"
        ],
        "primary_site": "Skin",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0702-C43_diagnosis2",
          "HCM-BROD-0702-C43_diagnosis"
        ],
        "updated_datetime": "2023-02-22T07:39:25.979291-06:00",
        "case_id": "71dd56c3-0048-4ec9-89cb-53b5e966255a",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "413604a8-a20d-4da4-9051-083ecc44754f"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0702-C43-10A-01"
        ]
      },
      {
        "id": "78b05ced-6b5e-44e0-b8e2-07ab9b69e440",
        "disease_type": "Nevi and Melanomas",
        "analyte_ids": [
          "ae460685-2852-4bc7-87f9-b324c76d0558",
          "d57bf160-5c4b-4282-bd91-2c8a23356c23",
          "64309926-7ff2-4de6-b07d-0abbd6eee3c2",
          "63966681-3f4f-404f-9211-eae9b235ebcd",
          "f928c0d6-8f0a-40cc-8d67-59555b0b1e1e"
        ],
        "submitter_id": "HCM-BROD-0716-C43",
        "submitter_analyte_ids": [
          "HCM-BROD-0716-C43-06B-01R",
          "HCM-BROD-0716-C43-85N-01D",
          "HCM-BROD-0716-C43-06B-01D",
          "HCM-BROD-0716-C43-10A-01D",
          "HCM-BROD-0716-C43-85N-01R"
        ],
        "aliquot_ids": [
          "4acb7fa8-42b7-43c9-bf49-65871fe91266",
          "87afffdd-0116-4577-944d-826f0539cf96",
          "6410ddef-a4aa-4e43-8206-2f6abf8b0cc9",
          "f2e67431-ac19-4ddb-874a-541cb159a49e",
          "c34938b3-41d9-4ddd-83c0-6e4e44a2b823"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0716-C43-06B-01R-A83V-41",
          "HCM-BROD-0716-C43-85N-01R-A83V-41",
          "HCM-BROD-0716-C43-06B-01D-A83U-36",
          "HCM-BROD-0716-C43-85N-01D-A83U-36",
          "HCM-BROD-0716-C43-10A-01D-A83U-36"
        ],
        "created_datetime": "2021-02-18T16:46:32.260492-06:00",
        "diagnosis_ids": [
          "05a718a8-d606-4fd7-aac3-0a2ff8f0960a",
          "3c60679c-040d-4809-aa03-24d976c7d4c2"
        ],
        "sample_ids": [
          "71fe789c-8b80-4bfd-b8c4-00868d725e46",
          "784ff8ad-8f7f-42fd-be62-bcbd3e1e20b1",
          "af10a9ce-db9d-48f7-b72f-022ca11f7654"
        ],
        "submitter_sample_ids": [
          "HCM-BROD-0716-C43-85N",
          "HCM-BROD-0716-C43-10A",
          "HCM-BROD-0716-C43-06B"
        ],
        "primary_site": "Skin",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0716-C43_diagnosis2",
          "HCM-BROD-0716-C43_diagnosis"
        ],
        "updated_datetime": "2023-02-22T07:39:25.979291-06:00",
        "case_id": "78b05ced-6b5e-44e0-b8e2-07ab9b69e440",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "7d850d7a-83a0-4436-9c96-658f98505bab",
          "754e0275-3ee9-42e2-8b11-3e0e3cf932d7"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0716-C43-10A-01",
          "HCM-BROD-0716-C43-06B-01"
        ]
      },
      {
        "id": "7ec96446-8a5e-400a-966b-986452c972e4",
        "lost_to_followup": null,
        "days_to_lost_to_followup": null,
        "disease_type": "Nevi and Melanomas",
        "analyte_ids": [
          "5254b24f-0249-4e25-9976-ddb94bbb1ee2",
          "936fe6c5-b03f-4523-bd25-640021351611",
          "66e6a1ba-cd18-4080-8706-98b218051c66"
        ],
        "submitter_id": "HCM-BROD-0227-C43",
        "submitter_analyte_ids": [
          "HCM-BROD-0227-C43-10A-01D",
          "HCM-BROD-0227-C43-85M-01D",
          "HCM-BROD-0227-C43-85M-01R"
        ],
        "aliquot_ids": [
          "2f1e377a-74d2-425d-841e-dec6329c9cc3",
          "38b885f6-7017-4910-bf17-6e233f084122",
          "2a96d2b9-f90a-4ef5-ab05-5554bb148a9a"
        ],
        "submitter_aliquot_ids": [
          "HCM-BROD-0227-C43-85M-01D-A78T-36",
          "HCM-BROD-0227-C43-85M-01R-A78V-41",
          "HCM-BROD-0227-C43-10A-01D-A78T-36"
        ],
        "created_datetime": "2019-10-14T10:42:12.652588-05:00",
        "diagnosis_ids": [
          "538d631d-9d1b-4fac-b5df-36547334d590",
          "0ce3d300-8b58-4f01-98ee-d4c4ff4b664b"
        ],
        "sample_ids": [
          "b0e2c28b-6ac9-41ed-85fc-0a91f0049193",
          "f7b79a6c-610e-450b-98a9-04d1095d40ad"
        ],
        "submitter_sample_ids": [
          "HCM-BROD-0227-C43-85M",
          "HCM-BROD-0227-C43-10A"
        ],
        "primary_site": "Skin",
        "submitter_diagnosis_ids": [
          "HCM-BROD-0227-C43_diagnosis",
          "HCM-BROD-0227-C43_diagnosis2"
        ],
        "updated_datetime": "2021-01-06T22:55:10.531130-06:00",
        "case_id": "7ec96446-8a5e-400a-966b-986452c972e4",
        "index_date": "Diagnosis",
        "state": "released",
        "portion_ids": [
          "0ed2f9c1-da3e-4cd2-ae8b-c8feea4a1e7a"
        ],
        "submitter_portion_ids": [
          "HCM-BROD-0227-C43-10A-01"
        ]
      }
    ],
    "pagination": {
      "count": 10,
      "total": 2898,
      "size": 10,
      "from": 0,
      "sort": "",
      "page": 1,
      "pages": 290
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
