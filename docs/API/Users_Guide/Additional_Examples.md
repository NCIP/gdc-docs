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
