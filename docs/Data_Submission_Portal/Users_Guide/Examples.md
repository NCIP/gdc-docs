# Examples
This guide will detail step-by-step procedures for data submission and how they relate to the GDC Data Model and structure.    

## Submitting a BAM Alignment: Data Model Basics

Pictured below is the submittable subset of the GDC Data Model.  This will work as a roadmap the user can follow for submission. The nodes that make up the completed portion of the submission process are highlighted in __green__ and the current target node for submission is highlighted in __red__. Because BAM files are aligned, they fall into the "Submitted Aligned Reads" category of the GDC.  Before submission can begin, the Program and Project must be approved and set by the GDC, which is why they are highlighted from the start.

[![GDC Data Model 1](images/DataModel-1.jpg)](images/DataModel-1.jpg "Click to see the full image.")

### Case Submission

The next step in the Data Model is the `case`, which must be registered beforehand with the GDC. The first step to submitting a `case` is to consult the [Data Dictionary](https://gdc-docs.nci.nih.gov/Data_Dictionary/viewer/#data-dictionary-viewer), which details the fields that are associated with a `case`, the fields that are required to submit a `case`, and the values that can populate each field. Although it is not explicitly stated in the Dictionary entry, each and every entity requires a connection to another entity and a `submitter_id`.

[![Dictionary Case](images/Dictionary_Case.png)](images/Dictionary_Case.png "Click to see the full image.")

Submitting a [__Case__](https://gdc-docs.nci.nih.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=case) etity requires:

* __`submitter_id`:__ A unique key to identify the `case`
* __`projects.code`:__ A link to the `project`

__Note:__ The submitter ID is different from the UUID, which is assigned to each entity automatically.

The `case` entity can be added in JSON or TSV format. A template for any entity in either of these formats can be found in the Data Dictionary at the top of each page.  

```JSON
{
    "type": "case",
    "submitter_id": "GDC-INTERNAL-000055",
    "projects": {
        "code": "INTERNAL"
    }
}
```
```TSV
type  submitter_id  projects.code
case  GDC-INTERNAL-000055 INTERNAL   
```

__Note:__ JSON and TSV formats handle links between entities (`case` and `project`) differently.  JSON includes the `code` field nested within `projects` while TSV appends `code` with a period.  

[![GDC Data Model 2](images/DataModel-2.jpg)](images/DataModel-2.jpg "Click to see the full image.")

### Uploading the Case Submission

The file detailed above can be uploaded in one of two ways:

__Data Submission Portal:__ Each file can be uploaded using the Data Submission Portal Upload Data Wizard. Log into the Submission Portal, choose a project, and choose UPLOAD on the [Dashboard](Dashboard.md).  The Wizard will then pop up, allowing the user to drag-and-drop each file into the window, or browse to retrieve the file.  Next, choose VALIDATE and the system will determine if the file is valid for submission.  When the file has been validated, two buttons will appear:  COMMIT and DISCARD.  Choose COMMIT to upload the file to its project or DISCARD to remove the file. See the [Data Wizard Upload Guide](Upload_Data/#step-2-upload-data-wizard) for more details.   

__API:__ The API has a much broader range of functionality than the Data Wizard. Entities can be created, updated, and deleted through the API. Generally uploading a file through the API will look like:

```Shell
curl --header "X-Auth-Token: $token" --request POST --data @CASE.json https://gdc-api.nci.nih.gov/v0/submission/GDC/INTERNAL/_dry_run?async=true
```
__CASE.json:__ the JSON-formatted `case` file above.  

Next, the file can either be committed through the Data Submission Portal like before, or another API query can be performed that will upload the file to the project.

```Shell
curl --header "X-Auth-Token: $token" --request POST https://gdc-api.nci.nih.gov/v0/submission/GDC/INTERNAL/transactions/467/commit?async=true
```

See the [API Submission User Guide](API/Users_Guide/Submission/#creating-and-updating-entities) for details on submission and the rest of the functionalities of the API.



### Sample Submission

A `sample` submission has the same general structure as `case` submission as it will require a unique key and a link to the `case`.  However, `sample` entities require one additional value:  `sample_type`. This peripheral data is required because it is necessary for any data to be interpreted. For example, we would need to know whether the `sample` came from tumor or normal tissue.  


[![Dictionary Sample](images/Dictionary_Sample.png)](images/Dictionary_Sample.png "Click to see the full image.")

Submitting a [__Sample__](https://gdc-docs.nci.nih.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=sample) entity requires:

* __`submitter_id`:__ A unique key to identify the `sample`
* __`cases.submitter_id`:__ The unique key that was used for the `case`, links the `sample` to the `case`
* __`sample_type`:__ The source or cellular type of the `sample`

```JSON
{
    "type": "sample",
    "cases": {
        "submitter_id": "GDC-INTERNAL-000055"
    },
    "sample_type": "Blood Derived Normal",
    "submitter_id": "Blood-00001_api55"
}
```
```TSV
type	cases.submitter_id	submitter_id	sample_type
sample	GDC-INTERNAL-000055	Blood-00001_api55	Blood Derived Normal  
```

[![GDC Data Model 3](images/DataModel-3.jpg)](images/DataModel-3.jpg "Click to see the full image.")

### Portion, Analyte, and Aliquot Submission

Submitting a [__Portion__](https://gdc-docs.nci.nih.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=portion) entity requires:

* __`submitter_id`:__ A unique key to identify the `portion`
* __`samples.submitter_id`:__ The unique key that was used for the `sample`, links the `portion` to the `sample`

```JSON
{
    "type": "portion",
    "submitter_id": "Blood-portion-000055",
    "samples": {
        "submitter_id": "Blood-00001_api55"
    }
}

```
```TSV
type	submitter_id	samples.submitter_id
portion	Blood-portion-000055	Blood-00001_api55
```

Submitting an [__Analyte__](https://gdc-docs.nci.nih.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=analyte) entity requires:

* __`submitter_id`:__ A unique key to identify the `analyte`
* __`portions.submitter_id`:__ The unique key that was used for the `portion`, links the `analyte` to the `portion`
* __`analyte_type`:__ The protocol-specific molecular type of the `analyte`

```JSON
{
    "type": "analyte",
    "portions": {
        "submitter_id": "Blood-portion-000055"
    },
    "analyte_type": "DNA",
    "submitter_id": "Blood-analyte-000055"
}

```
```TSV
type	portions.submitter_id	analyte_type	submitter_id
analyte	Blood-portion-000055	DNA	Blood-analyte-000055
```

Submitting an [__Aliquot__](https://gdc-docs.nci.nih.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=aliquot) entity requires:

* __`submitter_id`:__ A unique key to identify the `aliquot`
* __`analytes.submitter_id`:__ The unique key that was used for the `analyte`, links the `aliquot` to the `analyte`

```JSON
{
    "type": "aliquot",
    "submitter_id": "Blood-00021-aliquot55",
    "analytes": {
        "submitter_id": "Blood-analyte-000055"
    }
}

```
```TSV
type	submitter_id	analytes.submitter_id
aliquot	Blood-00021-aliquot55	Blood-analyte-000055
```

__Note:__ Aliquot entities can be directly linked to Sample entities.

[![GDC Data Model 4](images/DataModel-4.jpg)](images/DataModel-4.jpg "Click to see the full image.")

### Read Group Submission
Because information about genomic reads is necessary for downstream analysis, the `read_group` entity requires more fields than the other Biospecimen entities (`sample`, `portion`, `analyte`, `aliquot`).

Submitting a [__Read Group__](https://gdc-docs.nci.nih.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=read_group) entity requires:

* __`submitter_id`:__ A unique key to identify the `read_group`
* __`aliquot.submitter_id`:__ The unique key that was used for the `aliquot`, links the `read_group` to the `aliquot`
* __`experiment_name`:__ The name of the experiment
* __`is_paired_end`:__ If the reads are paired-end (Boolean value: true/false)
* __`library_name`:__ The name of the library  
* __`library_strategy`:__ The experimental strategy of the library
* __`platform`:__ Name of the platform used to sequence the data
* __`read_group_name`:__ The name of the `read_group`
* __`read_length`:__ The length of the reads (integer)
* __`sequencing_center`:__ The name of the center in which the sequencing was performed  

```JSON
{
    "type": "read_group",
    "submitter_id": "Blood-00001-aliquot_lane1_barcodeACGTAC_55",
    "experiment_name": "Resequencing",
    "is_paired_end": true,
    "library_name": "Solexa-34688",
    "library_strategy": "WXS",
    "platform": "Illumina",
    "read_group_name": "205DD.3-2",
    "read_length": 75,
    "sequencing_center": "BI",
    "aliquots":
        {
            "submitter_id": "Blood-00021-aliquot55"
        }    
}

```
```TSV
type	submitter_id	experiment_name	is_paired_end	library_name	library_strategy	platform	read_group_name	read_length	sequencing_center	aliquots.submitter_id
read_group	Blood-00001-aliquot_lane1_barcodeACGTAC_55	Resequencing	true	Solexa-34688	WXS	Illumina	205DD.3-2	75	BI	Blood-00021-aliquot55
```

[![GDC Data Model 5](images/DataModel-5.jpg)](images/DataModel-5.jpg "Click to see the full image.")

### Submitted Aligned-Reads Submission

Before the BAM file can be submitted, the GDC requires that the user provides information about the `submittable_data_file` (in this case, the BAM file itself).  This includes file-specific data needed to validate the file and assess which analyses should be performed.  

Submitting a [__Submitted Aligned-Reads__](https://gdc-docs.nci.nih.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=submitted_aligned_reads) entity requires:

* __`submitter_id`:__ A unique key to identify the `submitted_aligned_reads`
* __`read_groups.submitter_id`:__ The unique key that was used for the `read_group`, links the `submitted_aligned_reads` to the `read_group`
* __`data_category`:__ A broad categorization of the data file contents
* __`data_format`:__ The data file format
* __`data_type`:__ The specific contents of the data file (must be "Aligned Reads")
* __`experimental_strategy`:__ The sequencing strategy used to generate the file  
* __`file_name`:__ The name of the file
* __`file_size`:__ The size of the file in bytes
* __`md5sum`:__ The 128-bit hash value expressed as a 32 digit hexadecimal number

```JSON
{
    "type": "submitted_aligned_reads",
    "submitter_id": "Blood-00001-aliquot_lane1_barcodeACGTAC_55.bam",
    "data_category": "Raw Sequencing Data",
    "data_format": "BAM",
    "data_type": "Aligned Reads",
    "experimental_strategy": "WGS",
    "file_name": "test.bam",
    "file_size": 38,
    "md5sum": "aa6e82d11ccd8452f813a15a6d84faf1",
    "read_groups": [
        {
            "submitter_id": "Blood-00001-aliquot_lane1_barcodeACGTAC_55"
        }
    ]
}
```
```TSV
type	submitter_id	data_category	data_format	data_type	experimental_strategy	file_name	file_size	md5sum	read_groups.submitter_id#1
submitted_aligned_reads	Blood-00001-aliquot_lane1_barcodeACGTAC_55.bam	Raw Sequencing Data	BAM	Aligned Reads	WGS	test.bam	38	aa6e82d11ccd8452f813a15a6d84faf1	Blood-00001-aliquot_lane1_barcodeACGTAC_55
```

__Note:__ Because there can be many `read_groups` included in one `submitted_aligned_reads` file, the '\#1' is appended to the `read_groups.submitter_id` field in the TSV.  This relationship can be expressed with a list object (square brackets) in a JSON formatted file.   

### Uploading the BAM File to the GDC

The BAM file can be uploaded now that the GDC has information about it.  This can be performed using either the GDC Data Transfer Tool or the API.  

__GDC Data Transfer Tool:__ A file can be uploaded using its UUID (which can be retrieved from the portal or API) once it is registered. The following command can be used to upload the file:

```Shell
gdc-client upload --project-id GDC-INTERNAL --identifier a053fad1-adc9-4f2d-8632-923579128985 -t $token -f $path_to_file
```   

Additionally a manifest can be downloaded from the Submission Portal and passed to the Data Transfer Tool, this will allow for the upload of more than one `submittable_data_file`:

```Shell
gdc-client upload -m manifest.yml -t $token
```
__API Upload:__  `submittable_data_files` can be uploaded through the API by using the `/submission/program/project/files` endpoint.  The following command would be a typical file upload:  

```Shell
curl --request PUT --header "X-Auth-Token: $token" https://gdc-api.nci.nih.gov/v0/submission/GDC/INTERNAL/files/6d45f2a0-8161-42e3-97e6-e058ac18f3f3 -d@test.bam

```

For more details on how to upload a `submittable_data_file` to a project see the [API Users Guide](API/Users_Guide/Submission/) and the [Data Transfer Tool Users Guide](Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/).  


## Adding Clinical or Metadata Files to an Entity.

Typically a submission project will include additional information about a `case` such as `demographic`, `diagnosis`, or `exposure` data. The project also may include additional information about the experimental procedures used to produce the data.  The clinical data will be associated with the `case` entity and metadata can be associated with the `read_group` and `submitted_aligned_reads` entities.  See the diagram below for how these data types can be associated with each entity.

[![GDC Data Model Meta](images/Data_Model_Meta.jpg)](images/Data_Model_Meta.jpg "Click to see the full image.")

### Submitting a Demographic Entity to a Case

The `demographic` entity contains information that characterizes the `case` entity, which is the patient in most instances.  

Submitting a `demographic` entity requires:

* __`submitter_id`:__ A unique key to identify the `demographic` entity
* __`cases.submitter_id`:__ The unique key that was used for the `case`, links the `demographic` entity to the `case`
* __`ethnicity`:__ An individual's self-described identity as Hispanic or Latino
* __`gender`:__ The gender of the individual
* __`race`:__ The race of the individual based on values used by the U.S. Census Bureau
* __`year_of_birth`:__ The calendar year in which an individual was born.  

```JSON
{
    "type": "demographic",
    "submitter_id": "GDC-INTERNAL-000055-DEMOGRAPHIC-1",
    "cases": {
        "submitter_id": "GDC-INTERNAL-000055"
    },
    "ethnicity": "not hispanic or latino",
    "gender": "male",
    "race": "asian",
    "year_of_birth": "1946"
}
```
```TSV
type	cases.submitter_id	ethnicity	gender	race	year_of_birth
demographic	GDC-INTERNAL-000055	not hispanic or latino	male	asian	1946
```

### Submitting an Exposure Entity to a Case

Submitting an [Exposure](https://gdc-docs.nci.nih.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=exposure) entity does not require any information besides a link to the `case` and a `submitter_id`.  This entity contains information such as:  

* __`alcohol_history`:__ If the individual has consumed at least 12 drinks of any kind of alcoholic beverage in their lifetime
* __`alcohol_intensity`:__ An individual's self-reported level of alcohol use
* __`bmi`:__ Body mass divided by body height squared
* __`cigarettes_per_day`:__ The average number of cigarettes smoked per day
* __`height`:__ The height of the individual in cm  
* __`weight`:__ The weight of the individual in kg
* __`years_smoked`:__ The number of years an individual has been smoking

```JSON
{
    "type": "exposure",
    "submitter_id": "GDC-INTERNAL-000055-EXPOSURE-1",
    "cases": {
        "submitter_id": "GDC-INTERNAL-000055"
    },
    "alcohol_history": "yes",
    "bmi": 27.5,
    "cigarettes_per_day": 20,
    "height": 190,
    "weight": 100,
    "years_smoked": 5
}
```
```TSV
type	submitter_id	cases.submitter_id	alcohol_history	bmi	cigarettes_per_day	height	weight	years_smoked
exposure	GDC-INTERNAL-000055-EXPOSURE-1	GDC-INTERNAL-000055	yes	27.5	20	190	100	5
```

### Submitting an Experiment Metadata Entity to a Read Group

The `experiment_metadata` entity contains information about the experiment that was performed to produce each `read_group`. Unlike the previous two entities outlined, only information about the `experiment_metadata` file itself (SRA XML) is indexed and the `experiment_metadata` file is submitted in the same way that a BAM file would be submitted.

Submitting an __Experiment Metadata__ entity requires:

* __`submitter_id`:__ A unique key to identify the `experiment_metadata` entity
* __`read_groups.submitter_id`:__ The unique key that was used for the `read_group`, links the`experiment_metadata` entity to the `read_group`
* __`data_category`:__ Broad categorization of the data file
* __`data_format`:__ Format of the data file (must be SRA XML)
* __`data_type`:__ Specific contents of the data file
* __`file_name`:__ The name of the file  
* __`file_size`:__ The size of the file  
* __`md5sum`:__ 128-bit hash value expressed as a 32 digit hexadecimal number   

```JSON
{
    "type": "experiment_metadata",
    "submitter_id": "Blood-00001-aliquot_lane1_barcodeACGTAC_55-EXPERIMENT-1",
    "cases": {
        "submitter_id": "Blood-00001-aliquot_lane1_barcodeACGTAC_55"
    },
    "data_category": "Sequencing Data",
    "data_format": "SRA XML",
    "data_type": "Experiment Metadata",
    "file_name": "Experimental-data.xml",
    "file_size": 65498,
    "md5sum": "d79997e4de03b5a0311f0f2fe608c11d",
}
```

## Strategies for Submitting in Bulk


Each submission in the previous two sections was broken down by component to demonstrate the GDC Data Model structure.

### Registering a BAM File: One Step

Registering a BAM file (or any other type) can be performed in one step by including all of the entities, from `case` to `submitted_aligned_reads`, in one file.  See the example below:

```JSON
[{
    "type": "case",
    "submitter_id": "GDC-INTERNAL-000055",
    "projects": {
        "code": "INTERNAL"
    }
},
{
    "type": "sample",
    "cases": {
        "submitter_id": "GDC-INTERNAL-000055"
    },
    "sample_type": "Blood Derived Normal",
    "submitter_id": "Blood-00001_api55"
},
{
    "type": "portion",
    "submitter_id": "Blood-portion-000055",
    "samples": {
        "submitter_id": "Blood-00001_api55"
    }
},
{
    "type": "analyte",
    "portions": {
        "submitter_id": "Blood-portion-000055"
    },
    "analyte_type": "DNA",
    "submitter_id": "Blood-analyte-000055"
},
{
    "type": "aliquot",
    "submitter_id": "Blood-00021-aliquot55",
    "analytes": {
        "submitter_id": "Blood-analyte-000055"
    }
},
{
    "type": "read_group",
    "submitter_id": "Blood-00001-aliquot_lane1_barcodeACGTAC_55",
    "experiment_name": "Resequencing",
    "is_paired_end": true,
    "library_name": "Solexa-34688",
    "library_strategy": "WXS",
    "platform": "Illumina",
    "read_group_name": "205DD.3-2",
    "read_length": 75,
    "sequencing_center": "BI",
    "aliquots":
        {
            "submitter_id": "Blood-00021-aliquot55"
        }    
},
{
    "type": "submitted_aligned_reads",
    "submitter_id": "Blood-00001-aliquot_lane1_barcodeACGTAC_55.bam",
    "data_category": "Raw Sequencing Data",
    "data_format": "BAM",
    "data_type": "Aligned Reads",
    "experimental_strategy": "WGS",
    "file_name": "test.bam",
    "file_size": 38,
    "md5sum": "aa6e82d11ccd8452f813a15a6d84faf1",
    "read_groups": [
        {
            "submitter_id": "Blood-00001-aliquot_lane1_barcodeACGTAC_55"
        }
    ]
}]
```

All of the entities are placed into a JSON list object:

`[entity-1, entity-2, entity-3]`

The entities need not be in any particular order as they are validated together.

__Note:__ For this type of submission, a tab-delimited format is not recommended due to the inability of this format to accommodate multiple 'types' in one row.  

### Registering Numerous Cases

The GDC understands that submitters will have projects that comprise more entities than would be reasonable to individually parse into JSON formatted files. Additionally, many investigators store large amounts of data in a tab-delimited format (TSV).  For instances like this, we recommend parsing all entities of the same type into separate TSVs and submitting them on a type-basis.  

For example, a user may want to submit 100 Cases associated with 100 `samples`, 100 `portions`, 100 `analytes`, 100 `aliquots`, and 100 `read_groups`. Constructing and submitting 100 JSON files would be tedious and difficult to organize. Submitting one `case` TSV, one `sample` TSV, and the rest would not only require six TSVs but can easily be formatted in programs such as Microsoft Excel or Google Spreadsheets.  

See the following example TSV files:

* [Cases.tsv](Cases.tsv)
* [Samples.tsv](Samples.tsv)
* [Portions.tsv](Portions.tsv)
* [Analytes.tsv](Analytes.tsv)
* [Aliquots.tsv](Aliquots.tsv)
* [Readgroups.tsv](Readgroups.tsv)
