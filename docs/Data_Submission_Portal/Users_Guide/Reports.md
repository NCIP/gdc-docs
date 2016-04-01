# Reports

## Overview

The GDC Data Submission Portal provides access to the different reports listed below. Access to the reports is granted to users based on their project permissions.

The reports are available from the [Homepage](Homepage.md) and the [Project Dashboard](Dashboard.md).

[![GDC Submission Reports](images/GDC_Submission_Reports.png)](images/GDC_Submission_Reports.png "Click to see the full image.")

## Case Overview Report

This report provides an overview of all cases associated with one or more projects and is generated daily.

The report is downloadable in TSV format.

```TSV
Program	Project	dbGaP Study ID	# Cases Released	# Cases Unreleased	# Cases Expected	Date Expected for Release	Initial Upload Date	Project Last Update Date	# Cases Registered in dbGaP	# Cases with Clinical	# Cases with Sample	# Cases with Aliquot	# Cases with Portion	# Cases with Analyte	# Cases with WGS Data	# Cases with WXS Data	# Cases with RNA-Seq Data	# Cases with ChIP-Seq Data	# Cases with MiRNA-Seq Data	# Cases with Bisulfite-Seq Data	# Cases with Validation Data	# Cases with Amplicon Data	# Cases with Other Data
TCGA	TEST	TEST	256	0	Not Available	Not Available	2015-11-10 10:40:15.689597-06:00	2015-11-10 10:40:15.689597-06:00	256	0	11	10	10	10	0	4	0	0	0	0	0	0	0
TCGA	DEV1	phs000178	1	15	Not Available	Not Available	2015-11-18 09:27:36.990449-06:00	2015-11-18 09:27:36.990449-06:00	16	0	16	16	16	16	0	0	0	0	0	0	0	0	0
TCGA	DEV2	phs000178	1	7	Not Available	Not Available	2015-11-17 20:14:02.194718-06:00	2015-11-17 20:14:02.194718-06:00	8	0	8	8	0	0	0	0	0	0	0	0	0	0	0
TCGA	DEV3	phs000178	64	0	Not Available	Not Available	2015-12-11 15:41:29.760763-06:00	2015-12-11 15:41:29.760763-06:00	64	0	25	25	0	0	0	2	0	0	0	0	0	0	0
```

Previous versions of the report are available for comparison. To download a previous version, the user should click on the "Previous Snapshots" button. This will open a calendar to allow the user to select the date of the report for download.

## Data Validation Report

The Data Validation Report provides a live view of quality metrics coming from data validation performed by the GDC when validating a file upload through the GDC Data Transfer Tool.

The report can be downloaded in TSV format.

```TSV
program	project	dbgap_accession_number	file_format_fail	file_size_fail	md5sum_fail	processed	processing	registered	submitted	uploaded	uploading	validated	validating
TCGA	DEV3		0	1	0	0	0	3	0	0	0	1	0
```
