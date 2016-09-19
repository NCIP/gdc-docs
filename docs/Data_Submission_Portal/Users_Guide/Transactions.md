# Transactions

## Overview

The transactions tab lists all transactions which have happened on the project since its creation.

[![GDC Submission Transactions](images/GDC_Submission_Transactions.png)](images/GDC_Submission_Transactions.png "Click to see the full image.")

The types of transactions are the following:

* __Upload__: The user uploads data to the project workspace.
* __Review__: The user reviews the project before submitting data to the GDC.
* __Open__: The user re-opens the project if it was under review. This will allow the upload of new data to the project workspace.
* __Submit__: The user submits data to the GDC. This trigger the data harmonization process.
* __Release__: The user releases data to the GDC Data Portal and other GDC data access tools.

## Transactions List View

The transactions list view displays the following information:

|Column|Description|
| --- | --- |
| ID | Identifier of the transaction |
| Type | Type of the transaction (see list transaction types in the previous section)|
| Step | XXXXXXXXXXXXXXX |
| DateTime | Date the transaction was initiated |
| User | XXXXXXXXXXXXXXX |
| Status | XXXXXXXXXXXXXXX |
| Commit/Discard | Two buttons that appear when data has been uploaded using the 'dry run' feature.  This allows for validated data to be incorporated into the project or discarded. |

## Transaction Filters

Choosing the radio buttons at the top of the table allows the transactions to be filtered by those that are in progress, to be committed, succeeded, failed, or discarded. The drop-down menu also allows for the transactions to be filtered by type.  


## Transactions Details

Clicking on a transaction will open the details panel. Data in this panel is organized into multiple sections including actions, details, types, and documents as described below.

[![GDC Submission Transactions](images/GDC_Submission_Transactions_Details.png)](images/GDC_Submission_Transactions_Details.png "Click to see the full image.")

Navigation between the sections can be performed by either scrolling down or by clicking on the section icon displayed on the left side of the details panel.

### Actions

The Actions section allows a user to perform an action for transactions that provide actions. For example, if a user uploads read groups and file metadata, a corresponding manifest file will be available for download from the transaction. This manifest is used to upload the actual files through the [GDC Data Transfer Tool](https://gdc.nci.nih.gov/access-data/gdc-data-transfer-tool).

[![GDC Submission Transactions Details Action](images/GDC_Submission_Transactions_Details_Action.png)](images/GDC_Submission_Transactions_Details_Action.png "Click to see the full image.")

### Details

The Details section provides details about the transaction itself, such as its project, type, and number of affected cases.

[![GDC Submission Transactions Details](images/GDC_Submission_Transactions_Details_Details.png)](images/GDC_Submission_Transactions_Details.png "Click to see the full image.")

### Types

The Types section lists the type of files submitted and the number of affected cases and entities.

[![GDC Submission Transactions Types](images/GDC_Submission_Transactions_Details_Types.png)](images/GDC_Submission_Transactions_Details_Types.png "Click to see the full image.")

### Documents

The Documents section lists the files submitted during the transaction.
The user can __download the original files from the transaction__ or a report detailing the original files in tab-delimited format.

[![GDC Submission Transactions Documents](images/GDC_Submission_Transactions_Details_Documents.png)](images/GDC_Submission_Transactions_Details_Documents.png "Click to see the full image.")
