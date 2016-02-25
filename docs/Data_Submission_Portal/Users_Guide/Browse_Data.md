# Browse Data

## Overview

The _"Browse"_ tab provides access to all of a project's content. Most content is driven by the [GDC Dictionary](../../Dictionary/index.md) and the interface is dynamically generated to accomodate the content.

Please refer to the [GDC Dictionary](../../Dictionary/index.md) for more details about dictionary-generated fields, columns and filters.

All static components are detailed in this User Guide.

[![GDC Submission Cases Default View](images/GDC_Submission_Cases_Default.png)](images/GDC_Submission_Cases_Default.png "Click to see the full image.")

## Main interface elements

### Upload

An [upload](Upload_Data.md) button is available at the top-left section of the left panel. Depending of the project status, this will allow data to be [uploaded](Upload_Data.md) to the GDC. 

### Filters

A wide set of filters are available for the user to select the type of entity to be displayed. Most filters are dynamically created based on the [GDC Dictionary](../../Dictionary/index.md), but some have been added to provide additional insight about the data.

Hard-coded filters are:

|Filter|Description|
| --- | --- |
| All Cases | Display all cases associated with the project |
| Missing Clinical Data | Display only cases with missing Clinical Data |
| Missing Samples Data | Display only cases with missing Samples Data|
| Submitted Files | List files uploaded to GDC |
| Transactions | List all transactions associated to the project |
| Annotations | List all annotations associated to the project |

#### Case Filters

Cases can be accessed from the menu through multiple filters.

|Filter|Description|
| --- | --- |
| All Cases | Display all cases associated with the project |
| Missing Clinical Data | Display only cases with missing Clinical Data |
| Missing Samples Data | Display only cases with missing Samples Data|

The number of cases corresponding to each filter is displayed on the right side of the filter name.

#### Submitted Files

ADD CONTENT

#### Transactions

More details about transactions can be found in the [Transactions](Transactions.md) section of the documentation.

#### Annotations

ADD CONTENT


### List View

The list view is a paginated list of all entities corresponding to the selected filter.

On the top-left section of the screen, the user can download data about all entities associated with the selected filter.

### Details Panel

Clicking on a case will open the details panel. Data in this panel is broken down into multiple sections and driven by the [GDC Dictionary](../../Dictionary/index.md)

[![GDC Submission Case Details](images/GDC_Submission_Cases_Details.png)](images/GDC_Submission_Cases_Details.png "Click to see the full image.")

Navigation between those sections can be done either by scrolling down or by clicking on the section icon on the left side of the details panel.

[![GDC Submission Cases Details Navigation](images/GDC_Submission_Cases_Details_Navigation.png)](images/GDC_Submission_Cases_Details_Navigation.png "Click to see the full image.")


#### Related Entities

Table listing all entities, grouped by type, related to the selected case.

This table contains the following columns.

|Column|Description|
| --- | --- |
| Category | Category of the Entity (Clinical, Biospecimen, Experiment Data)  |
| Type | Type of entity (based on dictionary)  |
| Count | Number of occurences of an entity of this type |

Clicking on the count will open a list page listing those entities.


#### Transactions

Lists the 10 most recent transactions associated with this entity ordered by date. Clicking on a transaction will open it in the list page.
