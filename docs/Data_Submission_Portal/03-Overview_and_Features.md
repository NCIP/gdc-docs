# Overview

[![GDC Data Submission Portal Sitemap](images/GDC_Submission_Portal_Sitemap.png)](images/GDC_Submission_Portal_Sitemap.png "Click to see the full image.")

The above image shows the GDC Data Submission Portal sitemap. The Portal comprises different sections available through a toolbar on the upper level of the page and a navigation panel on the left-hand side of the page.

[![GDC Data Submission Toolbar](images/GDC_Submission_Toolbar.png)](images/GDC_Submission_Toolbar.png "Click to see the full image.")

The top toolbar provides access to the following elements:

* __GDC Identifying Image__: Clicking on GFDC Identifying image will redirect to the submission dashboard.
* __Search__: A context sensitive search provides a mechanism to search through different pages.
* __Profile__: The user profile icon provides access to the token and to a logout feature.

# Submission Portal Layout

[![GDC Data Submission Portal Layout](images/GDC_Submission_Portal_Layout.png)](images/GDC_Submission_Portal_Layout.png "Click to see the full image.")

The GDC Data Submission Portal main window, visible when navigating through GDC Submission Portal, has been divided into three separate panels:

* __Navigation panel__: A panel composed of links to different sections of the application.
* __List panel__: A table listing items corresponding to the current selection.
* __Details panel__: Details and previous activities related to a particular entity selected in the results panel.

[![GDC Data Submission Portal Results Panel](images/GDC_Submission_Results_Panel.png)](images/GDC_Submission_Results_Panel.png "Click to see the full image.")

# Navigation Panel

[![GDC Data Submission Portal Navigation Panel](images/GDC_Submission_Navigation.png)](images/GDC_Submission_Navigation.png "Click to see the full image.")

The navigation panel provides access to the following elements (detailed in subsequent sections of the documentation):

* __SUBMIT__: Access to the submission wizard. This button is only available when a single project is selected.
* __Dashboard__: Shortcut to the project's dashboard.
* __Cases / All Cases__: Access to the list of cases in the project.
* __Biospecimen / Samples__: List of tissues or specimens obtained from patients
* __Biospecimen / Portions__: Portion of a sample identified as one single entity
* __Biospecimen / Analytes__: Experimental test material extracted from a portion of sample
* __Biospecimen / Aliquots__: Portion of an analyte used for experimental assay
* __Data / Data Read Groups__:  A read group is a type of data bundle that contains molecular data. A Data Bundle is set of files with associated metadata.
* __Reports / Case Overview__: Reports related to submission
* __Transactions__: List all actions happening on a project (submit data and release).

# Results Panel

The results panel lists items in a table based on the element selected in the navigation panel or in the results of a search.

Columns in the table vary depending of the selected element in the navigation panel.

[![GDC Data Submission Portal Navigation Panel](images/GDC_Submission_Results_Panel.png)](images/GDC_Submission_Results_Panel.png "Click to see the full image.")

The table view supports pagination and columns sorting and, if applicable, alerts are displayed to notify that additional attention is required, such as the those shown below.

[![Cases with no clinical data](images/GDC_Submission_Cases_with_no_Clinical_Data.png)](images/GDC_Submission_Cases_with_no_Clinical_Data.png "Click to see the full image.")

# Details and Activity Panel

The third panel, on the right-hand side of the page, provides more details about a selected entity in two sub-sections, _'Details'_ and _'Activity'_.

The exact content of those sub-sections depends on the type of entity (e.g. case, sample, transaction) of the selection in the navigation panel.

[![Transaction details](images/GDC_Submission_Transaction_Details.png)](images/GDC_Submission_Transaction_Details.png "Click to see the full image.")

The Activity section lists all past transactions associated with an entity, clicking on a transaction will redirect the user to this particular transaction.

# Filtering and Searching

The search menu available at the top of the screen is context-sensitive.

For example, when visiting the _'Samples'_ page, the search menu will display the prefix _'Within Samples'_ and only search within Samples.

[![Transaction details](images/GDC_Submission_Search_within_Samples.png)](images/GDC_Submission_Search_within_Samples.png "Click to see the full image.")
