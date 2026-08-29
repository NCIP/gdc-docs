# Cohort Builder

The Cohort Builder is a good starting point for users looking to gather information for a specific disease, project, or group of patients. Building a cohort allows users to download files, perform analyses, and query metadata for the same group of cases in multiple sections of the GDC Data Portal. This section will cover the process of building a cohort and downstream actions will be documented in their respective sections.

> **NOTE:** Filters within the Cohort Builder are applied to the cases in your cohort. If you wish to target specific types of files for download, use the filters within the Repository.

The Cohort Builder can be accessed in one of the following ways:

* Selecting the Cohort Builder link in the GDC Data Portal header

[![Cohort Builder In Header](images/ToolLinksInHeader.png)](images/ToolLinksInHeader.png "Click to see the full image.")

* Selecting the play button on the Cohort Builder card in the Analysis Center

[![Cohort Builder Play Buttons](images/CohortBuilderInAnalysisCenter.png)](images/CohortBuilderInAnalysisCenter.png "Click to see the full image.")

## Cohort Builder Panel

[![Cohort Builder Panel](images/CohortBuilderPanel.png)](images/CohortBuilderPanel.png "Image of Cohort Builder Panel. Click to see the full image.")

The Cohort Builder tool will be displayed as a panel in the Analysis Center and is used to filter the current cohort to a specific set of cases. The current cohort is always displayed in the [main toolbar](getting_started.md#main-toolbar) and can be changed from the main toolbar.

At the left side of the panel are a series of broad filter categories can be selected. Each filter category contains a set of specific filters within cohort builder cards that can be used to narrow your cohort to the desired set.

## Cohort Builder Cards
Each card within the Cohort Builder can be used to apply the corresponding filters on the current cohort. As filters are applied, they will be displayed on the [Query Expressions](getting_started.md#query-expressions) section.

Additional features can be accessed at the top right of each card's header to facilitate filtering:

* **Search**: the search icon can be selected to reveal or hide a search field for entering text to search within the values of the current card. This feature is only available when the values are enums.
* **Flip Card**: cards can be flipped to reveal or hide a summary chart. This feature is only available when the values can be meaningfully displayed as bar graphs.
* **Reset Card**: this button will reset any filtering that has been applied within the card

[![Cohort Builder Card Features](images/search_flip_reset_card.png)](images/search_flip_reset_card.png "Image of Cohort Builder Card Features. Click to see the full image.")

In addition, filters in each card can be sorted, either alphabetically or by the number of cases based on current filters, by selecting one of the two headers directly underneath the card title. The default sort is alphabetical order.

[![Cohort Builder Card Filter Sort](images/alpha_desc_card_filter_order.png)](images/alpha_desc_card_filter_order.png "Image of Cohort Builder Card Filter Sort. Click to see the full image.")

The first six (or fewer) filters are shown for each card, but can be expanded to show 20 filters at once by clicking the "+" button which also indicates the number of additional filters not in view. The expanded view can be toggled off by clicking the resulting "show less" button.

[![Cohort Builder Card Filter Expand](images/card_expand.png)](images/card_expand.png "Image of how to expand cohort builder card filters. Click to see the full image.")

### Biospecimen Filters
The filters within the "Biospecimen" category allow for cases that have certain types of biospecimens. For example, filtering in the "Tissue Type" card for "normal" will ensure that cases within your cohort have a normal tissue type. These filters may be useful for studies that require only cases for which a certain type of biospecimen is available.

### Available Data Filters

Toward the bottom of the list of filter categories, "Available Data" can be selected.  These filters differ from the other default filters as they allow for cases that have certain types of associated data files. For example, filtering in the "Experimental Strategy" card for "RNA-Seq" will only display cases in the active cohort that have associated RNA-Seq files. These filters may be useful for studies that require only cases for which a certain type of analysis was performed.

> **NOTE:** Some of the "Biospecimen" and "Available Data" filters are also available in the Repository, where they may be used to specify the types of files for browsing and download. For example, selecting "WXS" in the "Experimental Strategy" card will ensure that all cases in your cohort have WXS data. This, however, does not mean that your cohort will not also have other types of data, such as RNA-Seq or WGS. Likewise, selecting "normal" in the "Tissue Type" card narrows your cohort to those cases that have normal samples, but does not exclude cases with other types of samples. If the goal is to view or download only the files of a specific type, use the filters in the Repository rather than the filters located here in the Cohort Builder.

## Custom Filters ##

If a filter cannot be found within one of the categories, use the "Add a Custom Filter" button in the "Custom Filters" category to access any filters that are not displayed. Browse through the list of additional filters, or use the search box to search for filters by name. Once a filter is selected, it is then added to the "Custom Filters" category.  A custom filter can be removed from this category by choosing the "X" at the top right of the filter card.

[![Custom Filter](images/CustomFilter.png)](images/CustomFilter.png "Image of Custom Filter search box. Click to see the full image.")

Filters that exist in the GDC but do not have any cases that have a value for the filter can also be removed from the "Custom Filters" list by selecting the "Only show properties with values" box.

[![Custom Filter Only Show Props Vals](images/only_show_props_vals.png)](images/only_show_props_vals.png "Image of toggle for properties with values for custom filters in cohort builder. Click to see the full image.")

## Cohort Builder Search

The Cohort Builder includes the ability to search across all the cards within it. This feature is located on the right of the Cohort Builder header.

[![Cohort Builder Header](images/CohortBuilderHeader.png)](images/CohortBuilderHeader.png "Image of Cohort Builder Header. Click to see the full image.")

As a search term is entered, the Cohort Builder Search feature will display a list of properties that contain matching results. When a result is moused over, additional information is displayed to its left, including a description of the property and a list of values that match the search term.

[![Cohort Builder Search Results](images/CohortBuilderSearchResults.png)](images/CohortBuilderSearchResults.png "Image of Cohort Builder Search Results. Click to see the full image.")

When a result is selected, the card corresponding to the selected result will be displayed. If there are values that match the search term, the card's search field will be automatically populated with the search term.

## Closing the Cohort Builder

Once a custom cohort is built and filtering is complete, users can close the Cohort Builder and use the custom cohort with other tools.

To close the Cohort Builder panel and display all the tools within the Analysis Center, click on the "X" button on the left of the Cohort Builder header.

[![Cohort Builder Header](images/CohortBuilderHeader.png)](images/CohortBuilderHeader.png "Image of Cohort Builder Header. Click to see the full image.")

Alternatively, users can select the Analysis Center link or any of the other links on the GDC Data Portal header to close the Cohort Builder.

[![Header with Links to Tools](images/ToolLinksInHeader.png)](images/ToolLinksInHeader.png "Click to see the full image.")

Changes made to the cohort with the Cohort Builder will persist through the other sections of the GDC Data Portal.

Users can then perform the following actions:

* Download files associated with the cohort from the [Repository](Repository.md)
* Analyze data from the cohort in the [Analysis Center](analysis_center.md)

## Cohort Types

Depending on how they are modified or created, cohorts can have different types of filters and thus behave differently after a data release with regard to the cases they contain. The following are the types of filters cohorts can have:

* __Custom Queries__
* __Specific List of Cases__

### Custom Queries

A very common way to modify a cohort is by using the filters in the Cohort Builder. Using the filters in the Cohort Builder to build a cohort will create a cohort with custom queries (see note below for an exception). Examples of cohorts with custom queries would be:

1. All cases in the `TCGA-BRCA` project.
1. Cases with a primary site of `bronchus and lung` and at least `30` pack years smoked.

Cohorts based on custom queries will change depending on the data available. For example, if a data release adds cases to the TCGA-BRCA project, the first cohort example will include the new cases automatically and increase in size. 

The query expression section will display these custom queries with information about the properties and values that were applied as filters to the cohort.

[![Query Expressions section with TCGA-BRCA cohort](images/QueryExpressionsTCGABRCA.png)](images/QueryExpressionsTCGABRCA.png "Click to see the full image.")

[![Query Expressions section with lung cancer cohort with at least 30 pack years smoked](images/QueryExpressionsLungPackYears.png)](images/QueryExpressionsLungPackYears.png "Click to see the full image.")

**NOTE:** The Case ID filter in the Cohort Builder will result in a cohort based on a specific list of cases.

### Specific List of Cases

Cohorts can also be based on a list of specific cases. These cases do not necessarily share common properties and can comprise any group of released cases. Data releases that happen after these cohorts are created will not add additional cases to the cohort, but could subtract cases if there were redactions. A common way to create cohorts based on a list of specific cases is to use the Import New Cohort function in the Cohort Bar. Another common way is to create the cohorts from one of the many analysis tools available in the GDC (e.g. Clinical Data Analysis, Mutation Frequency, or Set Operations).

The query expression section will display a case UUID with the Case ID property if the filter has only 1 specific case. Otherwise, it will display the number of cases in the list.

[![Query Expressions section with 1 specific case](images/QueryExpressions1Case.png)](images/QueryExpressions1Case.png "Click to see the full image.")

[![Query Expressions section with many specific cases](images/QueryExpressionsManyCases.png)](images/QueryExpressionsManyCases.png "Click to see the full image.")

> **NOTE:** 
If an imported cohort was originally created by exporting a cohort with custom queries, it will still result in a cohort with specific cases. The export function saves a list of cases, but does not preserve the custom queries used to filter for those cases.