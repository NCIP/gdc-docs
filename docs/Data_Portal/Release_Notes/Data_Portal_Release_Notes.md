# Data Portal Release Notes

| Version | Date |
|---|---|
| [v2.3.0](Data_Portal_Release_Notes.md#release-230) | October 29, 2024 |
| [v2.2.0](Data_Portal_Release_Notes.md#release-220) | June 26, 2024 |
| [v2.1.0](Data_Portal_Release_Notes.md#release-210) | April 30, 2024 |
| [v2.0.0](Data_Portal_Release_Notes.md#release-200) | February 8, 2024 |
| [v1.30.4](Data_Portal_Release_Notes.md#release-1304) | May 11, 2023 |
| [v1.30.0](Data_Portal_Release_Notes.md#release-1300) | July 8, 2022 |
| [v1.29.0](Data_Portal_Release_Notes.md#release-1290) | August 23, 2021 |
| [v1.28.0](Data_Portal_Release_Notes.md#release-1280) | May 17, 2021 |
| [v1.25.1](Data_Portal_Release_Notes.md#release-1251) | August 14, 2020 |
| [v1.25.0](Data_Portal_Release_Notes.md#release-1250) | July 2, 2020 |
| [v1.24.1](Data_Portal_Release_Notes.md#release-1240) | March 10, 2020 |
| [v1.23.1](Data_Portal_Release_Notes.md#release-1231) | December 10, 2019 |
| [v1.23.0](Data_Portal_Release_Notes.md#release-1230) | November 6, 2019 |
| [v1.22.0](Data_Portal_Release_Notes.md#release-1220) | July 31, 2019 |
| [v1.21.0](Data_Portal_Release_Notes.md#release-1210) | June 5, 2019 |
| [v1.20.0](Data_Portal_Release_Notes.md#release-1200) | April 17, 2019 |
| [v1.19.0](Data_Portal_Release_Notes.md#release-1190) | February 20, 2019 |
| [v1.18.0](Data_Portal_Release_Notes.md#release-1180) | December 18, 2018 |
| [v1.17.0](Data_Portal_Release_Notes.md#release-1170) | November 7, 2018 |
| [v1.16.0](Data_Portal_Release_Notes.md#release-1160) | September 27, 2018 |
| [v1.15.0](Data_Portal_Release_Notes.md#release-1150) | August 23, 2018 |
| [v1.14.0](Data_Portal_Release_Notes.md#release-1140) | June 13, 2018 |
| [v1.13.0](Data_Portal_Release_Notes.md#release-1130) | May 21, 2018 |
| [v1.12.0](Data_Portal_Release_Notes.md#release-1120) | February 15, 2018 |
| [v1.11.0](Data_Portal_Release_Notes.md#release-1110) | December 21, 2017 |
| [v1.10.0](Data_Portal_Release_Notes.md#release-1100) | November 16, 2017 |
| [v1.9.0](Data_Portal_Release_Notes.md#release-190) | October 24, 2017 |
| [v1.8.0](Data_Portal_Release_Notes.md#release-180)| August 22, 2017 |
| [v1.6.0](Data_Portal_Release_Notes.md#release-160) | June 29, 2017 |
| [v1.5.2](Data_Portal_Release_Notes.md#release-152) | May 9, 2017 |
| [v1.4.1](Data_Portal_Release_Notes.md#release-141) | October 31, 2016 |
| [v1.3.0](Data_Portal_Release_Notes.md#release-130) | September 7, 2016 |
| [v1.2.0](Data_Portal_Release_Notes.md#release-120) | August 9th, 2016 |
| [v1.1.0](Data_Portal_Release_Notes.md#release-110) | June 1st, 2016 |
| [v1.0.1](Data_Portal_Release_Notes.md#release-101) | May 18, 2016 |

---
## Release 2.3.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  October 29, 2024

### New Features and Changes
* __Cohort Builder__:
    * Custom filters now display their parent category name. <!--PEAR-1083-->
    * Filter cards in classification categories have been moved to the General Diagnosis category or the new Disease Specific Classifications category. <!--PEAR-1989-->
    * The Years-Days toggle has been removed for Age at Index. <!--PEAR-1994-->
    * Cards with number range filters are better aligned with other cards in the same row. <!--PEAR-1881-->
    * Search results have been improved to display more relevant searches first. <!--PEAR-2013-->
    * The Best Overall Response card has been moved to the first card in the Treatment category. <!--PEAR-2123-->
    * Descriptions are now available for Other Clinical Attribute properties when adding custom filters. <!--PEAR-2226-->
    * UICC Clinical and Pathologic Stage filter cards have been added to the General Diagnosis category, and Specimen Type has been added to the Biospecimen category. <!--PEAR-2205/2042-->
    * The Enneking MSTS Stage card and Composition card have been removed from the default cards. <!--PEAR-2205/2042-->
    * Cards for filtering cohort by specific cases, mutated genes, and SSMs have been added. <!--PEAR-1269-->
    * Improved number range inputs in Cohort Builder by removing autofill behavior, adding Min/Max labels, and validating user inputs with error messages for out-of-range values. <!--PEAR-2030-->
* __General UX/UI Improvements__:
    * Filter cards for entering text and entering/uploading sets now have an appropriate maximum height. <!--PEAR-1311-->
    * Additional loading indicators have been added throughout the portal to indicate that information is still in the process of rendering. <!--PEAR-748-->
    * Total counts are now consistently displayed above tables. <!--PEAR-1896/2093-->
    * Row selection is now appropriately disabled for table rows containing 0 items in a set. <!--PEAR-1951-->
    * Styling for the tool cards in the Analysis Center has been standardized. <!--PEAR-1957-->
    * The search bar in the left panel within the __Clinical Data Analysis__ tool now remains fixed at the top of the page. <!--PEAR-1963-->
    * The message "No data for this field" will only be displayed when information for a filter card has been loaded. <!--PEAR-2036-->
    * Vertical alignment has been improved for tables that are displayed next to each other. <!--PEAR-2049--> 
    * Filter panels located on the left side of the __Projects__, __Repository__, and __Mutation Frequency__ tools will now extend up to the height of the tables in the tools. <!--PEAR-2110-->
    * Styling for survival plots has been improved for consistency. <!--PEAR-2176-->
    * Download icons have been standardized. <!--PEAR-2186-->
    * Text size has been increased for instructions in modals for selecting cohorts. <!--PEAR-2175-->
* __File Summary Page__:
    * The __Reference Genome__ section is no longer displayed for files that have not been processed with the reference genome. <!--PEAR-1967-->
    * The Case ID column is now displayed by default in the Annotations table. Additionally, the Case UUID column is no longer displayed by default. <!--PEAR-1968-->
    * Pagination has been added to the Read Groups table. <!--PEAR-2192-->
    * Sample Type has been removed from the Associated Cases/Biospecimens table and replaced with Tissue Type and Tumor Descriptor. <!--PEAR-2042-->
* __Repository__:
    * Stability improvements have been added. <!--PEAR-2022-->
    * The placement and design of the buttons to add custom filters and reset them have been updated. <!--PEAR-2059-->
* __Clinical Data Analysis__:
    * The y-axis of histograms will now only display integers for case counts. <!--PEAR-2087-->
    * The rounding of numbers displayed in the tool has been improved. <!--PEAR-1040-->
* __Case Summary Page__:
    * Information about Other Clinical Attributes has been added to the Clinical section. Additionally, deprecated properties have been removed from the Follow-Ups table. <!--PEAR-1983-->
    * Sample Type, Sample Type ID, and Composition have been removed from the Biospecimen tree's Samples table. Additionally, the table has been updated with the addition of Specimen Type. <!--PEAR-2042-->
* __ProteinPaint__:
    * Allows users to customize consequence colors and restore to default.
* __Gene Expression Clustering__:
    * The tool now displays the top 1000 variably expressed genes and first 1000 cases as the default plot.
    * The tool and Gene Expression API are now ~80% more performant.
    * Users can now change the plot color scheme to Blue-White-Red, Green-Black-Red, Green-Black-Red, and Blue-Black-Yellow, under the Clustering tab.
    * Z-score values are now capped to not exceed absolute values, under the Clustering tab.
    * Implemented support for adding user-saved custom gene sets, under the Genes tab.
    * Added support for adding “Overall Survival” as an annotation variable.
    * The tool now supports clicking the gene dendrogram to select genes and launch “Gene Set Overrepresentation Analysis.”
    * When screening user-defined gene sets, use a close-to-zero min_median_log2_uqfpkm parameter to keep more genes expressed at low level.
    * Allows the drag/drop of mutation and dictionary variable rows not used for clustering.
    * Improved performance of the top variably expressed genes query from the GDC API by directly submitting the case filters without first retrieving a list of cases.
* __OncoMatrix__:
    * Implemented support for adding user-saved custom gene sets, under the Genes tab.
    * Added support of gene expression rows along with mutation and dictionary variables.
    * When plotting gene expression genes as a variable, users can now edit the genes display as a z-score.
    * For cases that do not have expression data, the displayed gene expression variable track now displays them as blank values.
    * Added support for “Overall Survival” variable.
    * Allows the display of only dictionary variables with all genes removed.
    * Users can now click on a mutated gene and edit and customize its variant grouping.
* The tooltips for the __Survival Plot__ now display the time to death and the interval of last follow-up in both years and months. The downloaded TSV now includes the time value in years, months, and days, and the downloaded JSON now includes the time value in days. <!--PEAR-1961/2060-->
* The ability to reset all filters in the __Projects__, __Repository__, and __Mutation Frequency__ tools to their defaults has been added. <!--PEAR-1431-->
* Filters in the __Projects__, __Repository__, and __Mutation Frequency__ tools no longer reset when the composition of the active cohort has been changed. <!--PEAR-1856-->
* Filter cards in the __Projects__, __Repository__, and __Mutation Frequency__ tools can now be expanded and collapsed singly or all at once. <!--PEAR-2029-->
* Except for the Most Frequent Somatic Mutations table in the __Case Summary Page__, downloadeded JSON and TSV files now reflect the information displayed in the associated tables whenever search filters have been applied. <!--PEAR-1865/2190--> 
* A modal will now be displayed to inform users of any issues that occurred when saving sets and cohorts, and when exporting sets. <!--PEAR-1971/2141-->
* __Quick Search__'s accuracy has been improved to account for files that are no longer available. <!--PEAR-2082-->
* When genes or mutations are entered or uploaded for filtering in __Mutation Frequency__, other filters within the tool will be cleared. <!--PEAR-2133-->
* The ability to display a banner notifying users of government shutdowns has been added. <!--PEAR-2161-->
* __Slide Image Viewer__'s performance has been improved. <!--PEAR-1771-->

### Bugs Fixed Since Last Release
* __Section 508 Accessibility__:
    * Aria roles now contain the expected children. <!--PEAR-1669-->
    * Responsiveness for __Mutation Frequency__, all summary pages, and all table headers has been improved. <!--PEAR-1927/2130/2090-->
    * An equivalent alternative to the body plot on the home page is now available. <!--PEAR-1937-->
    * Aria labels have been made consistent with the displayed text in the __Query Expressions__ section. <!--PEAR-2117-->
    * Fixed input labels and color contrasts in __ProteinPaint__, __Gene Expression Clustering__, and __OncoMatrix__.
* __Cohort Builder__:
    * Fixed inconsistent behavior for number range cards when removing filters.  <!--PEAR-1379-->
    * The tool now ensures that cards are displayed by default when loaded, resolving issues where no cards appeared after using the browser's back button or clicking the Cohort Builder link. <!--PEAR-1808-->
    * Cards in the Treatment tab now display the correct case counts after a selection has been made. <!--PEAR-1830-->
    * Fixed issue where entering "0" in range cards would not persist after applying, affecting all range cards that accept "0" as a valid entry. <!--PEAR-2098-->
    * Updated the minimum value for the "Age at Diagnosis" range to "0" years, replacing the incorrect value of "-90" years. <!--PEAR-2106-->
    * Fixed an issue where number range filters can be added to a cohort multiple times and not be properly removed. <!--PEAR-2181-->
* __Cohort Bar__:
    * Resolved inconsistencies in the Discard Changes button and cohort status indicators, ensuring clearer behavior for unsaved cohorts and improved messaging for users. <!--PEAR-1499-->
    * The Metadata download now correctly includes entries when molecular filters are applied to the cohort, resolving the issue where the file was previously blank. <!--PEAR-2203-->
* __File Summary Page__:
    * Fixed an issue where the incorrect file version table could appear on the File Summary page after performing multiple searches. <!--PEAR-2057-->
    * Resolved an issue in the portal where navigating from a file page summary to its source files page summary and hitting the back button did not load the complete content of the initially searched file. <!--PEAR-2120-->
    * Addressed the issue where the action button icon in the Source Files table was not displayed on screens narrower than 1280px. <!--PEAR-2089-->
* __Authentication__:
    * Implemented error modal/banner for users without controlled data access when attempting to view "Download Token" in the GDC data portal. <!--PEAR-1024-->
    * Fixed continuous loading spinner in the header for users without access to controlled projects upon login or page navigation. <!--PEAR-1351-->
* __ProteinPaint__:
    * Improved the gene search results latency when a user presses the Enter key immediately, to prevent showing invalid search error message.
* __Gene Expression Clustering__:
    * Fixed socket hangup error.
    * Improved the error message when <3 genes are submitted for gene clustering.
    * The tooltip is now displayed when clicking an expression data cell.
    * The tool will now render even when 1 or more submitted genes have no expression data for any sample, instead of showing a computation error.
* __Sequence Reads__:
    * In the table listing available BAM files, replaced deprecated sample_type variable with tissue_type and tumor_descriptor.
    * In the initial search input, the number of available BAM files now maxes out at first 1000 files for applicable cohorts.
* Fixed issue causing an application error when searching for redacted Entity UUIDs in __Quick Search__. <!--PEAR-2032-->
* Selected values that do not match the search criteria will no longer be displayed amongst the search results in filter cards. <!--PEAR-1848-->
* Implemented fix to ensure search bar filters are included when creating or modifying gene and mutation sets. <!--PEAR-2050-->
* Fixed an issue where adding custom number range cards in the __Repository__ resulted in an infinite spinner on the cards. <!--PEAR-2114-->
* The Unexpected Error modal issue caused by rapidly clicking options in the __Customize Columns__ feature has been resolved for all tables. <!--PEAR-1999-->
* Fixed the issue where the reset button tooltip in the Customize Columns feature appeared behind other elements. <!--PEAR-2056-->
* Fixed an issue where clicking on an operator for number range filters in the __Query Expressions__ section did not remove the expected operand. <!--PEAR-2108/2109-->
* The CNV counts in the __Cancer Distribution__ table are now consistent with the associated counts in Mutation Frequency when the gene summary page is loaded from Mutation Frequency. <!--PEAR-1649-->
* Removed "undefined" text in the __Survival Plot__ of __Mutation Frequency__ when no data is available in the Mutations tab, ensuring consistent messaging with the Genes tab. <!--PEAR-1524-->
* Fixed an issue in __Set Operations__ where cohorts were not properly displayed after saving an "Unsaved_Cohort" during comparison, ensuring correct cohort selection and comparison behavior. <!--PEAR-2009-->
* Users can now change the number of rows displayed in the table when selecting an existing cohort as the basis of a new cohort. <!--PEAR-2198-->
* Fixed sorting functionality for "Submitted Gene Identifier" columns (Symbol, Ensembl ID, Entrez ID) in __Manage Sets__ to correctly sort by numbers/alphabet. <!--PEAR-2045-->
* In the __Case Summary Page__, the issue with the Create Cohort button label case count in the Most Frequent Somatic Mutations table not matching the actual cohort case count has been resolved to ensure accurate numerator and denominator calculations. <!--PEAR-2142-->
* Fixed issue where y-axis labels on __Clinical Data Analysis__ histograms were cut off; labels are now fully visible. <!--PEAR-2055-->
* The Survival Analysis section in __Cohort Comparison__ now correctly displays the message “No Survival data available for this Cohort Comparison” when there is insufficient data for the survival plot. <!--PEAR-1941-->
* Spinners on the __Cart__ page now display only for the specific download option selected, and no spinners will appear if the download does not start or has completed, ensuring consistent behavior. <!--PEAR-1722-->
* Minor text and styling fixes. <!--PEAR-2034/2122/2099/2056/2023/2077/2095/867/1908-->

### Known Issues and Workarounds
* __Section 508 Accessibility__:
    * There are known Section 508 accessibility issues that the GDC plans to address in subsequent releases. If a user encounters a Section 508 barrier, please contact GDC Support (support@nci-gdc.datacommons.io) for assistance. Known Section 508 issues are identified below.
        * There are keyboard focus and navigation issues in analysis tools that use popup windows/overlays for custom user selections. Impacted analysis tools include BAM Slicing, Sequence Reads, Gene Expression Clustering, OncoMatrix, and ProteinPaint.
        * Heatmaps within the Sequence Reads tool do not contain concise alternative text or equivalent alternatives. 
        * In the Gene Expression Clustering tool and OncoMatrix, there are no headers for genes, clusters, and/or cases in the heatmap.
        * In the Gene Expression Clustering tool, color is used to convey gene expression values but there are no patterns to convey the same information as color. Color is also used in ProteinPaint and the Sequence Reads tool to convey consequence type but there are no distinguishing patterns.
        * Some text can be difficult to read on a small screen at a 200% zoom level.
* __Survival Plot__:
    * In Mutation Frequency, the downloaded image may display a survival curve when none is plotted within the portal. <!--SV-2356-->
    * When the survival plot is zoomed in and an image is downloaded, the curves within the image may extend beyond the y-axis. <!--SV-2348-->
* __Gene Expression Clustering__:
    * The tool allows deleting the gene expression group and displays an uninformative error message after submitting the deletion.
* In __ProteinPaint__, the "Gene Expression" option is non-functional when filtering samples in a sub-track.
* Adding family_histories.relationship_age_at_diagnosis or follow_ups.other_clinical_attributes.undescended_testis_corrected_age_range as a custom filter in the __Cohort Builder__ will result in endless spinners being displayed on these filter cards. To remove these cards, close the browser tab and return to the portal. <!--SV-2532-->
* Using multiple browser tabs with the portal when adding or removing files from the __Cart__ may result in the Cart not being updated as expected. <!--SV-2412-->
* In the __files, cases, and annotations tables__, the case ID search field is case-sensitive. If the search does not return the expected results, try changing the input to uppercase as case IDs are most commonly uppercased.
* __Cohorts__ filtered by mutated genes and SSMs not in those genes will result in 0 cases since the mutations have to belong to those particular genes in order to match cases for the results. As a workaround, first filter the cohort by the mutated genes and export the cohort using the Export Cohort feature in the Cohort Bar. Then, reimport the cohort using the Import New Cohort feature before applying the SSM filters. <!--SV-2331/PEAR-1616-->
* The __Slide Image Viewer__ will display a black image temporarily if a user zooms in on a slide then switches to another slide. <!--SV-2370-->
* The TSV of the __Most Frequent Somatic Mutations table in the case summary page__ does not reflect the displayed information in the table if a search filter has been applied. <!--PEAR-2143-->
* Repeated and consecutive uses of the browser's back and/or forward buttons to return to a previously viewed page may result in a different page being displayed than the one indicated in the browser address bar. <!--SV-2552-->

## Release 2.2.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  June 26, 2024

### New Features and Changes

* __GDC 1.0__:
    * GDC 1.0 has been officially retired and can no longer be reached.
* __Annotations__:
    * Annotations tables have been added to the __project, case, and file summary pages__. <!--PEAR-1853/1839/1855-->
    * Links to the case summary page have been removed for annotations that concern a redaction of the case. <!--PEAR-1942-->
* __ProteinPaint__:
    * A new option to toggle lollipops pointing up or down is now available.
* __OncoMatrix__:
    * Advanced sorting options for power users have been added.
    * Implemented a prototype for adding gene expression rows.
* __Gene Expression Clustering__:
    * Allows re-sort cases by dictionary variable, gene mutation, or expression level.
* __Sequence Reads__:
    * Adds ability to visualize truncated BAM slice when the slice file size exceeds 20MB and streaming is terminated.
* __Cohorts__ created from analysis tools now consistently consist of a specific list of cases that will remain unchanged after a data release. This includes cohorts created from the gene and mutation summary pages. <!--PEAR-1911-->
* The files tables in the __Cart__ and the __Repository__ now allow searching for files based on the associated cases' submitter ID and UUID. <!--PEAR-1900-->
* In the __case summary page__, values whose units are days, e.g. Days to Death or Days to Birth, are now displayed in years and days as appropriate for the user's convenience. <!--PEAR-1467-->
* __Quick Search__ results and the headers of all __summary pages__ have been updated with new designs and icons. <!--PEAR-1301-->
* Word wrapping has been improved for __Quick Search__ results to avoid unexpected word breaks. <!--PEAR-1897-->
* The Best Overall Response card in the Treatments category of the __Cohort Builder__ has been moved to a new position in the category. <!--PEAR-1988-->
* The text referencing the deletion of custom sets in the __Manage Sets__ page has been updated. <!--PEAR-1939-->
* The Access column has been added to the Source Files and Download Analyses Files tables in the __file summary page__. <!--PEAR-1872-->
* Text within the downloaded histogram image has been updated for greater clarity in the __Clinical Data Analysis__ tool. <!--PEAR-1805-->
* Filter panels in the __Projects__, __Mutation Frequency__, and __Repository__ tools have been standardized and now consistently allow scrolling to occur independently of the tables on the right. <!--PEAR-1579-->

### Bugs Fixed Since Last Release
* __Section 508 Accessibility__:
    * Aria labels have been added to the tables in Set Operations. <!--PEAR-1936-->
    * The Venn diagrams in __Set Operations__ and __Cohort Comparison__ now have the appropriate alt text and roles. <!--PEAR-1936-->
    * The Venn diagram button in __Cohort Comparison__ has been updated with both an aria label and an informative label. <!--PEAR-1936-->
    * The statistics table and its TSV for Box and QQ plots in __Clinical Data Analysis__ now contain data for Q1 and Q3. <!--PEAR-1935-->
    * Alt text has been added to both the Box plot and the QQ plot in __Clinical Data Analysis__. <!--PEAR-1935-->
    * Responsiveness for the header, footer, home page, and also the Projects and Repository tools has been improved so that these areas are accessible at a 200% zoom level. <!--PEAR-1925/1912-->
* __Cases Table__:
    * The downloaded TSV now contains the expected tabs. <!--PEAR-1947-->
    * The correct number of annotations will now be displayed for each case. <!--PEAR-1965-->
    * The Customize Columns options are no longer cut off at the bottom. <!--PEAR-1997-->
    * Search now correctly displays results even if the same search input is removed and then reapplied quickly. <!--PEAR-1740-->
* __Cohorts__:
    * Cohorts containing FM-AD cases will now update correctly when users with dbGaP access to FM-AD (phs001179) log in or out. <!--SV-2389-->
    * Cohorts created based on CNV losses or gains will now have the correct composition when filtered by additional mutated genes. <!--PEAR-1597-->
* __Cohort Builder__:
    * The buckets for Age At Index will no longer display incorrect ranges and counts. <!--PEAR-1964-->
    * Cards now display at the correct width when either the browser window is small or the zoom level is increased. <!--PEAR-1945-->
* __Mutation Frequency__:
    * Gene/mutation sets created from the tables in the Mutation Frequency tool will now contain the expected genes/mutations even if the cohort has Available Data filters or Biospecimen filters. <!--SV-2314-->
    * Users will no longer be able to create several cohorts in quick succession from Mutation Frequency without waiting for previous actions to be completed. <!--PEAR-1922-->
* __OncoMatrix__:
    * Made OncoMatrix react to divide-by term edits from the label click menu.
    * Fixed the continuous term scale for density plots.
    * Gene expression variable may not show expression data for all applicable cases, especially with a large cohort size.
* __Gene Expression Clustering__:
    * Fixed the continuous term scale for density plots.
* __Sequence Reads__:
    * The tool now displays the correct number of available BAM files when a cohort filter is in use.
* __Cohort MAF__:
    * Added the "tumor_bam_uuid" column.
* Limited the CSS reset to avoid conflict with embedded styles, by using scoped normalize CSS rules.
* Fixed conflicting CSS that can alter portal styling.
* Fixed an issue where __Sample Sheet__ downloads can be incomplete due to missing sample type information. <!--PEAR-1972-->
* Addressed an issue where changes to default filters in the __Cohort Builder__ and the __Repository__ may not be reflected after a release. <!--PEAR-1960-->
* The __Query Expressions__ section now correctly displays a maximum of 3 rows by default. Additionally, the button to display more than 3 rows at a time is enabled only when the cohort query exceeds 3 rows. <!--PEAR-1578-->
* The loading spinner is no longer displayed above the other areas of the Analysis Center when the __Cohort Comparison__ tool is loading. <!--SV-2360-->
* Default values in the Custom Bins modal within the __Clinical Data Analysis__ tool are now properly updated when the user toggles between displaying continuous values in days and years. <!--PEAR-1981-->
* The right side of the chart on the __Home Page__ is no longer cut off at smaller browser sizes. <!--PEAR-1167-->
* Tooltips are no longer displayed when there is no description available for filter properties. <!--SV-2425/PEAR-869-->

### Known Issues and Workarounds
* __Section 508 Accessibility__:
    * There are known Section 508 accessibility issues that the GDC plans to address in subsequent releases. If a user encounters a Section 508 barrier, please contact GDC Support (support@nci-gdc.datacommons.io) for assistance. Known Section 508 issues are identified below.
        * There are keyboard focus and navigation issues in analysis tools that use popup windows/overlays for custom user selections. Impacted analysis tools include BAM Slicing, Sequence Reads, Gene Expression Clustering, OncoMatrix, and ProteinPaint.
        * Heatmaps within the Sequence Reads tool do not contain concise alternative text or equivalent alternatives. Additionally, an equivalent alternative to the body plot on the home page is not available.
        * In the Gene Expression Clustering tool and OncoMatrix, there are no headers for genes, clusters, and/or cases in the heatmap.
        * In the Gene Expression Clustering tool, color is used to convey gene expression values but there are no patterns to convey the same information as color. Color is also used in ProteinPaint and the Sequence Reads tool to convey consequence type but there are no distinguishing patterns.
        * Some text can be difficult to read on a small screen at a 200% zoom level.
* __Survival Plot__:
    * The survival plot in Cohort Comparison does not display text indicating that there is insufficient survival data to plot. <!--SV-2357-->
    * In Mutation Frequency, the downloaded image may display a survival curve when none is plotted within the portal. <!--SV-2356-->
    * When the survival plot is zoomed in and an image is downloaded, the curves within the image may extend beyond the y-axis. <!--SV-2348-->
* __Cart__:
    * Spinners on the Download Cart and Download Associated Data buttons may be displayed longer than expected. This is a visual issue and does not affect the use of these buttons. <!--SV-2343-->
    * Using multiple browser tabs with the portal when adding or removing files from the cart may result in the cart not being updated as expected. <!--SV-2412-->
* In the __files, cases, and annotations tables__, the case ID search field is case-sensitive. If the search does not return the expected results, try changing the input to uppercase as case IDs are most commonly uppercased.
* __Cohorts__ filtered by mutated genes and SSMs not in those genes will result in 0 cases since the mutations have to belong to those particular genes in order to match cases for the results. As a workaround, first filter the cohort by the mutated genes and export the cohort using the Export Cohort feature in the Cohort Bar. Then, reimport the cohort using the Import New Cohort feature before applying the SSM filters. <!--SV-2331/PEAR-1616-->
* The __Slide Image Viewer__ will display a black image temporarily if a user zooms in on a slide then switches to another slide. <!--SV-2370-->
* The annotations table in the __file summary page__ does not include the Case ID column. This column is planned to be added in a future update. <!--PEAR-1899-->
* In __ProteinPaint__, the "Gene Expression" option is non-functional when filtering samples in a sub-track.
* In __Gene Expression Clustering__, the tooltip is not displayed when clicking an expression data cell.
* The custom range inputs for the __Age at Index__ card in the __Cohort Builder__ are not behaving as expected. As a workaround, use the predefined ranges available. Alternatively, use the custom range inputs on the Days tab to query for ages in years.

## Release 2.1.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  April 30, 2024

### New Features and Changes
* __Annotations Browser__
    * The Annotations Browser and annotation summary page have been implemented. <!--PEAR-1838/1181/1182/1160-->
* __Repository__:
    * Files can now be filtered by Tissue Type, Tumor Descriptor, Specimen Type, and Preservation Method. <!--PEAR-1514-->
    * Metadata and Sample Sheet downloads have been added to the Repository. <!--PEAR-1863-->
* __Cohort Level MAF__:
    * A cohort level MAF analysis tool has been added to the Analysis Center.
* __BAM Slicing Download and Sequence Reads__:
    * In BAM Slicing Download, call GDC API directly from client without going through ProteinPaint backend. No limits are applied on slicing region size or BAM slice file size.
    * In Sequence Reads Visualization, user can slice a BAM with a range lower than 300Kb, and if the resulting BAM slice is under 100Mb. Slicing and caching a BAM slice bigger than 100Mb will abort and user will be notified to reduce region size and try again. Before creating new cache file, find out old enough ones to delete to free up storage.
    * For both apps: The table listing available cases and bam files can be filtered by assay types.
* __Gene Expression Clustering__:
    * Enable gene variant legend group filter.
    * Support creating a single-case cohort.
    * Supported more clustering and distance calculation methods.
* __OncoMatrix__:
    * Enable downloading data.
    * Hide synonymous mutations by default.
    * Improve the matrix sorting options to easily toggle sorting by cnv and/or consequence.
    * Add Mutation and CNV control buttons, and hide CNV by default.
    * Create a mutations/consequences legend group for mutations.
    * Enable selecting individual mutation classes upon clicking the Mutation/CNV button.
    * Support creating a single-case cohort.
    * Display hints about persisted matrix gene set and option to unhide CNV and mutations when there is no matrix data to render.
    * Group similar mutation class colors together when sorting matrix samples and if CNVs are displayed.
    * Add "Single" style to render consequence data, as alternative to Stacked and OncoPrint styles.
* __ProteinPaint__:
    * Allow visualizing SSM in any genomic locus, besides "protein" mode.
    * Support creating a single-case cohort.
* The performance of the __Clinical Data Analysis__ tool has been improved, especially when large cohorts are used with QQ plots. <!--PEAR-1536-->
* __Quick Search__ now returns results for the latest versions of files when searching for older versions of those files. <!--PEAR-1804-->
* The X button on the __Unexpected Error__ dialog box has been removed. <!--SV-2367-->
* Buttons for launching demos have been removed from the selection view of __Cohort Comparison__ and __Set Operations__. <!--SV-2328/2327-->
* Responsiveness improvements have been made to the __Analysis Center__ and the __Cohort Bar__. <!--PEAR-1836-->
* The UX/UI for the __Cohort Builder__ has been improved. <!--PEAR-1547-->
* The __case summary page__ has been enhanced with a table listing all the files associated with the case. Additionally, a link to the table is now available in the header of the summary page, and information has been added to the File Counts summary tables to lead users to the new files table. The clinical and biospecimen supplements tables have also been removed from the case summary page. <!--PEAR-1822/1849/1833/PEAR-1832-->
* Set names for sets of the same type are now enforced to be unique when editing names in __Manage Sets__. <!--PEAR-1359-->
* Number range cards in the __Cohort Builder__ no longer display the custom range option when there is no data. <!--PEAR-661-->
* The Cohort Buider image on the __home page__ has been updated to reflect the latest design. <!--PEAR-1593-->
* The tooltip on the __Mutation Frequency__ card in the Analysis Center has been updated. <!--PEAR-1877-->

### Bugs Fixed Since Last Release
* __Section 508 Accessibility__:
    * Small aria-label inconsistencies have been addressed. <!--PEAR-1715-->
    * Keyboard focus is now returned to the triggering element when modals are closed. <!--PEAR-1658-->
    * Screen readers will now read out the contents of toast messages. <!--PEAR-1765-->
    * Toggles in the Clinical Data Analysis tool now have the correct number of labels. <!--PEAR-1764-->
    * Table header checkboxes are now correctly labelled. <!--PEAR-1754-->
    * Modal icons now have appropriate null alt text. <!--PEAR-1753-->
    * Assistive technologies no longer behave incorrectly with some controls due to incorrect, missing, or redundant labels, attributes, or roles. <!--PEAR-1672-->
    * Aria labels have been added to Cancer Gene Census annotation icon in Mutation Frequency. <!--PEAR-740-->
    * The Survival icon is now appropriately hidden from the accessibility tree for the benefit of screen readers. <!--PEAR-739-->
* __Cohorts__:
    * Using "Save As" to replace a cohort with itself will no longer result in an error notification despite the replacement being successful. <!--SV-2363-->
    * Saving a cohort that was previously saved now displays the correct message. <!--PEAR-1651-->
    * Cohorts will now display data in Mutation Frequency, Cohort Builder, and the summary charts even when removing gene/mutation filters from a cohort temporarily results in 0 cases. <!--SV-2414-->
    * Cohorts now contain the correct cases when created from the cases table by using the "Existing Cohort With Selected Cases" and "Existing Cohort Without Selected Cases" options with a cohort containing gene or mutation filters. <!--PEAR-1043-->
    * When saving a cohort, the confirmation notification will no longer be automatically dismissed before the saving dialog has closed. <!--SV-2366-->
* __Cohort Builder__:
    * Cohort Builder cards for number ranges now display an informative message rather than a spinner when there is no data for the facet. <!--PEAR-1646-->
    * Removing a custom Cohort Builder card no longer incorrectly removes the associated filters from the current cohort. <!--PEAR-1612-->    
    * Filters related to numeric values in the Cohort Builder now correctly displays the numbers entered. <!--SV-2383-->
* __Case Summary Page__:
    * The Biospecimen tree in the case summary page is no longer hidden when the bioId provided in the URL does not exist. <!--PEAR-1202-->
    * The error that sometimes occurs when viewing the __Follow-Ups/Molecular Tests__ tab in the case summary page has been resolved. <!--SV-2431-->
* __Mutation Frequency__:
    * The survival plot in __Mutation Frequency__ no longer flickers when the cohort has 0 cases. <!--SV-2331/PEAR-1701-->
    * Attempting to download a TSV of all the mutations in the GDC no longer results in an error due to the length of time needed to generate the TSV. <!--SV-2388-->
* __All ProteinPaint-based Tools__:
    * In GDC query, do not supply empty "case_filters{content[]}" that will slow down API. Lollipop and OncoMatrix are now faster when there's no cohort.
    * Updated mutation class definitions and rank for protein_altering_variant. Affects all tools that can show mutation data.  
    * Deprecated term "sample_type" is dropped from GDC dictionary.
* __BAM Slicing Download and Sequence Reads__:
    * When downloading GDC BAM slice (no caching), do not limit request region max size.
    * Reloading page while streaming/downloading GDC BAM slice to client will not crash server.
    * App UI requires hitting Enter to search by GDC file or case, and will no longer auto search (on pressing any key) to avoid showing duplicate SSM table.
    * BAM track bug fix to handle reads with no sequence.
    * BAM track bug fix for hide/show toggling at track menu.
* __Disco Plot__:
    * Bug fix for disco plot launched from sunburst showing AAchange in sandbox header rather than undefined.
    * Pass the cohort filter to the lollipop track from the matrix and disco plot label click.
* __Gene Expression Clustering__:
    * Enable gene variant legend group filter.
    * Support creating a single-case cohort.
    * Supported more clustering and distance calculation methods.
* __OncoMatrix__:
    * Fix position errors after OncoMatrix/hierCluster zooming in/out caused by outdated imgBox.
    * Do not allow hiding all the alteration groups.
    * Disable the geneset submit button when there there is less than a minNumGenes option (3 for hier cluster, 1 for matrix).
    * Add to OncoMatrix mutation/cnv buttons all available mutation/cnv classes in all GDC instead of only within current cohort.
    * Change the definition of truncating/protein-changing mutation, change OncoMatrix mutation classes sorting order.
    * Fix the detection of sorting-related updates in the matrix app, as distinct from the Gene Expression Clustering.
    * Pass the cohort filter to the lollipop track from the matrix and disco plot label click.
* __ProteinPaint__:
    * Sample summary table will scroll if too tall.
    * Bug fix to convert "case." to "cases." in case_filters[] for GDC mds3 sunburst clicking to load sample table.
    * Do not force the sample table to be positioned relative to screen bottom after a sunburst click.
    * Prevent double-clicking on a sunburst ring so that same sample table will not appear duplicated.
    * Bug fix for Lollipop category total sample count to respond/shrink with cohort change.
* Tokens are no longer refreshed when the __User Profile__ is viewed. <!--PEAR-1818-->
* __Quick Search__ now correctly displays results even if the same search input is applied twice quickly. <!--SV-2410-->
* In __Set Operations__, saving gene and mutation sets will now be successful if the saving dialog is manually dismissed after the Save button is clicked. <!--SV-2368-->
* Users will no longer be able to download more than 5 GB of files in total at a time via the browser from the __cart__. <!--SV-2342-->
* Table buttons in __Clinical Data Analysis__ no longer overlay the survival plot on smaller screens when many survival plots are displayed at the same time. <!--PEAR-1600-->
* The correct file size total will now be displayed in the __Repository__ when filtering is applied within the tool and the active cohort contains Available Data filters. <!--SV-2376-->
* Downloading the **Clinical/Biospecimen TSV or JSON** before the cohort has fully loaded will no longer result in an error. <!--SV-2402-->

### Known Issues and Workarounds

* __Section 508 Accessibility__:
    * There are known Section 508 accessibility issues that the GDC plans to address in subsequent releases. If a user encounters a Section 508 barrier, please contact GDC Support (support@nci-gdc.datacommons.io) for assistance. Known Section 508 issues are identified below.
        * There are keyboard focus and navigation issues in analysis tools that use popup windows/overlays for custom user selections. Impacted analysis tools include BAM Slicing, Sequence Reads, Gene Expression Clustering, OncoMatrix, and ProteinPaint.
        * Heatmaps within the Sequence Reads tool do not contain concise alternative text or equivalent alternatives. Additionally, equivalent alternatives to the Box plots, QQ plots, Venn diagrams, and the body plot on the home page are not available.
        * In the Gene Expression Clustering tool and OncoMatrix, there are no headers for genes, clusters, and/or cases in the heatmap.
        * In the Gene Expression Clustering tool, color is used to convey gene expression values but there are no patterns to convey the same information as color. Color is also used in ProteinPaint and the Sequence Reads tool to convey consequence type but there are no distinguishing patterns.
        * Some text can be difficult to read on a small screen at a 200% zoom level.
* __Cohorts__:
    * Cohorts are under active development and their behavior may change in the first several months after the release of GDC Portal 2.0. As this process may result in the loss of saved cohorts on the portal, we highly recommend [exporting cohorts](/Data_Portal/Users_Guide/getting_started.md#main-toolbar) locally.
    * Cohorts created based on CNV losses or gains may not have the correct composition when filtered by additional mutated genes. As a workaround, first filter by the mutated genes before creating cohorts based on CNV losses and gains.
    * Cohorts filtered by mutated genes and SSMs not in those genes may unexpectedly result in 0 cases. <!--SV-2331/PEAR-1616-->
    * Cohorts containing FM-AD cases may not update correctly when users with dbGaP access to FM-AD (phs001179) log in or out. As a workaround, logging in before creating cohorts with FM-AD cases is recommended. <!--SV-2389-->
* __Survival Plot__:
    * The survival plot in Cohort Comparison does not display text indicating that there is insufficient survival data to plot. <!--SV-2357-->
    * In Mutation Frequency, the downloaded image may display a survival curve when none is plotted within the portal. <!--SV-2356-->
    * When the survival plot is zoomed in and an image is downloaded, the curves within the image may extend beyond the y-axis. <!--SV-2348-->
* __Cart__:
    * Spinners on the Download Cart and Download Associated Data buttons may be displayed longer than expected. This is a visual issue and does not affect the use of these buttons. <!--SV-2343-->
    * Using multiple browser tabs with the portal when adding or removing files from the cart may result in the cart not being updated as expected. <!--SV-2412-->
* The aggregated MAF generated using the __Cohort Level MAF__ tool is missing the tumor_bam_uuid column. The tumor_sample_uuid and case_id should be used for reproducibility until the tumor_bam_uuid has been added.
* In both the __Sequence Reads__ and __BAM Slicing Download__ tools, the number of available BAM files may be overcounted when a cohort filter is in use.
* Gene/mutation sets created from the tables in the __Mutation Frequency__ tool may contain 0 genes/mutations if the cohort has Available Data filters or Biospecimen filters. <!--SV-2314-->
* The TSV of the __cases table__ may not contain the expected tabs. <!--DEV-2324-->
* In the Repository and cases table, the case ID search field is case-sensitive. If the search does not return the expected results, try changing the input to uppercase as case IDs are most commonly uppercased.
* When the __Cohort Comparison__ tool is loading, the loading spinner may be displayed above the other areas of the Analysis Center. <!--SV-2360-->
* The __Slide Image Viewer__ will display a black image temporarily if a user zooms in on a slide then switches to another slide. <!--SV-2370-->

## Release 2.0.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  February 8, 2024

### New Features and Changes

GDC 2.0 is a major update to the original GDC Data Portal introduced in 2016. This latest version adopts a "cohort-centric" workflow, in which users build custom sets of cases to analyze, and introduces several new analysis tools. New features of GDC 2.0 include:

* A cohort-centric workflow in which a cohort is first built and then analyzed using tools on the Data Portal. All of these functionalities can be reached from the Analysis Center.  
    * This includes a toolbar, that can be used to view or modify an existing cohort while using any analysis tool.
* Core tools that compose the main functionalities of the GDC Data Portal:
    * __Cohort Builder:__ Build a cohort of cases using clinical and biospecimen properties
    * __Repository:__ Download files based on a specific cohort
    * __Projects:__ Browse, filter, and create cohorts based on GDC projects
* Analysis tools that analyze specific cohorts:
    * __Mutation Frequency:__ Analyze somatic mutations that were called in the WXS and Targeted Sequencing pipelines and their associated genes
    * __Clinical Data Analysis:__ Analyze and visualize clinical data associated with your cohort
    * __Cohort Comparison:__ Analyze the properties of multiple cohorts
    * __Set Operations:__ Display a Venn diagram and compare/contrast cohorts or gene/mutation sets
    * __BAM Slicing Download:__ Download a specific region of a BAM file created by the GDC
    * __ProteinPaint:__ Visualize somatic mutations on a specific linear gene or chromosomal region
    * __Gene Expression Clustering:__ Visualizes gene expression clustering for a specific cohort
    * __Sequence Reads:__ Visualize the reads within a specific BAM file
    * __OncoMatrix:__ Visualize the most commonly mutated genes across a cohort

### Bugs Fixed Since Last Release

Not applicable as this is the initial release of GDC 2.0.

### Known Issues and Workarounds

* __Section 508 Accessibility__:
    * There are known Section 508 accessibility issues that the GDC plans to address in subsequent releases. If a user encounters a Section 508 barrier, please contact GDC Support (support@nci-gdc.datacommons.io) for assistance. Known Section 508 issues are identified below.
        * There are keyboard focus and navigation issues in analysis tools that use popup windows/overlays for custom user selections. Impacted analysis tools include BAM Slicing, Sequence Reads, Gene Expression Clustering, OncoMatrix, and ProteinPaint.
        * Heatmaps within the Sequence Reads tool do not contain concise alternative text or equivalent alternatives. Additionally, equivalent alternatives to the Box plots, QQ plots, Venn diagrams, and the body plot on the home page are not available.
        * In the Gene Expression Clustering tool and OncoMatrix, there are no headers for genes, clusters, and/or cases in the heatmap.
        * In the Gene Expression Clustering tool, color is used to convey gene expression values but there are no patterns to convey the same information as color. Color is also used in ProteinPaint and the Sequence Reads tool to convey consequence type but there are no distinguishing patterns.
        * Some text can be difficult to read on a small screen at a 200% zoom level.
        * Keyboard focus is not returned to the triggering element when modals are closed.
        * Assistive technologies may not behave correctly with some controls due to incorrect, missing, or redundant labels, attributes, or roles.
* __Cohorts__:
    * Cohorts are under active development and their behavior may change in the first several months after the release of GDC Portal 2.0. As this process may result in the loss of saved cohorts on the portal, we highly recommend [exporting cohorts](/Data_Portal/Users_Guide/getting_started.md#main-toolbar) locally.
    * Cohorts created based on CNV losses or gains may not have the correct composition when filtered by additional mutated genes. As a workaround, first filter by the mutated genes before creating cohorts based on CNV losses and gains.
    * Cohorts filtered by mutated genes and SSMs not in those genes may unexpectedly result in 0 cases. <!--SV-2331/PEAR-1616-->
    * When saving a cohort, the confirmation notification may be automatically dismissed before the saving dialog has closed. <!--SV-2366-->
    * Using "Save As" to replace a cohort with itself will result in an error notification despite the replacement being successful. <!--SV-2363-->
    * Cohorts containing FM-AD cases may not update correctly when users with dbGaP access to FM-AD (phs001179) log in or out. As a workaround, logging in before creating cohorts with FM-AD cases is recommended. <!--SV-2389-->
    * If removing gene/mutation filters from a cohort temporarily results in 0 cases, cohorts may not display data in Mutation Frequency, Cohort Builder, and the summary charts. As a workaround, remove the gene and mutation filters, then add them back. <!--SV-2414-->
* __Survival Plot__:
    * The survival plot in Cohort Comparison does not display text indicating that there is insufficient survival data to plot. <!--SV-2357-->
    * The survival plot in Mutation Frequency may flicker when the cohort has 0 cases. <!--SV-2331/PEAR-1701-->
    * In Mutation Frequency, the downloaded image may display a survival curve when none is plotted within the portal. <!--SV-2356-->
    * When the survival plot is zoomed in and an image is downloaded, the curves within the image may extend beyond the y-axis. <!--SV-2348-->
* __Cart__:
    * Spinners on the Download Cart and Download Associated Data buttons may be displayed longer than expected. This is a visual issue and does not affect the use of these buttons. <!--SV-2343-->
    * More than 5 GB of files in total may be downloaded at a time via the browser if the user first attempts to download controlled access data without being logged in, then logs in via the information dialog displayed before continuing with the download. <!--SV-2342-->
    * Using multiple browser tabs with the portal when adding or removing files from the cart may result in the cart not being updated as expected. <!--SV-2412-->
* __Mutation Frequency__:
    * Gene/mutation sets created from the tables in Mutation Frequency may contain 0 genes/mutations if the cohort has Available Data filters or Biospecimen filters. <!--SV-2314-->
    * Attempting to download a TSV of all the mutations in the GDC may result in an error due to the length of time needed to generate the TSV. As a workaround, limit the number of mutations downloaded to 1.5 million. <!--SV-2388-->
* __Main Toolbar__:
    * Attempting to download the **Clinical/Biospecimen TSV or JSON** before the cohort has fully loaded may result in an error. <!--SV-2402-->
    * The TSV of the **cases table** may not contain the expected tabs. <!--DEV-2324-->
* __OncoMatrix__:
    * Manually deleting all genes will result in an error message "Error: Cannot read properties of undefined (reading 'lst')". The user can close and re-open OncoMatrix for use.
    * Dragging genes only works once. After one gene is dragged to a new position, no genes can be dragged to new positions.
    * In OncoMatrix, the cases may not get re-sorted as expected after a certain sequence of actions
* __ProteinPaint__:
    * A nested filter may be constructed for a Lollipop subtrack, e.g. sex=male AND ( primarysite=aa OR disease=bb ), but cannot be translated into GDC cohort filters. The translation code has a preliminary implementation that only works for "flat" filters without nesting.
    * Cohorts cannot be created using the Create Cohort button in ProteinPaint for a single sample
    * In ProteinPaint, the total number of samples in a category breakdown and the total number of samples in the sunburst ring are not based on a user's current cohort
    * In ProteinPaint, when clicking on a sunburst ring, the sample table is not showing up
    * In ProteinPaint, a Disco plot launched from the sunburst ring can show "undefined" in the plot header
    * A ProteinPaint plot launched from OncoMatrix and Gene Expression Clustering does not observe the current cohort and displays mutated cases for all GDC
* In the __Gene Expression Clustering__ tool, if any part of the dendrogram is selected and the current cohort is modified, then the new dendrogram will render with scattered subtrees selected.
* The "A" in the Allele Summary text is cut off in the __Sequence Reads__ tool.
* __Quick Search__ may not display results if the same search input is applied twice quickly. As a workaround, temporarily change the input before reentering the intended search. <!--SV-2410-->
* Filters related to numeric values may display a smaller number than what the user entered within the __Cohort Builder__. This is a visual issue and does not affect the filters applied to the cohort. <!--SV-2383-->
* When the __Cohort Comparison__ tool is loading, the loading spinner may be displayed above the other areas of the Analysis Center. <!--SV-2360-->
* The __Repository__ tool may display an incorrect file size total of 0 bytes when filtering is applied within the tool and the active cohort contains Available Data filters. <!--SV-2376-->
* The __Slide Image Viewer__ will display a black image temporarily if a user zooms in on a slide then switches to another slide. <!--SV-2370-->
* In __Set Operations__, the saving of gene and mutation sets may be unsuccessful if the saving dialog is manually dismissed after the Save button is clicked. <!--SV-2368-->
* Clicking the X button on the __Unexpected Error dialog box__ does not dismiss it. The workaround is to click the OK button. <!--SV-2367-->

## Release 1.30.4

* __GDC Product__: GDC Data Portal
* __Release Date__:  May 11, 2023

### New Features and Changes

* The GDC Legacy Archive has officially been retired.
    * The Legacy Archive Portal can no longer be reached.
    * Any API call to query files from the Legacy Archive will no longer work.
    * Downloads for files from the Legacy Archive will work normally with manifests that were generated previously.


### Bugs Fixed Since Last Release

* The Clinical TSV download in the case entity page and cart now contain TSVs for pathology detail, follow up, and molecular test entities.   <!--SV-2004-->
* Fixed bug in which demographic information would not download as a TSV when diagnosis information was not available. <!--DEV-1238-->

### Known Issues and Workarounds

*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  Negative numbers may be displayed for the Missing value category in the Treatment node within a Clinical Analysis.  This occurs with projects that have multiple treatment nodes per case. All other values should be accurate. <!--SV-1604-->
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->

## Release 1.30.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  July 8, 2022

### New Features and Changes

* None

### Bugs Fixed Since Last Release

* Fixed error in which the manifest could not be downloaded directly from the GDC Data Portal in certain instances. <!--SV-2112-->

### Known Issues and Workarounds

*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  Negative numbers may be displayed for the Missing value category in the Treatment node within a Clinical Analysis.  This occurs with projects that have multiple treatment nodes per case. All other values should be accurate. <!--SV-1604-->
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    * The footer says version 1.9, but it is actually 1.13
    *	Filtering by vital_status does not function in the Legacy Archive due to updates in how this property has been indexed.  A workaround is to perform the case level filtering in the GDC Data Portal and copy the filter string for use in the Legacy Archive or the legacy API. <!--SV-1508-->
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->


## Release 1.29.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  August 23, 2021

### New Features and Changes <!--REQ-431-->

* None

### Bugs Fixed Since Last Release

* Fixed error in which the data summaries in the clinical analysis and cart pages were only partially displayed when viewed in Chrome v.91.0.4472. <!--SV-1940-->

### Known Issues and Workarounds

*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  Negative numbers may be displayed for the Missing value category in the Treatment node within a Clinical Analysis.  This occurs with projects that have multiple treatment nodes per case. All other values should be accurate. <!--SV-1604-->
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    * The footer says version 1.9, but it is actually 1.13
    *	Filtering by vital_status does not function in the Legacy Archive due to updates in how this property has been indexed.  A workaround is to perform the case level filtering in the GDC Data Portal and copy the filter string for use in the Legacy Archive or the legacy API. <!--SV-1508-->
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->

## Release 1.28.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  May 17, 2021

### New Features and Changes <!--REQ-422-->

* New columns were added to the "molecular test" table at the bottom of the case entity page to display additional molecular test fields.

### Bugs Fixed Since Last Release

* None

### Known Issues and Workarounds

*  When accessing the data portal with Chrome v.91.0.4472, users may experience some display errors.  This includes the data summary in the clinical analysis and cart pages being only partially displayed. <!--SV-1940-->
*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  Negative numbers may be displayed for the Missing value category in the Treatment node within a Clinical Analysis.  This occurs with projects that have multiple treatment nodes per case. All other values should be accurate. <!--SV-1604-->
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    * The footer says version 1.9, but it is actually 1.13
    *	Filtering by vital_status does not function in the Legacy Archive due to updates in how this property has been indexed.  A workaround is to perform the case level filtering in the GDC Data Portal and copy the filter string for use in the Legacy Archive or the legacy API. <!--SV-1508-->
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->

## Release 1.25.1

* __GDC Product__: GDC Data Portal
* __Release Date__:  August 14, 2020

### New Features and Changes <!--REQ-408-->

* API improvements were made to increase portal performance.

### Bugs Fixed Since Last Release

* None

### Known Issues and Workarounds

*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  Negative numbers may be displayed for the Missing value category in the Treatment node within a Clinical Analysis.  This occurs with projects that have multiple treatment nodes per case. All other values should be accurate. <!--SV-1604-->
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    * The footer says version 1.9, but it is actually 1.13
    *	Filtering by vital_status does not function in the Legacy Archive due to updates in how this property has been indexed.  A workaround is to perform the case level filtering in the GDC Data Portal and copy the filter string for use in the Legacy Archive or the legacy API. <!--SV-1508-->
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->

## Release 1.25.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  July 2, 2020

### New Features and Changes <!--REQ-395-->

* Suppressed Experimental Strategy filter on the Exploration page as this currently filters for files with a particular strategy, not for cases.  This may cause confusion amongst users.  The filter will be re-instated in a future release once the logic is available to filter more appropriately for cases tied to a specific strategy. <!--PRTL-3119-->
* Updated the filter control panel styling across the Portal to have clearer titles (e.g. "Search Cases" instead of "Cases" in the quick search box). <!--PRTL-2987-->
* Made minor updates to the styling of the filter query display at the top of the Exploration page (spacing, borders). <!--PRTL-2988-->
* Added an expand/collapse control to the quick search bar of Clinical tab on the Exploration page, to be consistent with other Exploration tabs. <!--PRTL-2990-->
* Added a clear title above the counts in each filter control panel across the Portal (e.g. "# Cases", "# Genes", etc.). <!--PRTL-2991-->
* Moved various action buttons above the results table on the Repository Page to more accessible locations. <!--PRTL-2994-->
* Improved load time of the initial custom filter list on the Repository Page, when clicking "Add a Filter Filter" or "Add a Case/Biospecimen Filter". <!--PRTL-3057-->

### Bugs Fixed Since Last Release

* Fixed a bug in the Age at Diagnosis table on the Cohort Comparison page, where the # of cases in the table was not consistent with the # of cases shown when clicking the link to the Exploration page. <!--PRTL-3032-->
* Fixed minor positional accuracy issue of the lollipop data points on the Protein Viewer. <!--PRTL-2967-->
* Fixed bug on the Protein Viewer where, if clicking to switch between different lollipop data points, details of the previous lollipop was not closing. <!--PRTL-2968-->
* Fixed bug where the quick search bar on the Exploratin Page's Genes filter tab was not expanding/collapsing properly. <!--PRTL-2989-->
* Fixed bug in the pop-up warning message when adding or removing items from the Cart, where long filenames were spilling outside the border of the pop-up. <!--PRTL-3031-->
* Fixed typo in the "View Cases in Exploration" button on the Repository page. <!--PRTL-3122-->
* Fixed typo in the pop-up user consent message when downloading controlled files from the Cart. <!--PRTL-3123-->

### Known Issues and Workarounds

*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  Negative numbers may be displayed for the Missing value category in the Treatment node within a Clinical Analysis.  This occurs with projects that have multiple treatment nodes per case. All other values should be accurate. <!--SV-1604-->
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    * The footer says version 1.9, but it is actually 1.13
    *	Filtering by vital_status does not function in the Legacy Archive due to updates in how this property has been indexed.  A workaround is to perform the case level filtering in the GDC Data Portal and copy the filter string for use in the Legacy Archive or the legacy API. <!--SV-1508-->
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->

## Release 1.24.1

* __GDC Product__: GDC Data Portal
* __Release Date__:  March 10, 2020

### New Features and Changes <!--REQ-395-->

* Removed unnecessary comma and y-axis value from title of the mutation details pop-up in the Protein Viewer. <!--PRTL-2969-->
* Added Tobacco Smoking Status field to the Exposures tab on the Case entity page. <!--PRTL-2952-->
* Added a link to the Cart where users can access instructions for downloading the GDC Genome Build reference files. <!--PRTL-2929-->
* Added logic to prevent duplicate fetching of data for Clinical Analysis survival plots and optimize rendering. <!--PRTL-2910-->
* Added a button to clear searches for certain Portal search controls that were previously missing this ability. <!--PRTL-2900-->
* Reduced whitespace between Oncogrid and its control panel to optimize spacing and layout. <!--PRTL-2899-->
* Made entire Clinical Analysis results page responsive (card columns now scale & stack in response to the size of the browser window). <!--PRTL-2882-->
* Replaced Clinical Analysis function for printing clinical cards to a single PDF file, with more flexible functionality to instead download all the cards in SVG and/or PNG format. <!--PRTL-2870-->
* Added message to notify users when they try to access the Portal using Microsoft Internet Explorer, indicating which browsers are officially supported. <!--PRTL-2868-->
* Added arrow icon to sortable columns across the Portal to indicate the current sort direction. <!--PRTL-2330-->

### Bugs Fixed Since Last Release

* Fixed bug where clicking a primary site on the Human Body Image was not re-directing to the Exploration page. <!--PRTL-2979-->
* Fixed layout issue where long Annotation Notes were exceeding the border of the text box. <!--PRTL-2944-->
* Fixed layout issue where the Repository header and action buttons were scaling and wrapping incorrectly if the browser window is shrunk beyond a certain threshold. <!--PRTL-2938-->
* Fixed layout issue where the responsive Clinical Analysis Cards were clipping improperly as the browser window is shrunk beyond a certain threshold. <!--PRTL-2928-->
* Fixed bug where the Clinical Tab on the Exploration page was crashing when entering a custom range of Years for the Age at Diagnosis facet. <!--PRTL-2913-->
* Fixed various minor cosmetic and color issues in PNG, SVG downloads of the Clinical Analysis survival plots. <!--PRTL-2887-->
* Fixed bug where the x-axis in PNG, SVG downloads of histograms across the Portal was being bolded incorrectly. <!--PRTL-2885-->
* Fixed bug where the expand/collapse symbols in the UI were incorrectly being exported in the TSV download of the Projects table. <!--PRTL-2883-->
* Fixed bug where Oncogrid's modal for customizing colors could not be scrolled below the fold if it was shrunk beyond a certain threshold. <!--PRTL-2881-->
* Fixed incorrect DTT hyperlink in the GDC Apps menu. <!--PRTL-2867-->
* Fixed bug where the "dbSNP rs ID" facet could not be minimized in the Exploration page's Mutations facet tab. <!--PRTL-2840-->
* Fixed layout issue where the Portal's header incorrectly overlaps some content when a notification banner is displayed. <!--PRTL-2753-->
* Fixed some minor layout & styling issues in the Exploration page's facets panel. <!--PRTL-2614-->
* Fixed bug where the Case ID on the Exploration page's Cases facet tab was not searchable in certain scenarios. <!--PRTL-2587-->
* Fixed bug where the Expand/Collapse button state was not changing properly when being used in the Biospecimen section of the Case entity page. <!--PRTL-2575-->
* Fixed incorrect capitalization of "dbGaP" in the Summary section of the Project entity page. <!--PRTL-2436-->
* Fixed layout issue where the Advanced Search query box on the Repository page could expanded beyond the margins of the box's border. <!--PRTL-2271-->

### Known Issues and Workarounds

*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  Negative numbers may be displayed for the Missing value category in the Treatment node within a Clinical Analysis.  This occurs with projects that have multiple treatment nodes per case. All other values should be accurate. <!--SV-1604-->
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    * The footer says version 1.9, but it is actually 1.13
    *	Filtering by vital_status does not function in the Legacy Archive due to updates in how this property has been indexed.  A workaround is to perform the case level filtering in the GDC Data Portal and copy the filter string for use in the Legacy Archive or the legacy API. <!--SV-1508-->
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->

## Release 1.23.1

* __GDC Product__: GDC Data Portal
* __Release Date__:  December 10, 2019

### New Features and Changes <!--REQ-396-->

* Updated display of x-axis units on the homepage Human Body chart to more easily display increased case counts for newly-added projects <!--PRTL-2925-->

### Bugs Fixed Since Last Release

* None

### Known Issues and Workarounds

*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  Negative numbers may be displayed for the Missing value category in the Treatment node within a Clinical Analysis.  This occurs with projects that have multiple treatment nodes per case. All other values should be accurate. <!--SV-1604-->
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    * The footer says version 1.9, but it is actually 1.13
    *	Filtering by vital_status does not function in the Legacy Archive due to updates in how this property has been indexed.  A workaround is to perform the case level filtering in the GDC Data Portal and copy the filter string for use in the Legacy Archive or the legacy API. <!--SV-1508-->
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->

## Release 1.23.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  November 6, 2019

### New Features and Changes <!--REQ-390-->

* Added Clinical Data Analysis feature that allows Users to: <!--FEAT-520, FEAT-521-->
    * Explore clinical data via the new Clinical Tab on the Exploration page.
    * Build custom Case sets based on that clinical data for later analysis.
    * Create an analysis to examine the clinical variables in a Case set, using various tools including histograms, survival plots, box plots, QQ plots, and custom binning.
    * Download the data (as TSV, JSON) and plots (as PNG, SVG) of each clinical variable in an anlysis.
    * Save an analysis to local storage to resume later (as long as storage is not cleared).
* Added links to CIViC annotations on the Gene and Mutation entity pages. <!--PRTL-2543-->
* Updated the default Top Mutated Genes histogram on the Exploration page to display only COSMIC Genes by default. <!--PRTL-2371-->
* Added Follow-Ups tab and nested Molecular Tests to Case entity page. <!--PRTL-2704-->
* Added text to BAM slicing modal to instruct Users how to access unmapped reads. <!--PRTL-2620-->

### Bugs Fixed Since Last Release

* Fixed font in exported PNGs, SVGs to be consistent with the Portal UI. <!--PRTL-2653-->
* Made custom Case and File filters in the Repository page case insensitive. <!--PRTL-2571-->
* Fixed bug where pfam domains in Protein Viewer could not be clicked in Firefox. <!--PRTL-2511-->
* Fixed bug where TSV download button could not be clicked in MS Edge. <!--PRTL-2357-->
* Fixed controlled access alert pop-up in the Cart so that the modal disappears correctly once the User has successfully logged in and initiated the download. <!--PRTL-2429-->

### Known Issues and Workarounds

*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  Negative numbers may be displayed for the Missing value category in the Treatment node within a Clinical Analysis.  This occurs with projects that have multiple treatment nodes per case. All other values should be accurate. <!--SV-1604-->
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    * The footer says version 1.9, but it is actually 1.13
    *	Filtering by vital_status does not function in the Legacy Archive due to updates in how this property has been indexed.  A workaround is to perform the case level filtering in the GDC Data Portal and copy the filter string for use in the Legacy Archive or the legacy API. <!--SV-1508-->
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->

## Release 1.22.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  July 31, 2019

### New Features and Changes <!--REQ-387-->

* Replaced existing Clinical, Biospecimen columns on the Projects page with 4 columns: Clinical, Clinical Supplement, Biospecimen, Biospecimen Supplement. The Clinical and Biospecimen columns now link directly to the project page, and their counts indicate the total cases in the project. The Clinical Supplement and Biospecimen Supplement columns work the same as the old Clinical and Biospecimen columns - They link to the Repository page with Files filtered based on the Project and Data Category (Clinical or Biospecimen). <!--PRTL-2528-->
* Added a new icon to the GDC Apps menu, which links to the GDC Publications website page. <!--PRTL-2547-->
* Added the Synchronous Malignancy field to the Diagnoses / Treatments tab on the Case entity page. <!--PRTL-2582-->
* Added the Pack Years Smoked field to the Exposures tab on the Case entity page. <!--PRTL-2584-->
* Increased length of x-axis labels on histograms to 10 characters so that projects with names that are typically standard 10 chars will display fully (e.g. most TCGA projects like TCGA-BRCA). <!--PRTL-2598-->

### Bugs Fixed Since Last Release

* Fixed bug where the PNG, SVG files for the Overall Survival Plot could not be downloaded. <!--PRTL-2528-->

### Known Issues and Workarounds

*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Filtering by vital_status does not function in the Legacy Archive due to updates in how this property has been indexed.  A workaround is to perform the case level filtering in the GDC Data Portal and copy the filter string for use in the Legacy Archive or the legacy API. <!--SV-1508-->
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->   

## Release 1.21.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  June 5, 2019

### New Features and Changes <!--REQ-383-->

* Changed all Survival Plots to display the Duration (x-axis) in years instead of days. <!--PRTL-2404-->
* Updated data references to clinical properties throughout the Portal to match the underlying changes in the GDC data dictionary. <!--PRTL-2459-->

### Bugs Fixed Since Last Release

* Fixed bug where X-axis labels in histograms were cut off when displayed. <!--PRTL-1896-->
* Renamed the 'Experimental Strategies' facet on the Projects page to singular form. <!--PRTL-2262-->
* Fixed bug where columns with a % value of infinity (due to division by zero) show as 'NaN%'.  Replaced instead with a label of '--'. <!--PRTL-2384-->
* Fixed bug where the download button in the cart access banner was still disabled after a user logged in from the banner.  Instead, the experience is now improved so that after login, the banner is closed and the user must explicitly click 'Download' again. <!--PRTL-2393-->
* Fixed bug where if a new user logs into the Portal and views their profile, the app crashes if the user has no projects assigned yet. <!--PRTL-2529-->
* Fixed bug where Survival Rate numbers in the Survival Plot plot y-axis did not scale properly and overlapped into the axis lines. <!--PRTL-2530-->

### Known Issues and Workarounds

*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->   

## Release 1.20.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  April 17, 2019

### New Features and Changes <!--REQ-382-->

* Upgraded the Portal to use the latest React Javascript library (version 16.8) <!--PRTL-2440-->

### Bugs Fixed Since Last Release

* None

### Known Issues and Workarounds

*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    

## Release 1.19.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  February 20, 2019

### New Features and Changes <!--REQ-381-->

* Added support for viewing of controlled-access mutations in the Data Portal
* Added a new data access notification to remind logged-in users with access to controlled data that they need to follow their data use agreement.  The message is fixed at the top of the Portal.<!--PRTL-2400, PRTL-2434-->
* Added the ability to search for previous versions of files.  If the user enters the UUID of a previous version that cannot be found, the Portal returns the UUID of the latest version available. <!--PRTL-2387-->
* Renamed the Data Category for "Raw Sequencing Data" to "Sequencing Reads" throughout the portal where this appears, to be consistent with the Data Dictionary. <!--PRTL-118-->
* Added a link in the Portal footer to the GDC support page. <!--PRTL-2383-->

### Bugs Fixed Since Last Release

* Fixed bug where Survival Plot button never stops loading if plotting mutated vs. non-mutated cases for a single Gene. <!--PRTL-2398-->
* Fixed inconsistent button styling when downloading controlled Downstream Analyses Files from File Entity page. <!--PRTL-2395-->
* Removed unnecessary Survival column from Arrange Columns button on Case Entity, Gene Entity pages. <!--PRTL-2281-->
* Removed unnecessary whitespace from pie charts on Repository page. <!--PRTL-1923-->
* Added missing File Size unit to Clinical Supplement File, Biospecimen Supplement File tables on Case Entity page. <!--PRTL-2070-->
* Fixed bug where clicking on Case Counts in Projects Graph tab was going to the Repository Files tab instead of the Cases tab. <!--PRTL-2272-->
* Fixed bug where the counts shown beside customer filters on the Repository Cases tab were not updating when filtering on other facets. <!--PRTL-2412-->
* Fixed bug where clicking the # of Affected Cases denominator on the Gene page's Most Frequent Somatic Mutations table displayed an incorrect number of Cases.

### Known Issues and Workarounds

*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    

## Release 1.18.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  December 18, 2018

### New Features and Changes <!--REQ-335-->

* A new data access message has been added when downloading controlled data.  Users must agree to abide by data access control policies when downloading controlled data. <!--PRTL-2324,PRTL-2370,-->
* In the Mutation free-text search in Exploration, mutation display now includes the UUID, genomic location, and matched search term for easier mutation searching. <!--PRTL-1996, PRTL-2367-->
* The ability to sort on ranked columns has been made available. <!--PRTL-2274,PRTL-2336,PRTL-2365-->

### Bugs Fixed Since Last Release

*  In some cases, text was being cut off on the Project page visualization tab.  Text is no longer cut off. <!--PRTL-2290-->
*  HGNC link on Gene page broke as the source format url changed; The format was updated and the link is now functional <!--PRTL-2380-->
*  In the biospecimen details on the Case page, the cart icon would disappear once clicked.  It now is always visible. <!--PRTL-2282-->

### Known Issues and Workarounds

*  Pre-release Data Portal login is not supported on Internet Explorer or the last version of Edge (42).  Edge 41 does login successfully.
*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    

## Release 1.17.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  November 7, 2018

### New Features and Changes <!--REQ-334-->

* Copy Number Variation (CNV) data derived from GISTIC results are now available in the portal:
	* View number of CNV events on a gene in a cohort in the Explore Gene table tab <!-- PRTL-2259, PRTL-2328 -->
	* Explore CNVs associated with a gene on the Gene Entity Page <!-- PRTL-2252,PRTL-2273,PRTL-2344 -->
	* Explore CNVs concurrently with mutations on the Oncogrid with new visualization <!-- PRTL-2244,PRTL-2251,PRTL-2256,PRTL-2257,PRTL-2275,PRTL-2276, PRTL-2214, PRTL-2314, PRTL-2315, PRTL-2317, PRTL-2319, PRTL-2325, PRTL-2326, PRTL-2327,PRTL-2266  -->

### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds

*  Custom Facet Filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  SIFT and PolyPhen annotations are missing from the export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    

## Release 1.16.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  September 27, 2018

### New Features and Changes <!--REQ-333-->

* Updated Human Body Image to aggregate all current primary sites to available Major Primary Sites <!--PRTL-2248-->


### Bugs Fixed Since Last Release

* Fixed link on cart download error popup <!--PRTL-2284-->
* Updated Cancer Distribution table to have dropdown menus for primary_site and disease_type <!--PRTL-2286-->
* Updated Y-axis label on `Top Mutated Cancer Genes in Selected Projects Graph` <!--PRTL-2264-->
* Updated Set Operation Image to remove stray text <!--PRTL-2214-->

### Known Issues and Workarounds

*  Advanced Search
    * For advanced search and custom file facet filtering there are some properties that will appear as options that are no longer supported (e.g. file_state). <!--API-530-->
*  Custom facet filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  SIFT and PolyPhen annotations are missing from the Export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    


## Release 1.15.0

* __GDC Product__: GDC Data Portal
* __Release Date__:  August 23, 2018

### New Features and Changes <!--REQ-329-->

* File Versions are now visible in the "File Versions" section on the File Entity Page. <!--PRTL-2015-->
* "View Files in Repository" and "View Cases in Repository" button methods were updated to work faster. <!--PRTL-2003 PRTL-2106 -->

### Bugs Fixed Since Last Release

* Fixed warning messages that prompted users to login even when already logged in.  Error warnings now correctly prompt users to reference dbGAP for data access if already signed in. <!--PRTL-1937-->
* Fixed error where you could click Go on Case ID wildcard facet before inputting any data. <!--PRTL-2083-->
* Fixed cart header to be a consistent color for the whole table. <!--PRTL-1943-->
* Fixed error where you could save a set with no name or items, which resulted in an infinite spinner. <!--PRTL-2218--> <!--PRTL-1858-->
* Fixed table width issue when FM-AD was selected as a filter. <!--PRTL-2054-->
* Updated broken help link on Advanced Query. <!--PRTL-2218-->

### Known Issues and Workarounds

*  Advanced Search
    * For advanced search and custom file facet filtering there are some properties that will appear as options that are no longer supported (e.g. file_state). <!--API-530-->
*  Custom facet filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  SIFT and PolyPhen annotations are missing from the Export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    


## Release 1.14.0

* __GDC Product__: GDC Data Portal
* __Release Date__: June 13, 2018

### New Features and Changes <!--REQ-322-->

* Added new Experimental Strategies Diagnostic Slide Image, Bisulfite-Seq, ChIP-Seq, and ATAC-Seq to Case and Project entity pages.

### Bugs Fixed Since Last Release

* Fixed download of clinical and biospecimen data from the Repository when Case table rows are selected.


### Known Issues and Workarounds

*  Custom facet filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  SIFT and PolyPhen annotations are missing from the Export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  When user is logged in and try to download a controlled file he does not have access to, he's prompted to log in. He should be promted to request access. <!-- PRTL-1937 -->
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).


## Release 1.13.0

* __GDC Product__: GDC Data Portal
* __Release Date__: May 21, 2018

### New Features and Changes <!--REQ-322-->

*  Added new image viewer functionality for viewing tissue slide images <!--PRTL-2016-->


### Bugs Fixed Since Last Release

* Updated gene reference labels on gene entity page to adhere to preferred usage <!--PRTL-1911-->
* Fixed issue with user profile displaying all projects twice <!--PRTL-2035-->


### Known Issues and Workarounds

*  Custom facet filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  SIFT and PolyPhen annotations are missing from the Export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  When user is logged in and try to download a controlled file he does not have access to, he's prompted to log in. He should be promted to request access. <!-- PRTL-1937 -->
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    



Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).

## Release 1.12.0

* __GDC Product__: GDC Data Portal
* __Release Date__: February 15, 2018

### New Features and Changes <!--REQ-322-->

*  Provided the ability to export clinical and biospecimen data in a TSV format from the Case, Project, Exploration, Repository and Cart pages.<!--PRTL-1929-->
*  Removed from the Project entity page the sections about mutated genes, somatic mutations and affected Cases and replaced with a button "Explore data" that will open the Exploration page filtered on the project. Indeed the Exploration page provides the same information.  Added a breakdown of cases per primary site for a Project entity page with multiple primary sites (e.g. FM-AD). <!--PRTL-1903-->
*  Added display of coding DNA change and impacts for all the transcripts (instead of canonical transcript only) in the Mutation entity page - Consequences section. In the mutation table (e.g. in Repository), the impacts and consequences are displayed for the canonical transcript only. <!--PRTL-1927-->

### Bugs Fixed Since Last Release

*  Replaced the suggested set name when saving a set with selected items, e.g. for case set the suggested name is now "Custom Case selection". <!--PRTL-1911-->
*  Fixed the protein viewer to indicate when there are overlapping mutations. Mousing over the dot showing multiple mutations will open a right panel with the list of all the corresponding mutations.  <!--SV-750-->
*  Fixed Mutation entity page - Consequences table: the "Coding DNA Change" column is now populated for all the transcripts. <!--SV-751-->
*  Fixed download clinical and download biospecimen actions from TCGA-BRCA project.
*  Fixed facet behavior that did not reset back to showing all options after pressing reset-arrow. <!--PRTL-1928-->
*  Fixed error when user was trying to save a set with no value in the textbox "Save top:".  <!--PRTL-1909-->
*  Removed somatic mutation section from Case entity page for cases with no open-access mutation data (e.g. FM-AD or TARGET cases). <!--PRTL-1926-->
*  Fixed error where a blank page appears after unselecting `Cancer Gene Census` mutation facet. <!--PRTL-1933-->
*  Fixed duplicated date in sample sheet name (e.g. gdc_sample_sheet_YYYY-MM-DD_HH-MM.tsv.YYYY-MM-DD_HH-MM.tsv). <!--SV-942-->
*  Fixed error when annotations were not downloaded along with the file (in File entity page and Cart).  <!--SV-401-->

### Known Issues and Workarounds

*  Custom facet filters
    * Some definitions are missing from the property list when adding custom facet file or case filters. <!--SV-989-->
*  Visualizations
    *  SIFT and PolyPhen annotations are missing from the Export JSON of the mutation table. They are present in the export TSV. <!--PRTL-1990-->
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
*  Repository and Cart
    *  When user is logged in and try to download a controlled file he does not have access to, he's prompted to log in. He should be promted to request access. <!-- PRTL-1937 -->
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only. <!-- LGCY-81 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).


## Release 1.11.0

* __GDC Product__: GDC Data Portal
* __Release Date__: December 21, 2017

### New Features and Changes <!--REQ-322-->

* Updated UI to support SIFT and Polyphen annotations <!--PRTL-1404-->
* A `Sample Sheet` can now be created which allows easy association between file names and the case and sample submitter_id <!--PRTL-1872-->
* Updated Advanced Search page to include options to `Add All Files to Cart`, `Download Manifest`, and `View X Cases in Exploration` <!--PRTL-1796-->
* Provide clear message rather than blank screen if survival plots cannot be calculated for particular cohort comparison <!--PRTL-1842-->
* Display sample_type on associated entities section on file page <!--PRTL-1890-->
* Allows for special characters in case, gene, and mutation set upload (`-, :, >, .`) <!--PRTL-1847-->


### Bugs Fixed Since Last Release

*  Fixed error when trying to download large number of files from the Legacy Archive cart <!--LGCY-75-->
*  Fixed number of annotations displayed in Legacy Archive for particular entities <!--LGCY-74-->
*  Replaced missing bars to indicate proportion of applicable files and cases on project entity page in Cases and File Counts by Data Category table <!--PRTL-1725-->
*  Fixed project page display when projects are selected that contain no mutation data in the facet panel <!--PRTL-1754-->
*  Fixed error where exporting case sets as TSV included fewer cases than the total <!--PRTL-1888-->
*  Fixed error in exploration section when adding custom facets.  Previously selecting 'Only show fields with values' did not result in the expected behavior <!--PRTL-1901-->
*  Fixed error where number of associated entities for a file was showing an incorrect number <!--PRTL-1891-->

### Known Issues and Workarounds

*  Sample sheet will download with a file name including the date duplicated (e.g. gdc_sample_sheet_YYYY-MM-DD_HH-MM.tsv.YYYY-MM-DD_HH-MM.tsv)<!--SV-942-->
*  Custom facet filters
    * Definitions are missing from the property list when adding custom facet file or case filters <!--SV-916-->
*  Visualizations
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
    *  In the protein viewer there may be overlapping mutations.  In this case mousing over a point will just show a single mutation and the other mutations at this location will not be apparent.  <!--SV-750-->
*  Entity page
    *  On the mutation entity page, in the Consequences Table, the "Coding DNA Change" column is not populated for rows that do not correspond to the canonical mutation. <!-- SV-751 -->
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only.
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).



## Release 1.10.0

* __GDC Product__: GDC Data Portal
* __Release Date__: November 16, 2017

### New Features and Changes <!--REQ-322-->

* Support for uploading Case and Mutation sets in Exploration page <!--PRTL-1452, PRTL-1453-->
* Support for saving, editing, removing Case, Gene and Mutation sets in the Exploration page <!--PRTL-1472,PRTL-1464, PRTL-1473, PRTL-1468, PRTL-1465, PRTL-1470,PRTL-1469, PRTL-1466, PRTL-1471-->
* Added a Managed Sets menu where the user can see their saved sets <!--PRTL-1597-->
* Added an Analysis menu with two analyses: Set Operation and Cohort Comparison <!--PRTL-1599, PRTL-1600-->
* Added a User Profile page that shows all the projects and permissions assigned to the user: available in the username dropdown after the user logs in <!--PRTL-1458-->

### Bugs Fixed Since Last Release

*  Project page
    *  On the project page, the Summary Case Count link should open the case tab on the Repository page - instead it opens the file page <!--PRTL-1591-->

### Known Issues and Workarounds

*  Custom facet filters
    * Definitions are missing from the property list when adding custom facet file or case filters <!--SV-916-->
    * Selecting 'Only show fields with values' will show some fields without values in the Repository section.  This works correctly under the Exploration section. <!--SV-917-->
*  Visualizations
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
    *  In the protein viewer there may be overlapping mutations.  In this case mousing over a point will just show a single mutation and the other mutations at this location will not be apparent.  <!--SV-750-->
*  Entity page
    *  On the mutation entity page, in the Consequences Table, the "Coding DNA Change" column is not populated for rows that do not correspond to the canonical mutation. <!-- SV-751 -->
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only.
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).

## Release 1.9.0

* __GDC Product__: GDC Data Portal

* __Release Date__: October 24, 2017

### New Features and Changes <!--REQ-317-->

* Support for projects with multiple primary sites per project <!--PRTL-1478, PRTL-1676, PRTL-1668, PRTL-1675,PRTL-1685,PRTL-1695-->
* Support for slides that are linked to `sample` rather than `portion` <!--PRTL-1696-->

### Bugs Fixed Since Last Release

None

### Known Issues and Workarounds
*  Visualizations
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
    *  In the protein viewer there may be overlapping mutations.  In this case mousing over a point will just show a single mutation and the other mutations at this location will not be apparent.  <!--SV-750-->
*  Project page
    *  On the project page, the Summary Case Count link should open the case tab on the Repository page - instead it opens the file page <!--PRTL-1591-->
*  Entity page
    *  On the mutation entity page, in the Consequences Table, the "Coding DNA Change" column is not populated for rows that do not correspond to the canonical mutation. <!-- SV-751 -->
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only.
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).


## Release 1.8.0

* __GDC Product__: GDC Data Portal
* __Release Date__: August 22, 2017

### New Features and Changes

Major features/changes:

* A feature that links the exploration and repository pages was added. For example:
    - In the exploration page, cases with a specific mutation could be selected. This set could then be linked to the repository page to download the data files associated with these cases.
    - In the repository menu, the user can select cases associated with specific files. The set could then be linked to exploration page to view the variants associated with this set of cases.

*  Users can now upload a custom gene list to the exploration page and leverage the GDC search and visualization features for cases and variants associated with the gene set.

*  Filters added for the gene entity page. For example:
    - Clicking on a mutated gene from the project page will display mutations associated with the gene that are present in this project (filtered protein viewer, etc.).
    - Clicking on a mutated gene from the exploration page will display the mutations associated with the gene filtered by additional search criteria, such as "primary site is Kidney and mutation impact is high".

* UUIDs are now hidden from tables and charts to simplify readability. The UUIDs can still be exported and viewed in the tables using the "arrange columns" feature. In the mutation table, UUIDs are automatically exported.

* Mutation entity page - one consequence per transcript is shown (10 rows by default) in the consequence table. The user should display all rows before exporting the table.

### Bugs Fixed Since Last Release
*  Exploration
    *  Combining "Variant Caller" mutation filter with a case filter will display incorrect counts in the mutation facet. The number of mutations in the resulting mutation table is correct. <!-- API-307 -->
    *  Mutation table: it is difficult to click on the denominator in "#Affected Cases in Cohort" column displayed to the left side of the bar. The user should click at a specific position at the top of the number to be able to go to the corresponding link. <!-- PRTL-1377 -->


### Known Issues and Workarounds
*  Visualizations
    *  Data Portal graphs cannot be exported as PNG images in Internet Explorer. Graphs can be exported in PNG or SVG format from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet Explorer does not display chart legend and title when re-opening previously downloaded SVG files, the recommendation is to open downloaded SVG files with another program.
    *  In the protein viewer there may be overlapping mutations.  In this case mousing over a point will just show a single mutation and the other mutations at this location will not be apparent.  <!--SV-750-->
*  Project page
    *  On the project page, the Summary Case Count link should open the case tab on the Repository page - instead it opens the file page <!--PRTL-1591-->
*  Entity page
    *  On the mutation entity page, in the Consequences Table, the "Coding DNA Change" column is not populated for rows that do not correspond to the canonical mutation. <!-- SV-751 -->
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only.
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).


## Release 1.6.0

* __GDC Product__: GDC Data Portal
* __Release Date__: June 29, 2017

### New Features and Changes

There was a major new release of the GDC Data Portal focused on Data Analysis, Visualization, and Exploration (DAVE). Some important new features include the following:

*  New visual for the Homepage: a human body provides the number of Cases per Primary Site with a link to an advanced Cancer Projects search
*  The Projects menu provides the Top 20 Cancer Genes across the GDC Projects and the Case Distribution per Project
*  A new menu "Exploration" is an advanced Cancer Projects search which provides the ability to apply Case, Gene, and Mutation filters to look for:
    * List of Cases with the largest number of Somatic Mutations
    *	The most frequently mutated Genes
    *	The most frequent Variants
    *	Oncogrid view of mutation frequency
*  Visualizations are provided across the Project, Case, Gene and Mutation entity pages:
    *	List of most frequently mutated genes and most frequent variants
    *	Survival plots for patients with or without specific variants
    *	Survival plots for patients with or without variants in specific genes
    *	Lollipop plots of mutation frequency across protein domains
*  Links to external databases (COSMIC, dbSNP, Uniprot, Ensembl, OMIM, HGNC)
*  Quick Search for Gene and Mutation entity pages
*  The ability to export the current view of a table in TSV
*  Retired GDC cBioPortal

_For detailed updates please review the [Data Portal User Guide](../Users_Guide/getting_started/)._

### Bugs Fixed Since Last Release

*  BAM Slicing dialog box does not disappear automatically upon executing the BAM slicing function. The box can be closed manually. <!-- PRTL-282 -->
*  Very long URLs will produce a 400 error.  Users may encounter this after clicking on "source files" on a file page where the target file is derived from hundreds of other files such as for MAF files. <!-- SV-396 / PRTL-342-->
*  If bam slicing produces an error pop-up message it will be obscured behind the original dialog box. <!--SV-419-->
    *   Internet Explorer users are not able to use the "Only show fields with no values" when adding custom facets <!-- PGDC-2467 / PRTL-109 -->
    *   Exporting large tables in the Data Portal may produce a 500 error.  Filtering this list to include fewer cases or files should eliminate the error <!--API-223-->    

### Known Issues and Workarounds
*  New Visualizations
    *  Cannot export Data Portal graphs in PNG in Internet Explorer. Graphs can be exported to PNG or SVG from Chrome or Firefox browsers <!-- PRTL-1325 / PRTL-1114 -->. Internet would not display chart legend and title when re-opening previously downloaded SVG files, recommendation is to open downloaded SVG files with another software.
    *  In the protein viewer there may be overlapping mutations.  In this case mousing over a point will just show a single mutation and the other mutations at this location will not be apparent.  <!--SV-750-->
*  Exploration
    *  Combining "Variant Caller" mutation filter with a case filter will display wrong counts in the mutation facet. The number of mutations in the result mutation table is correct. <!-- API-307 -->
    *  Mutation table: it is difficult to click on the denominator in "#Affected Cases in Cohort" column displayed to the left side of the bar. The user should click at a specific position at the top of the number to be able to go to the corresponding link. <!-- PRTL-1377 -->
*  Entity page
    *  On the mutation entity page, in the Consequences Table, the "Coding DNA Change" column is not populated for rows that do not correspond to the canonical mutation. <!-- SV-751 -->
*  Repository and Cart
    *  The annotation count in File table of Repository and Cart does not link to the Annotations page anymore. The user can navigate to the annotations through the annotation count in Repository - Case table.
*  Legacy Archive
    *	Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.
    *	Exporting the Cart table in JSON will export the GDC Archive file table instead of exporting the files in the Cart only.
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).



## Release 1.5.2

* __GDC Product__: GDC Data Portal
* __Release Date__: May 9, 2017

### New Features and Changes

* Removed link to Data Download Statistics Report <!--PRTL-1081-->
* Updated version numbers of API, GDC Data Portal, and Data Release

### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds

*   General
    *   Exporting large tables in the Data Portal may produce a 500 error.  Filtering this list to include fewer cases or files should eliminate the error <!--API-223-->
    *   After successful authentication, the authentication popup does not close for Internet Explorer users running in "Compatibility View". Workaround is to uncheck "Display Intranet sites in Compatibility View" in Internet Explorer options. Alternatively, refreshing the portal will correctly display authentication status. <!-- PGDC-2403 / PRTL-133 -->
    *   BAM Slicing dialog box does not disappear automatically upon executing the BAM slicing function. The box can be closed manually. <!-- PRTL-282 -->
    *   Due to preceding issue, If bam slicing produces an error pop-up message it will be obscured behind the original dialog box. <!--SV-419-->
    *   Very long URLs will produce a 400 error.  Users may encounter this after clicking on "source files" on a file page where the target file is derived from hundreds of other files such as for MAF files.  To produce a list of source files an API call can be used with the search parameter "fields=analysis.input_files.file_name". <!-- SV-396 / PRTL-342-->
		*   Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.


Example

    https://api.gdc.cancer.gov/files/455e26f7-03f2-46f7-9e7a-9c51ac322461?pretty=true&fields=analysis.input_files.file_name




*   Cart
    *   Counts displayed in the top right of the screen, next to the Cart icon, may become inconsistent if files are removed from the server. <!-- PGDC-2403 / PRTL-133 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   Internet Explorer users are not able to use the "Only show fields with no values" when adding custom facets <!-- PGDC-2467 / PRTL-109 -->
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).





## Release 1.4.1

* __GDC Product__: GDC Data Portal
* __Release Date__: October 31, 2016

### New Features and Changes

* Added a search feature to help users select values of interest in certain facets that have many values. <!-- PRTL-21 -->
* Added support for annotation ID queries in quick search. <!-- PRTL-29 -->
* Added a warning when a value greater than 90 is entered in the "Age at Diagnosis" facet. <!-- PRTL-77 -->
* Added Sample Type column to file entity page. <!-- PRTL-42 -->
* Authentication tokens are refreshed every time they are downloaded from the GDC Data Portal. <!-- PRTL-278 -->
* Buttons are inactive when an action is in progress. <!-- PRTL-270 -->
* Improved navigation features in the overview chart on portal homepage. <!-- PRTL-2 -->
* Removed State/Status from File and Case entity pages <!-- PRTL-292 -->
* Removed the "My Projects" feature. <!-- PRTL-174 -->
* Removed "Created" and "Updated" dates from clinical and biospecimen entities. <!-- PRTL-3 -->

### Bugs Fixed Since Last Release

*  Advanced search did not accept negative values for integer fields. <!-- PRTL-283 -->
*  Moving from facet search to advanced search resulted in an incorrect advanced search query. <!-- PRTL-284 -->
*  Some facets were cut off in Internet Explorer and Firefox. <!-- PRTL-290 -->

### Known Issues and Workarounds

*   General
    *   Exporting large tables in the Data Portal may produce a 500 error.  Filtering this list to include fewer cases or files should eliminate the error <!--API-223-->
    *   After successful authentication, the authentication popup does not close for Internet Explorer users running in "Compatibility View". Workaround is to uncheck "Display Intranet sites in Compatibility View" in Internet Explorer options. Alternatively, refreshing the portal will correctly display authentication status. <!-- PGDC-2403 / PRTL-133 -->
    *   BAM Slicing dialog box does not disappear automatically upon executing the BAM slicing function. The box can be closed manually. <!-- PRTL-282 -->
    *   Due to preceding issue, If bam slicing produces an error pop-up message it will be obscured behind the original dialog box. <!--SV-419-->
    *   Very long URLs will produce a 400 error.  Users may encounter this after clicking on "source files" on a file page where the target file is derived from hundreds of other files such as for MAF files.  To produce a list of source files an API call can be used with the search parameter "fields=analysis.input_files.file_name". <!-- SV-396 / PRTL-342-->
		*   Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.


Example

    https://api.gdc.cancer.gov/files/455e26f7-03f2-46f7-9e7a-9c51ac322461?pretty=true&fields=analysis.input_files.file_name




*   Cart
    *   Counts displayed in the top right of the screen, next to the Cart icon, may become inconsistent if files are removed from the server. <!-- PGDC-2403 / PRTL-133 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   Internet Explorer users are not able to use the "Only show fields with no values" when adding custom facets <!-- PGDC-2467 / PRTL-109 -->
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).




## Release 1.3.0

* __GDC Product__: GDC Data Portal
* __Release Date__: September 7, 2016

### New Features and Changes

*   A new "Metadata" button on the cart page to download merged clinical, biospecimen, and file metadata in a single consolidated JSON file. **May require clearing browser cache** <!-- PRTL-177 -->
*   Added a banner on the Data Portal to help users find data <!-- PRTL-237 -->
*   Added support for "Enter" key on login button <!-- PRTL-136 -->
*   On the Data page, the browser will remember which facet tab was selected when hitting the "Back" button <!-- PRTL-147 -->
*   In file entity page, if there is a link to one single file, redirect to this file's entity page instead of a list page.  <!-- PRTL-75 -->


### Bugs Fixed Since Last Release

* 	Adding a mix of open and controlled files to the cart from any Case entity pages was creating authorization issues <!-- PRTL-226 -->
*  Opening multiple browser tabs and adding files in those browser tabs was not refreshing the cart in other tabs. <!-- PRTL-181 -->
*   When user logs in from the advanced search page, the login popup does not automatically close <!-- PRTL-88 -->
*   When removing a file from the cart and clicking undo, GDC loses track of permission status of the user towards this file and will ask for the user to log-in again. <!-- PGDC-2496 -->
*   Download File Metadata button produces incomplete JSON output omitting such fields as file_name and submitter_id.  The current workaround includes using the API to return file metadata. <!-- SV-359 / PRTL-177 -->    
*   Annotations notes do not wrap to the next line at the beginning or the end of a word, some words might be split in two lines <!-- PGDC-2474 / PRTL-182-->   
*   Sorting annotations by Case UUID causes error <!-- PRTL-95 -->   

### Known Issues and Workarounds

*   General
    *   When no filters are engaged in the Legacy Archive or Data Portal, clicking the Download Manifest button may produce a 500 error and the message "We are currently experiencing issues. Please try again later.".  To avoid this error the user can first filter by files or cases to reduce the number files added to the manifest.
    *   After successful authentication, the authentication popup does not close for Internet Explorer users running in "Compatibility View". Workaround is to uncheck "Display Intranet sites in Compatibility View" in Internet Explorer options. Alternatively, refreshing the portal will correctly display authentication status. <!-- PGDC-2403 / PRTL-133 -->
    *   BAM Slicing dialog box does not disappear automatically upon executing the BAM slicing function. The box can be closed manually. <!-- PRTL-282 -->
    *   Due to preceding issue, If bam slicing produces an error pop-up message it will be obscured behind the original dialog box. <!--SV-419-->
    *   Very long URLs will produce a 400 error.  Users may encounter this after clicking on "source files" on a file page where the target file is derived from hundreds of other files such as for MAF files.  To produce a list of source files an API call can be used with the search parameter "fields=analysis.input_files.file_name". <!-- SV-396 -->
    *   On the Legacy Archive, searches for "Case Submitter ID Prefix" containing special characters are not displayed correctly above the result list. The result list is correct, however. <!--SV-412 / LGCY-33-->

Example

    https://api.gdc.cancer.gov/files/455e26f7-03f2-46f7-9e7a-9c51ac322461?pretty=true&fields=analysis.input_files.file_name




*   Cart
    *   Counts displayed in the top right of the screen, next to the Cart icon, may become inconsistent if files are removed from the server. <!-- PGDC-2403 / PRTL-133 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   Internet Explorer users are not able to use the "Only show fields with no values" when adding custom facets <!-- PGDC-2467 / PRTL-109 -->
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibility mode. <!-- PGDC-2480 -->    


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).

## Release 1.2.0

* __GDC Product__: GDC Data Portal
* __Release Date__: August 9th, 2016

### New Features and Changes

*   Added a retry (1x) mechanism for API calls <!-- PGDC-2393 -->
*   Added support for ID fields in custom facets <!-- PGDC-1222 -->
*   Added Case Submitter ID to the Annotation entity page <!-- PGDC-747 -->
*   Added a link to Biospeciment in the Case entity page <!-- PGDC-2346 -->

### Bugs Fixed Since Last Release

*   General.
    *   Not possible to use the browser's back button after hitting a 404 page <!-- PGDC-2429 -->
    *   404 page missing from Legacy Archive Portal <!-- PGDC-2477 -->
    *   Table widget icon and export JSON icon should be different <!-- PGDC-2446 -->    
    *   Download SRA XML files from the legacy archive portal might not be possible in some context <!-- PGDC-2457 --> <!-- PGDC-2469 -->
*   Data and facets
    *   Default values for age at diagnosis is showing 0 to 89 instead of 0 to 90 <!-- PGDC-2478 -->
    *   Biospecimen search in the case entity page does not highlight (but does bold and filter) results in yellow when title case is not followed <!-- PGDC-2451 -->
    *   Table sorting icon does not include numbers <!-- PGDC-35 -->   
    *   '--' symbol is missing on empty fields (blank instead), additional missing fields identified since last release.  <!-- PGDC-2447 -->    
### Known Issues and Workarounds

*   General
    *   When no filters are engaged in the Legacy Archive or Data Portal, clicking the Download Manifest button may produce a 500 error and the message "We are currently experiencing issues. Please try again later.".  To avoid this error the user can first filter by files or cases to reduce the number files added to the manifest.
    *   After successful authentication, the authentication popup does not close for Internet Explorer users running in "Compatibility View". This only impact users at the NIH. Workaround is to uncheck "Display Intranet sites in Compatibility View" in Internet Explorer options. Alternatively, refreshing the portal will correctly display authentication status.
    *   When user login from the advanced search page, the login popup does not automatically close <!-- PRTL-88 -->
*   Cart
    *   When removing a file from the cart and clicking undo, GDC looses track of permission status of the user towards this file and will ask for the user to log-in again. <!-- PGDC-2496 -->     
    *   Counts displayed in the top right of the screen, next to the Cart icon, might get inconsistent if files are removed from the server. <!-- PGDC-2403 -->
    *   Download File Metadata button produces incomplete JSON output omitting such fields as file_name and submitter_id.  The current workaround includes using the API to return file metadata. <!-- SV-359 -->
*   Annotations
    *   Annotations notes do not wrap to the next line at the beginning or the end of a word, some words might be split in two lines <!-- PGDC-2474 -->
    *   Sorting annotations by Case UUID causes error <!-- PRTL-95 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   Internet Explorer users are not able to use the "Only show fields with no values" when adding custom facets <!-- PGDC-2467 -->
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibilty mode <!-- PGDC-2480 -->    


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).


## Release 1.1.0

* __GDC Product__: GDC Data Portal
* __Release Date__: June 1st, 2016

### New Features and Changes

*   This is a bug-fixing release, no new features were added.

### Bugs Fixed Since Last Release

*   General
    *  Fixed 508 compliance issues. <!-- PGDC-2497 --><!-- PGDC-2431 -->
    *  Disabled download manifest action on projects without files.  <!-- PGDC-2416 -->
    *  Updated the portal to indicate to the user that his session expired when he tries to download the authentication token. <!-- PGDC-2455 -->   
    *  Unselected "My project" filter after user logs-in. <!-- PGDC-2462 -->  
    *  Fixed missing padding when query includes "My Projects". <!-- PGDC-2420 -->
    *  Enforced "Add to cart" limitation to 10,000 files everywhere on the Data Portal. <!-- PGDC-2409 -->
*   Tables
    *  Improved usability of the "Sort" feature  <!-- PGDC-1771 -->
    *  Updated the "Add all files to cart" button to add all files corresponding to the current query (and not only displayed files). <!-- PGDC-2439 -->
    *  Fixed an issue where Platform would show "0" when selected platform is "Affymetrix SNP 6.0". <!-- PGDC-2419 -->
*   Data
    *  Corrected default values populated when adding a custom range facet. <!-- PGDC-2445 --> <!-- PGDC-2444 -->  <!-- PGDC-2225 -->    
    *  Fixed an issue preventing the user to sort by File Submitter ID in data tables.<!-- PGDC-2430 -->    
*   File Entity Page
    *  Improved "Associated Cases/Biospecimen" table for files associated to a lot of cases.  <!-- PGDC-1993 -->    
    *  Fixed an error when performing BAM Slicing. <!-- PGDC-2433 -->

### Known Issues and Workarounds

*   General.
    *   After successful authentication, the authentication popup does not close for Internet Explorer users running in "Compatibility View". This only impact users at the NIH. Workaround is to uncheck "Display Intranet sites in Compatibility View" in Internet Explorer options. Alternatively, refreshing the portal will correctly display authentication status.
    *   Download SRA XML files from the legacy archive portal might not be possible in some context <!-- PGDC-2457 --> <!-- PGDC-2469 -->  
    *   Not possible to use the browser's back button after hitting a 404 page <!-- PGDC-2429 -->
    *   404 page missing from Legacy Archive Portal <!-- PGDC-2477 -->
    *   Table widget icon and export JSON icon should be different <!-- PGDC-2446 -->     
*   Data and facets
    *   Default values for age at diagnosis is showing 0 to 89 instead of 0 to 90 <!-- PGDC-2478 -->
    *   Biospecimen search in the case entity page does not highlight (but does bold and filter) results in yellow when title case is not followed <!-- PGDC-2451 -->
    *   Table sorting icon does not include numbers <!-- PGDC-35 -->    
    *   '--' symbol is missing on empty fields (blank instead), additional missing fields identified since last release.  <!-- PGDC-2447 -->    
*   Cart
    *   When removing a file from the cart and clicking undo, GDC looses track of permission status of the user towards this file and will ask for the user to log-in again. <!-- PGDC-2496 -->     
    *   Counts displayed in the top right of the screen, next to the Cart icon, might get inconsistent if files are removed from the server. <!-- PGDC-2403 -->
*   Annotations
    *   Annotations notes do not wrap to the next line at the beginning or the end of a word, some words might be split in two lines <!-- PGDC-2474 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).
    *   Internet Explorer users are not able to use the "Only show fields with no values" when adding custom facets <!-- PGDC-2467 -->
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibilty mode <!-- PGDC-2480 -->    


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).


## Release 1.0.1

* __GDC Product__: GDC Data Portal
* __Release Date__: May 18, 2016

### New Features and Changes

*   This is a bug-fixing release, no new features were added.

### Bugs Fixed Since Last Release

*   Tables and Export
    *   Restore default table column arrangement does not restore to the default but it restores to the previous state <!-- PGDC-1769 -->
*   Cart and Download
    *   Make the cart limit warning message more explanatory <!-- PGDC-1952 -->   
    *   In some situations, adding filtered files to the cart might fail <!-- PGDC-1981 -->   
*   Layout, Browser specific and Accessibility
    *   When disabling CSS, footer elements are displayed out of order <!-- PGDC-1972 -->
    *   If javascript is disabled html tags are displayed in the warning message <!-- PGDC-1835 -->
    *   Layout issues when using the browser zoom in function on tables <!-- PGDC-116 -->
    *   Cart download spinner not showing at the proper place <!-- PGDC-2056 -->
    *   Not all facets are expanded by default when loading the app <!-- PGDC-2061 -->

### Known Issues and Workarounds

*   General
    *   If a user has previously logged into the Portal and left a session without logging out, if the user returns to the Portal after the user's sessionID expires, it looks as if the user is still authenticated. The user cannot download the token and gets an error message that would not close. The user should clear the cache to properly log out.
    *   '--' symbol is missing on empty fields (blank instead) <!-- PGDC-2418 -->
    *   Download manifest button is available for TARGET projects with 0 files, resulting in error if user clic on button <!-- PGDC-2416 -->
    *   After successful authentication, the authentication popup does not close for Internet Explorer users running in "Compatibility View". This only impact users at the NIH. Workaround is to uncheck "Display Intranet sites in Compatibility View" in Internet Explorer options. Alternatively, refreshing the portal will correctly display authentication status.
*   Data    
    *   When adding a custom range facet, default values are incorrectly populated <!-- PGDC-2445 --> <!-- PGDC-2444 -->  <!-- PGDC-2225 -->
    *   The portal might return incorrect match between cases and files when using field cases.samples.portions.created_datetime (custom facet or advanced search). Note: this is not a UI issue. <!-- PGDC-2440 -->
    *   Sorting File Submitter ID option on the file tab result in a Data Portal Error <!-- PGDC-2430 -->
*   Tables and Export
    *   Table sorting icon does not include numbers <!-- PGDC-35 -->    
*   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.cancer.gov/about-gdc/gdc-faqs).

Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).
