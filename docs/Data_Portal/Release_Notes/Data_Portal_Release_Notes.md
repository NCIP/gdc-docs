# Data Portal Release Notes

| Version | Date |
|---|---|
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

_For detailed updates please review the [Data Portal User Guide](../Data_Portal/Users_Guide/Getting_Started/)._

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
