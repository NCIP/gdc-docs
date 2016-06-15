# Data Portal Release Notes

## Release 0.3.24-spr5

* __GDC Product__: GDC Data Portal
* __Release Date__: April 18, 2016

### New Features and Changes

*   This is a bug-fixing release, no new features were added.

### Bugs Fixed Since Last Release

*   Tables and Export
    *   Table export return correct content (XML, TSV) although file extension is incorrect (always JSON) <!-- PGDC-1986 -->
    *   Nested values are not exported properly <!-- PGDC-910 SV-13 -->
    *   Project table summary chart does not display colors with Internet Explorer 11 <!-- PGDC-1953 -->
*   Cart and Download
    *   When saving a file for download, the download pop-up indicates "from:blob:" <!-- PGDC-115 -->
    *   User with update only can also download files also role should not allow <!-- PGDC-1782 -->
    *   User can add more files to the cart than supported if adding files one by one <!-- PGDC-1957 -->
    *   In very specific situations, the cart can display inconsistent data <!-- PGDC-1028 -->   
    *   In some situations the cart will display the authentication window for authenticated users when trying to download <!-- PGDC-1959 -->  
*   Search
    *   In the data page, some pie chart titles are too long <!-- PGDC-914 -->
    *   Facets displaying histogram do not display a mouseover tooltip if the value is very low <!-- PGDC-913 -->
*   Layout, Browser specific and Accessibility
    *   No favicon is displayed on Internet Explorer <!-- PGDC-1948 -->
    *   Download cart button does not function properly with Safari 9.0.3 <!-- PGDC-1978 -->
    *   Layout issue when browser is reduced to a small window size <!-- PGDC-16 -->


### Known Issues and Workarounds

*   General
    *   If a user has previously logged into the Portal and left a session without logging out, if the user returns to the Portal after the user's sessionID expires, it looks as if the user is still authenticated. The user cannot download the token and gets an error message that would not close. The user should clear the cache to properly log out. 
*   Tables and Export
    *   Table sorting icon does not include numbers <!-- PGDC-35 -->
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
*   Non UI-related tickets
    *   Investigations were done on this issue "Authenticated user with access to TARGET projects cannot download Target DCC files. They will get a pop-up message." was related  to ACL permissions (not UI related).
    *   Investigations were done on this issue "Missing pathology files and slide images for some TCGA datasets since they are not connected to the biospecimen chain yet" was related to Data not imported (not UI related).
    *   Associated entities is empty for some files (note: this is a data issue) <!-- PGDC-1887 -->


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).

## Release 0.3.24.2

* __GDC Product__: GDC Data Portal
* __Release Date__: April 4, 2016

### New Features and Changes

*   Implemented a homepage to provide high-level information about the GDC Data Portal
*   Created the GDC Legacy Archive to provide access to legacy data
*   Updated the GDC Data Portal to support harmonized data
*   Updated the login process to be consistent with the GDC Data Submission Portal login
*   Added the ability to export clinical and biospecimen data from Project and Case entity pages (replacing download clinical XML)
*   Improved GDC Data Portal search capabilities:
    *   Added customized facets, expanding the GDC Data Portal search capabilities to non-default facets
    *   Updated the list of fields suggested in the Advanced Search and Custom Facets
    *   Added prefix search on case submitter ID on the Facet panel
    *   Improved operator auto-suggestion based on field type in the Advanced Search
*   Improved the download experience by letting the browser handle the download (allow canceling)
*   Enhanced Cart capabilities:
    *   Added the ability to download metadata files (SRA XML and MAGE-TAB)
    *   Added the ability to download biospecimen and clinical data
*   Improved the Biospecimen tree in the Case entity page (users can search for ID and expand the tree)
*   Improved support for Section 508 standards and the Internet Explorer browser
*   Improved the table view in the Projects page (Added a "Total" row with links to the search page)
*   Improved the BAM Slicing UI (added an example and the ability to use tabulations for BED format in the text box)

### Bugs Fixed Since Last Release

*   Disease Type does not auto filter (Data Facet)
*   When exporting the annotation table to a file, the create date field is not displayed in the same format than the UI
*   This ticket is not applicable anymore: "In the Case Entity Page, the "Download Clinical XML" button is not applicable to TARGET projects. The label is incorrect but the functionality is accurate.". The Case entity page has a standard download clinical feature.


### Known Issues and Workarounds

*   General
    *   If a user has previously logged into the Portal and left a session without logging out, if the user returns to the Portal after the user's sessionID expires, it looks as if the user is still authenticated. The user cannot download the token and gets an error message that would not close. The user should clear the cache to properly log out. 
*   Tables and Export
    *   Table sorting icon does not include numbers <!-- PGDC-35 -->
    *   Table export return correct content (XML, TSV) although file extension is incorrect (always JSON) <!-- PGDC-1986 -->
    *   Nested values are not exported properly <!-- PGDC-910 SV-13 -->
    *   Project table summary chart does not display colors with Internet Explorer 11 <!-- PGDC-1953 -->
    *   Exporting a table containing a very large number of records might trigger an exception (when the server is timing out)<!-- PGDC-1962 -->
    *   Restore default table column arrangement does not restore to the default but it restores to the previous state <!-- PGDC-1769 -->
*   Cart and Download
    *   When saving a file for download, the download pop-up indicates "from:blob:" <!-- PGDC-115 -->
    *   User with update only can also download files also role should not allow <!-- PGDC-1782 -->
    *   User can add more files to the cart than supported if adding files one by one <!-- PGDC-1957 -->
    *   In very specific situations, the cart can display inconsistent data <!-- PGDC-1028 -->   
    *   Make the cart limit warning message more explanatory <!-- PGDC-1952 -->   
    *   In some situations the cart will display the authentication window for authenticated users when trying to download <!-- PGDC-1959 -->   
*   Search
    *   In the data page, some pie chart titles are too long <!-- PGDC-914 -->
    *   Facets displaying histogram do not display a mouseover tooltip if the value is very low <!-- PGDC-913 -->
*   Layout, Browser specific and Accessibility
    *   No favicon is displayed on Internet Explorer <!-- PGDC-1948 -->
    *   When disabling CSS, footer elements are displayed out of order <!-- PGDC-1972 -->
    *   Download cart button does not function properly with Safari 9.0.3 <!-- PGDC-1978 -->
    *   If javascript is disabled html tags are displayed in the warning message <!-- PGDC-1835 -->
    *   Layout issue when browser is reduced to a small window size <!-- PGDC-16 -->
    *   Layout issues when using the browser zoom in function on tables <!-- PGDC-116 -->
*   Non UI-related tickets
    *   Investigations were done on this issue "Authenticated user with access to TARGET projects cannot download Target DCC files. They will get a pop-up message." was related  to ACL permissions (not UI related).
    *   Investigations were done on this issue "Missing pathology files and slide images for some TCGA datasets since they are not connected to the biospecimen chain yet" was related to Data not imported (not UI related).
    *   Associated entities is empty for some files (note: this is a data issue) <!-- PGDC-1887 -->


Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).

## Release 0.2.18.3

* __GDC Product__: GDC Data Portal
* __Release Date__: November 23, 2015

### New Features and Changes

* BAM slicing UI in the File Entity page (only available for BAM files with BAI files)

### Bugs Fixed Since Last Release

* Disease Type does not auto filter (Data Facet)

### Known Issues and Workarounds

* Missing pathology files and slide images for some TCGA datasets since they are not connected to the biospecimen chain yet
* In some specific situation (complex queries using must_not term), the advanced search might returns incorrect results
* When exporting the annotation table to a file, the create date field is not displayed in the same format than the UI
* In the advanced search, adding a trailing space to the query causes the UI to display a syntax error although the query is valid
* Checksum missing for MAGE-TAB files
* User can download large files (> 1GB) from the Portal but it fails when extracting the tar.gz file. User should use the GDC Download Transfer Tool for large downloads
* Authenticated user with access to TARGET projects cannot download Target DCC files. They will get a pop-up message.
* In Case Entity Page, the "Download Clinical XML" button is not applicable to TARGET projects. The label is incorrect but the functionality is accurate.
* The advanced search suggests more fields than expected (e.g. "cases.aliquots" should not be displayed). User should remove the field from the query and use the long field name (e.g. "cases.samples.portions.analytes.aliquots")
* If a user has previously logged into the Portal and left a session without logging out, if he comes back to the Portal after his sessionID expires, it looks like the user is still authenticated. The user cannot download the token and gets an error message that would not close. The user should clear the cache to properly log out. 
* Sorting on size and removing files from cart does not work in IE 10
* Annotations page not section 508 compliant
* Facets are displayed above the results table when window is small
* Table sorting icon does not include numbers

Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).

## Release 0.2.15-oicr4

* __GDC Product__: GDC Data Portal
* __Release Date__: October 1, 2015

### New Features and Changes

*   n/a

### Bugs Fixed Since Last Release

*   Data access and subtype units are correct in the data download statistics report
*   Improved 508 compliance of the portal (including the data download statistics report)
*   Addressed an issue with the user getting an empty page in specific browsing situations
*   "My projects" filter has been re-enabled for users authenticated with eRACommons but without dbGap authorization.

### Known Issues and Workarounds

*   Missing pathology files and slide images for some TCGA datasets since they are not connected to the biospecimen chain yet
*   In some specific situation (complex queries using must_not term), the advanced search might returns incorrect results
*   When exporting the annotation table to a file, the create date field is not displayed in the same format than the UI
*   In the advanced search, adding a trailing space to the query causes the UI to display a syntax error although the query is valid
*   Checksum missing for MAGE-TAB files
*   User can download large files (> 1GB) from the Portal but it fails when extracting the tar.gz file. User should use the GDC Download Transfer Tool for large downloads
*   Authenticated user with access to TARGET projects cannot download Target DCC files. They will get a pop-up message.
*   In Case Entity Page, the "Download Clinical XML" button is not applicable to TARGET projects. The label is incorrect but the functionality is accurate.
*   The advanced search suggests more fields than expected (e.g. "cases.aliquots" should not be displayed). User should remove the field from the query and use the long field name (e.g. "cases.samples.portions.analytes.aliquots")
*   If a user has previously logged into the Portal and left a session without logging out, if he comes back to the Portal after his sessionID expires, it looks like the user is still authenticated. The user cannot download the token and gets an error message that would not close. The user should clear the cache to properly log out. 
*   Sorting on size and removing files from cart does not work in IE 10
*   Annotations page not section 508 compliant
*   Facets are displayed above the results table when window is small
*   Table sorting icon does not include numbers

Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).

## Release 0.2.15-oicr2

* __GDC Product__: GDC Data Portal
* __Release Date__: August 31, 2015

### New Features and Changes

*   Authentication, Authorization & Pop-up messages:
    *   If a user tries to download a related file that is controlled from the File Entity Page, he will get a pop-up message with appropriate guidance ("Login" or "No access to the file")
    *   If a user tries to download a controlled clinical file from the Case page, he will get a pop-up message to indicate that he does not have access to the file
    *   If a user authenticates to the Portal with an eRa account without dbGap authorization, he will get a warning message. Then he will have the choice to logout or to continue browsing the Portal with his eRa account. If he continues browsing the Portal with his eRa account, "My projects" filter will be hidden (temporary solution).
*   If User downloads large files (from the Cart or from the File table), the Portal displays a spinner to indicate the download is in progress
*   Data Download Report does not show the "Data Level" section anymore.

### Bugs Fixed Since Last Release

*   The add to Cart button in File Entity Page changes its display value following a click. User can then click on “Remove from Cart”
*   Total Case value in the Cart matches with the number of cases associated with the files in the Cart
*   When User is authenticated, "My project flag" in Case table indicates that the Cases belongs to his projects
*   In Projects table, if User clicks on the count on Files, it links to the Data page - File table
*   In File Entity Page, if there are no associated cases, it will display the message "No Cases Found." instead of "No Participants Found."

### Known Issues and Workarounds

*   Missing pathology files and slide images for some TCGA datasets since they are not connected to the biospecimen chain yet
*   In some specific situation (complex queries using must_not term), the advanced search might returns incorrect results
*   When exporting the annotation table to a file, the create date field is not displayed in the same format than the UI
*   In the advanced search, adding a trailing space to the query causes the UI to display a syntax error although the query is valid
*   Checksum missing for MAGE-TAB files
*   Data download statistics report is not Section 508 compliant
*   User can download large files (> 1GB) from the Portal but it fails when extracting the tar.gz file. User should use the GDC Download Transfer Tool for large downloads
*   Authenticated user with access to TARGET projects cannot download Target DCC files. They will get a pop-up message.
*   In Case Entity Page, the "Download Clinical XML"button is not applicable to TARGET projects. The label is incorrect but the functionality is accurate.
*   The advanced search suggests more fields than expected (e.g. "cases.aliquots" should not be displayed). User should remove the field from the query and use the long field name (e.g. "cases.samples.portions.analytes.aliquots")
*   Removing a filter after paginating table results could cause no results to be displayed. User should press Clear on the filters and start again.
*   If a user has previously logged into the Portal and left a session without logging out, if he comes back to the Portal after his sessionID expires, it looks like the user is still authenticated. The user cannot download the token and gets an error message that would not close. The user should clear the cache to properly log out. 

Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).

## Release 0.2.15-oicr1

* __GDC Product__: GDC Data Portal
* __Release Date__: August 12, 2015

### New Features and Changes

*   Renamed all references of High Performance Download Client (HPDC) to GDC Data Transfer Tool 
*   Implemented more links between the Projects page and the Data page
*   Improved usability and visual experience:
    *   Better handling of login failure
    *   User feedback when no results are available following a search
    *   Warning for unsupported browsers (Internet Explorer 9 is not supported)

### Bugs Fixed Since Last Release

*   Optimizations for downloading large files from the browser
*   Fixed various issues related to file search
*   Fixed Section 508 compliance issues
*   Clicking on "Total case" link in the Project List Page does not return results.

### Known Issues and Workarounds

*   Missing pathology files and slide images for some TCGA datasets since they are not connected to the biospecimen chain yet
*   In some specific situation (complex queries using must_not term), the advanced search might returns incorrect results
*   The add to cart button in the files data view does change its display value following a click (file is correctly added to cart though)
*   When exporting the annotation table to a file, the create date field is not displayed in the same format than the UI
*   In the advanced search, adding a trailing space to the query causes the UI to display a syntax error although the query is valid
*   Checksum missing for MAGE-TAB files
*   Data download statistics report is not Section 508 compliant

Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).

## Release 0.2.13

* __GDC Product__: GDC Data Portal
* __Release Date__: July 23, 2015

### New Features and Changes

*   Improved data portal searching. Three search mechanisms are available to the user:
    *   **Facet search**: Starting with all content available on the GDC, allow users to filter down their search by clicking on elements on the left of the screen. This feature is available in **Projects** and **Data** page. Range support was also added to facets in this release
    *   **Advanced Search**: Starting with all content available on the GDC, allow users to build a custom and complex query using all of GDC capabilities (any field with the parameter of their choice).
    *   **Quick Search**: When looking for specific portal element, allow users to launch the **Quick Search** by clicking the "?" or "CRTL+SPACE" and find high-level informations of some entities.
*   Updated styling to align with NCI new visual identify
*   Created a pie chart widget allowing user to easily switch between a pie chart and a table view. 
*   Improved Usability and visual experience:
    *   Added tooltips to various sections of the portal
    *   Added a range facet with barchart
    *   Add more charts (summary, cart)

### Bugs Fixed Since Last Release

*   Hooked-up reports to real data
*   Fixed various issues on GQL (Advanced Search)
*   Table export to export appropriate columns
*   Allow users to sort project list table

### Known Issues and Workarounds

*   Checksum missing for MAGE-TAB files

Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).

## Release 0.1.10

* __GDC Product__: GDC Data Portal
* __Release Date__: March 18, 2015

### New Features and Changes

*   Authentication and Authorization
    *   Allow users to authenticate to the portal using their to ERA Commons credentials
    *   Allow users to narrow down their search based on their own projects
    *   Allow users to download a token link (for use with the API) in the portal
    *   Allow users to download protected data (if they have permissions to do so)
*   Reports
    *   Allow users to view a data download statistics report
*   Improve usability and visual experience:
    *   Allow users to view projects data using a new type of graph ("githut" style)
    *   Implement more features into the table widget (sort per column, hide/show columns, re-arrange columns, export to JSON, XML, CSV, TSV, maintain user preferences)
    *   Display UI and API version at the bottom of the page
*   Search
    *   Allow users to select multiple terms in facets (OR operand)
    *   Improvements on advance search with auto-loading of possible fields
*   Cart
    *   Allow users to view more details through additional pie charts
    *   Allow users to download a manifest from the cart
    *   Improved the mechanism of adding data to cart in various sections of the portal.
*   Updated style/theme to match GDC Website
*   Display a NCI Warning banner to inform users about GDC policy

### Bugs Fixed Since Last Release

*   Improvements in 508 compliance

### Known Issues and Workarounds

*   TARGET data is currently not available
*   Data:
    *   Absence of "real" download statistics for the data download statistics report
    *   Missing checksum for magetab files in file entity page
*   Project list table not sortable

Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).

## Release 0.1.8

* __GDC Product__: GDC Data Portal
* __Release Date__: January 22, 2015

### New Features and Changes

*   Allow users to perform a project search and obtain a list of projects (Project Search and List Pages)
*   Allow users to retrieve project details (Project Entity Page)
*   Allow users to perform an annotations search and obtain a list of annotations (Annotations Search and List Pages)
*   Allow users to retrieve entity details (Annotation Entity Page)
*   Implement a basic search feature (Basic Search - Facets))
*   Implement an advance search feature (Advanced Search - Query Language)
*   Allow users to retrieve file details (File Entity Page)
*   Allow users to retrieve participant details (Participant Entity Page)
*   Allow users to view and add files to a cart (Cart)
*   Allow users to view reports (Reports - Data Download Statistics)
*   Allow users to export tables of search results (Export Tables)
*   Allow users to download files (Download)
*   Allow users to authenticate using eRA Commons (Authentication)

### Bugs Fixed Since Last Release

*   Initial Release - Not Applicable

### Known Issues and Workarounds

*   TARGET data is currently not available

Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).
