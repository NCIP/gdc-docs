# Data Portal Release Notes


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
    *   After successful authentication, the authentication popup does not close for Internet Explorer users running in "Compatibility View". Workaround is to uncheck "Display Intranet sites in Compatibility View" in Internet Explorer options. Alternatively, refreshing the portal will correctly display authentication status. <!-- PGDC-2403 / PRTL-133 -->
    *   BAM Slicing dialog box does not disappear automatically upon executing the BAM slicing function. The box can be closed manually. <!-- PRTL-282 -->
    *   Due to preceding issue, If bam slicing produces an error pop-up message it will be obscured behind the original dialog box. <!--SV-419-->
    *   Very long URLs will produce a 400 error.  Users may encounter this after clicking on "source files" on a file page where the target file is derived from hundreds of other files such as for MAF files.  To produce a list of source files an API call can be used with the search parameter "fields=analysis.input_files.file_name". <!-- SV-396 / PRTL-342-->
		*   Downloading a token in the GDC Legacy Archive does not refresh it. If a user downloads a token in the GDC Data Portal and then attempts to download a token in the GDC Legacy Archive, an old token may be provided. Reloading the Legacy Archive view will allow the user to download the updated token.


Example

    https://gdc-api.nci.nih.gov/files/455e26f7-03f2-46f7-9e7a-9c51ac322461?pretty=true&fields=analysis.input_files.file_name




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

    https://gdc-api.nci.nih.gov/files/455e26f7-03f2-46f7-9e7a-9c51ac322461?pretty=true&fields=analysis.input_files.file_name




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
