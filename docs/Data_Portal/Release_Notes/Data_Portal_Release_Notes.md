# Data Portal Release Notes

## Release 1.1.0

* __GDC Product__: GDC Data Portal
* __Release Date__: June 1st, 2016

### New Features and Changes

*   This is a bug-fixing release, no new features were added.

### Bugs Fixed Since Last Release

*   General
    *  '--' symbol is missing on empty fields (blank instead) <!-- PGDC-2418 -->
    *  Fixed 508 compliance issues <!-- PGDC-2497 --><!-- PGDC-2431 -->
    *  Download manifest button is available for TARGET projects with 0 files, resulting in error if user clic on button <!-- PGDC-2416 -->
    *  If a user tries to download a token after session expiry, an invite will offer the user to log in <!-- PGDC-2455 -->   
    *  My project filter should be unselected by default after user logs-in <!-- PGDC-2462 -->  
    *  Fixed missing padding when query includes "My Projects" <!-- PGDC-2420 -->
    *  The portal was allowing download manifest when there was no files in the project <!-- PGDC-2416 -->
    *  Add to cart limitation to 10,000 files was not enforced everywhere on the Data Portal <!-- PGDC-2409 -->
*   Data
    *  Tables: Improved usability of the "Sort" feature  <!-- PGDC-1771 -->
    *  Tables: adding all files to the cart was only adding files displayed on the current page instead of all files corresponding to the current query <!-- PGDC-2439 -->
    *  Tables: Platform show "0" when selected platform is "Affymetrix SNP 6.0" <!-- PGDC-2419 -->
    *  When adding a custom range facet, default values are incorrectly populated <!-- PGDC-2445 --> <!-- PGDC-2444 -->  <!-- PGDC-2225 -->    
    *  Sorting File Submitter ID option on the file tab result in a Data Portal Error <!-- PGDC-2430 -->    
*   File Entity Page
    *  Improved "Associated Cases/Biospecimen" table for files associated to a lot of cases  <!-- PGDC-1993 -->    
    *  BAM Slicing was retuning an error <!-- PGDC-2433 -->


### Known Issues and Workarounds

*   General
    *   If a user has previously logged into the Portal and left a session without logging out, if the user returns to the Portal after the user's sessionID expires, it looks as if the user is still authenticated. User needs to refresh the page to see current authentication status.
    *   After successful authentication, the authentication popup does not close for Internet Explorer users running in "Compatibility View". This only impact users at the NIH. Workaround is to uncheck "Display Intranet sites in Compatibility View" in Internet Explorer options. Alternatively, refreshing the portal will correctly display authentication status.
    *   Download SRA XML files from the legacy archive portal might not be possible in some context <!-- PGDC-2457 --> <!-- PGDC-2469 -->  
    *   Not possible to use the browser's back button after hitting a 404 page <!-- PGDC-2429 -->
    *   404 page missing from Legacy Archive Portal <!-- PGDC-2477 -->
    *   Table widget icon and export JSON icon should be different <!-- PGDC-2446 -->     
*   Data and facets
    *   Default values for age at diagnosis is showing 0 to 89 instead of 0 to 90 <!-- PGDC-2478 -->
    *   Biospecimen search in the case entity page does not highlight (but does bold and filter) results in yellow when title case is not followed <!-- PGDC-2451 -->
    *   Table sorting icon does not include numbers <!-- PGDC-35 -->    
    *   '--' symbol is missing on empty fields (blank instead), additional missing fields identified.  <!-- PGDC-2447 -->    
*   Cart
    *   When removing a file from the cart and clicking undo, GDC looses track of permission status of the user towards this file and will ask for the user to log-in again. <!-- PGDC-2496 -->     
    *   Counts displayed in the top right of the screen, next to the Cart icon, might get inconsistent if files are removed from the cart. <!-- PGDC-2403 -->
*   Annotations
    *   Annotations notes do not wrap to the next line at the beginning or the end of a word, some words might be split in two lines <!-- PGDC-2474 -->
*   Web Browsers
    *   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.nci.nih.gov/about-gdc/gdc-faqs).
    *   Internet Explorer users are not able to use the "Only show fields with no values" when adding custom facets <!-- PGDC-2467 -->
    *   The GDC Portals are not compatible with Internet Explorer running in compatibility mode. Workaround is to disable compatibilty mode <!-- PGDC-2480 -->    

Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).


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
    *   If a user has previously logged into the Portal and left a session without logging out, if the user returns to the Portal after the user's sessionID expires, it looks as if the user is still authenticated. The user cannot download the token and gets an error message that would not close. The user should clear the cache to properly log out.
    *   '--' symbol is missing on empty fields (blank instead) <!-- PGDC-2418 -->
    *   Download manifest button is available for TARGET projects with 0 files, resulting in error if user clic on button <!-- PGDC-2416 -->
    *   After successful authentication, the authentication popup does not close for Internet Explorer users running in "Compatibility View". This only impact users at the NIH. Workaround is to uncheck "Display Intranet sites in Compatibility View" in Internet Explorer options. Alternatively, refreshing the portal will correctly display authentication status.
*   Data    
    *   When adding a custom range facet, default values are incorrectly populated <!-- PGDC-2445 --> <!-- PGDC-2444 -->  <!-- PGDC-2225 -->
    *   The portal might return incorrect match between cases and files when using field cases.samples.portions.created_datetime (custom facet or advanced search). Note: this is not a UI issue. <!-- PGDC-2440 -->
    *   Sorting File Submitter ID option on the file tab result in a Data Portal Error <!-- PGDC-2430 -->
*   Tables and Export
    *   Table sorting icon does not include numbers <!-- PGDC-35 -->    
*   Browsers limit the number of concurrent downloads, it is generally recommended to add files to the cart and download large number of files through the GDC Data Transfer Tool, more details can be found on [GDC Website](https://gdc.nci.nih.gov/about-gdc/gdc-faqs).

Release details are maintained in the [GDC Data Portal Change Log](https://github.com/NCI-GDC/portal-ui/blob/master/CHANGELOG.md).
