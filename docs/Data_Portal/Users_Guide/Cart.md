# Cart and File Download

## Overview

While browsing through the GDC Data Portal, files can be downloaded either individually or bundled through a cart mechanism similar to online stores.


## GDC Cart

[![Cart](images/gdc-data-portal-cart.png)](images/gdc-data-portal-cart.png "Click to see the full image.")

### Cart Summary

The cart page shows a summary of all files currently in the cart:

* Number of files
* Number of cases associated with the files
* Total file size

The Cart page also displays two charts: 

* File count by project 
* File count by authorization level to those files: this number takes into account the user authorization to download a file upon his authentication to the GDC Data Portal.

### Cart table

The Cart table shows the list of all the files that were added to the Cart.

The user can click on the links to get more information on each file, it will link to:

* The file entity page
* The case entity page
* The project entity page


## Download Options

From the cart, the fowllowing download options are available to the end user:

* __Manifest__: Manifest used by the GDC Data Transfer Tool to download the files.
* __Cart__: All the files in the Cart downloaded directly through the browser. Users have to be cautious of the amount of data in the cart since this option will not optimize bandwidth neither provide resume capabilities.
* __Clinical__: GDC harmonized clinical data associated with the cases in the cart.
* __Biospecimen__: GDC harmonized biospecimen data associated with the cases in the cart.
* __SRA XML__: Experiment, analysis and run metadata files associated with the files in the cart
* __File Metadata__: Metadata of the files in the cart (file properties, associated entities and workflow information if applicable).

Although an entire cart could be downloaded from Web Browsers, those are not ideally equipped to download very large files, in particular due to the absence of a retry/resume mechanism. To improve user experience we implemented a limit of 5 GB for direct cart download. Over 5 GB we recommend using the GDC Data Transfer Tool.

__Note__: when downloading multiple files from the GDC Data Portal, those files are automatically bundled-up into one single Gzipped (.tar.gz) file.

### GDC Data Transfer Tool

The "Download Manifest" button will download a manifest that can be imported into the GDC Data Transfer Tool.

```manifest
id	filename	md5	size	state
4ea9c657-8f85-44d0-9a77-ad59cced8973	mdanderson.org_ESCA.MDA_RPPA_Core.mage-tab.1.1.0.tar.gz		2516051	live
b8342cd5-330e-440b-b53a-1112341d87db	mdanderson.org_SARC.MDA_RPPA_Core.mage-tab.1.1.0.tar.gz		4523632	live
c57673ac-998a-4a50-a12b-4cac5dc3b72e	mdanderson.org_KIRP.MDA_RPPA_Core.mage-tab.1.2.0.tar.gz		4195746	live
3f22dd8d-59c8-43a4-89cf-3b595f2e5a06	14-3-3_beta-R-V_GBL1112940.tif	56df0e4b4fc092fc3643bd2e316ac05b	6257840	live
7ce05059-9197-4d38-830f-04356f5f851a	14-3-3_beta-R-V_GBL11066140.tif	6abfee483974bc2e61a37b5499ae9a07	6261580	live
8e00d22a-ca6f-4da8-a1c3-f23144cb21b7	14-3-3_beta-R-V_GBL1112940.tif	56df0e4b4fc092fc3643bd2e316ac05b	6257840	live
96487cd7-8fa8-4bee-9863-17004a70b2e9	14-3-3_beta-R-V_GBL1112940.tif	56df0e4b4fc092fc3643bd2e316ac05b	6257840	live
```

The Manifest contains a list of the file UUIDs in the cart and can be used together with the GDC Data Transfer Tool to download all files.

Information on the GDC Data Transfer Tool is available in the [GDC Data Transfer Tool User's Guide](/node/8196/).


### Individual Files Download

Similar to the files page, each row contains a download button to download a particular file individually.

## Controlled Files

If a user tries to download a cart containing controlled files and without being authenticated, a pop-up will be displayed to offer the user either to download only open access files or to login into the GDC Data Portal through eRA Commons.

[![Cart Page](images/gdc-data-portal-download-cart.png)](images/gdc-data-portal-download-cart.png "Click to see the full image.")
