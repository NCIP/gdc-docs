# Overview

While browsing through the GDC Data Portal, files can be downloaded either individually or bundled through a cart mechanism similar to online stores.

The cart pages shows a summary of all files currently in the cart, display a file count by project and identified the authorization level of those files (controlled or open access).

# GDC Cart

ADD MORE DETAILS

# Download Options

From the cart, three different download options are available to the end user.

__Note__: when downloading multiple files from the GDC Data Portal, those files are automatically bundled-up into one single Gzipped (.tar.gz) file.

## GDC Data Transfer Tool

The "Manifest for Fast Download Tool" button will download a manifest that can be imported into the GDC Data Transfer Tool.

```
id	filename	md5	size	state
4ea9c657-8f85-44d0-9a77-ad59cced8973	mdanderson.org_ESCA.MDA_RPPA_Core.mage-tab.1.1.0.tar.gz	\N	2516051	live
b8342cd5-330e-440b-b53a-1112341d87db	mdanderson.org_SARC.MDA_RPPA_Core.mage-tab.1.1.0.tar.gz	\N	4523632	live
c57673ac-998a-4a50-a12b-4cac5dc3b72e	mdanderson.org_KIRP.MDA_RPPA_Core.mage-tab.1.2.0.tar.gz	\N	4195746	live
3f22dd8d-59c8-43a4-89cf-3b595f2e5a06	14-3-3_beta-R-V_GBL1112940.tif	56df0e4b4fc092fc3643bd2e316ac05b	6257840	live
7ce05059-9197-4d38-830f-04356f5f851a	14-3-3_beta-R-V_GBL11066140.tif	6abfee483974bc2e61a37b5499ae9a07	6261580	live
8e00d22a-ca6f-4da8-a1c3-f23144cb21b7	14-3-3_beta-R-V_GBL1112940.tif	56df0e4b4fc092fc3643bd2e316ac05b	6257840	live
96487cd7-8fa8-4bee-9863-17004a70b2e9	14-3-3_beta-R-V_GBL1112940.tif	56df0e4b4fc092fc3643bd2e316ac05b	6257840	live
```

The Manifest contains a list of the file UUIDs in the cart and can be used together with the GDC Data Transfer Tool to download all files.

Information on the GDC Data Transfer Tool is available in the [GDC Data Transfer Tool User’s Guide](/node/8196/).

## Cart Download

The "Download Cart (web browser)" button will initial a data download directly from the browser. Users have to be cautious of the amount of data in the cart since this option will not optimize bandwidth neither provide resume capabilities.

## Individual Files Download

Similar to the files page, each row contains a download button to download a particular file individually.

# Access

If a user tries to download a cart containing controlled data and without being authenticated, a pop-up will be displayed to offer the user either to download only open access files or to login into the GDC Data Portal through eRA Commons.

[![Cart Page](images/gdc-data-portal-download-cart.png)](images/gdc-data-portal-download-cart.png "Click to see the full image.")
