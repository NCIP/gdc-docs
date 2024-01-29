# Repository

## Introduction ##

The Repository tool is where the files associated with each case in the current cohort can be browsed and downloaded. It also offers users a variety of file filters for identifying files of interest.

The Repository tool can be reached in one of these two ways:
* choosing the Repository link in the GDC Data Portal header
* clicking the play button on the Repository card in the Analysis Center

[![Repository In Header](images/ToolLinksInHeader.png)](images/ToolLinksInHeader.png "Click to see the full image.")

## Choosing a Cohort ##

When searching for files to download, many users will have in mind a specific cohort whose associated files they wish to download. The set of files displayed in the Repository at any time will reflect the files that are associated with the active cohort. The current active cohort can be seen in the Main Toolbar (namely at the top of the page in the Analysis Center):

[![Main Cohort Toolbar](images/MainCohortToolbar.png)](images/MainCohortToolbar.png "Click to see the full image.")

The active cohort can be changed by clicking the cohort name and choosing a new one from the dropdown menu.  If the cohort of interest does not appear in the list, try [creating a new cohort.](quick_start.md#main-toolbar)

For users who want to browse all files that are available at the GDC, create a new cohort via the main toolbar and use it with the Repository tool.

## Filtering a Set of Files ##

Users may only be interested in browsing through or downloading a subset of files associated with the current cohort. For this purpose, a set of commonly-used default facet cards is provided in the left panel of the Repository tool to allow users to filter the files presented in the table on the right. The facet cards are as follow:

* **Experimental Strategy**: Experimental strategies used for molecular characterization of the cancer
* **WGS Coverage**: Range of coverage for WGS aligned reads
* **Data Category**: A high-level data file category, such as "Raw Sequencing Data" or "Transcriptome Profiling"
* **Data Type**: Data file type, such as "Aligned Reads" or "Gene Expression Quantification". Data Type is more granular than Data Category.
* **Data Format**: Format of the data file
* **Workflow Type**: Bioinformatics workflow used to generate or harmonize the data file
* **Platform**: Technological platform on which experimental data was produced
* **Access**: Indicator of whether access to the data file is open or controlled

 Values within each facet can be sorted alphabetically by choosing the "Name" header on the top left of each card. Alternatively, the "Files" header may be selected to sort the values by the number of files available.

Note that the categories displayed in the filters represent the values available for the active cohort.

[![Full Repository](images/FullRepo.png)](images/FullRepo.png "Click to see the full image.")

If a different filter needs to be used, a custom filter can be applied by choosing the "Add a Custom Filter" button at the top of the default filters. Each custom filter can then be searched and chosen within the pop-up window. Once a custom filter is selected, a new filter card will appear at the top of the default filters.  Custom filters can be removed from the Repository by choosing the X at the top right of each filter card.

[![Custom File Filter](images/CustomFileFilter.png)](images/CustomFileFilter.png "Click to see the full image.")

## Viewing Images ##

To view images associated with the active cohort, select the View Images button above the files table to launch the Slide Image Viewer. 

## Files Table

The table shows the list of all the files associated with the active cohort, subject to any filtering that has been applied in the Repository. By default, the table provides the following information for each file:

* **Access**: Displays whether the file is open or controlled access. Users must login to the GDC Portal and have the appropriate credentials to access these files.
* **File Name**: Name of the file. Clicking the link will bring the user to the File Summary Page.
* **Cases**: The number of cases associated with the file
* **Project**: The Project that the file belongs to. Clicking the link will bring the user to the Project Summary Page.
* **Data Category**: Type of data
* **Data Format**: The file format
* **File Size**: The size of the file
* **Annotations**: Whether there are any annotations

Additional information such as Data Type and Experimental Strategy can be displayed using the Customize Columms button above the table. The table can be sorted by clicking on the headers, and the search bar above the table can be used to locate specific files.

The JSON / TSV buttons will download the files' details (file name, file size, data category, access type, etc.) in JSON and TSV format, respectively.

## Downloading a Set of Files ##

When filtering has been completed, files are ready to be downloaded. Depending on the number and size of files, the GDC has several options and recommendations for downloading them. While any amount of data can be downloaded using the GDC Data Transfer Tool or the API, files can be downloaded directly from the Data Portal if the size is 5 GB or less in total and the number of files does not exceed 10,000. For any downloads larger than 5 GB or 10,000 files, it's recommended that the download be performed using the [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool). 

### Generating a Manifest File for the Data Transfer Tool ###

Select the Manifest button above the table to generate a manifest file required for batch download using the Data Transfer Tool. The manifest contains a list of the UUIDs corresponding to the files associated with the active cohort, subject to any filtering in the Repository. 

### Adding/Removing Files to the Cart for Download ###

Downloads can also be performed using the Cart by first adding a set of files to the Cart. This can be done using the following methods:
* Clicking the cart icon on the left of each file. This will toggle between adding to and removing the file from the cart.
* Selecting the Add All Files to Cart button. This will add all the files in the current cohort to the Cart, subject to any filtering that has been applied in the Repository.

[![Add Files To Cart](images/AddFilesToCart.png)](images/AddFilesToCart.png "Click to see the full image.")

### Cart

The Cart page can then be reached by clicking the Cart icon at the top right of the portal.

At the upper-right of the page is a summary of all files currently in the cart:
* Number of files
* Number of cases associated with the files
* Total file size

The Cart page displays the file count by project and authorization level, as well as a table of all files that have been added to the Cart.  Files can be removed from the Cart using the trash icons at the left of each file in the table or by selecting the "Remove from Cart" option at the top of the Cart page, which removes either all files or the unauthorized ones.

[![Cart Page](images/CartPage.png)](images/CartPage.png "Click to see the full image.")

### Cart Items Table

The Cart Items table shows the list of all the files that were added to the Cart and has the same functionality as the table in the Repository. By default, it displays the following information for each file:

* **Access**: Displays whether the file is open or controlled access. Users must login to the GDC Portal and have the appropriate credentials to access these files.
* **File Name**: Name of the file. Clicking the link will bring the user to the File Summary Page.
* **Cases**: The number of cases associated with the file
* **Project**: The Project that the file belongs to. Clicking the link will bring the user to the Project Summary Page.
* **Data Category**: Type of data
* **Data Format**: The file format
* **File Size**: The size of the file
* **Annotations**: Whether there are any annotations


Additional information can be displayed using the Customize Columms button above the table. Sort can be applied by clicking on the table headers, and the search bar provides additional options for locating specific files. Details of the files can be downloaded using the JSON and TSV buttons above the table.


### Downloading Files from the Cart

To download files in the Cart, select the Download Cart button and choose either:
* **Manifest**: Downloads a manifest for the files that can be passed to the GDC Data Transfer Tool. A manifest file contains a list of the UUIDs that correspond to the files in the cart.
* **Cart**: Download the files directly through the browser. Users have to be cautious of the amount of data in the cart since this option will not optimize bandwidth and will not provide resume capabilities. This option can only be used if the total size of the files in the Cart does not exceed 5 GB.

### Additional Data Download

Additional data can be downloaded from the Cart page using the Download Associated Data button at the top of the page and choosing one of the available options.

[![Cart Metadata](images/CartAssociatedData.png)](images/CartAssociatedData.png "Click to see the full image.")
* Clinical: TSV / Clinical: JSON - This includes all clinical information from the cases that are associated with the files (available as TSV or JSON)
* Biospecimen: TSV / Biospecimen: JSON - This includes all biospecimen information from the cases that are associated with the files (available as TSV or JSON).
* Sample Sheet -  A TSV with commonly-used elements associated with each file, such as sample barcode and sample type.
* Metadata - This includes all of the metadata associated with each and every file in the cart.  Note that this file is only available in JSON format and may take several minutes to download.

## File Summary Page ##

Users can launch the File Summary Page by clicking on a file name in the tables that appear on both the Repository and Cart pages. The File Summary Page provides information about a data file, including file properties like size, MD5 checksum, and data format; information on the type of data included; links to the associated cases and biospecimen; and information about how the data file was generated or processed.

The page also includes buttons to download the file, add it to the file cart, or (for BAM files) utilize the BAM slicing function.

[![File Summary Page](images/FileSummaryPage.png)](images/FileSummaryPage.png "Click to see the full image.")

In the lower section of the screen, the following tables provide more details about the file and its characteristics:

* __Associated Cases/Biospecimen__: List of cases or biospecimen the file is directly attached to
* __Analysis and Reference Genome__: Information on the workflow and reference genome used for file generation
* __Read Groups__: Information on the read groups associated with the file
* __Metadata Files__: Experiment metadata, run metadata and analysis metadata associated with the file
* __Downstream Analysis Files__: List of downstream analysis files generated by the file
* __File Versions__: List of all versions of the file
