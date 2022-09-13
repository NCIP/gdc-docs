# Getting Started

## The GDC Data Transfer Tool: An Overview

Raw sequence data, stored as BAM files, make up the bulk of data stored at the NCI Genomic Data Commons (GDC). The size of a single file can vary greatly. Most BAM files stored in the GDC are in the 50 MB - 40 GB size range, with some of the whole genome BAM files reaching sizes of 200-300 GB.

The GDC Data Transfer Tool, a command-line driven application, provides an optimized method of transferring data to and from the GDC and enables resumption of interrupted transfers.

## Downloading the GDC Data Transfer Tool

### System Recommendations

The system recommendations for using the GDC Data Transfer Tool are as follows:

* OS: Linux (Ubuntu 16.x or later), OS X (10.9 Mavericks or later), or Windows (8 or later)
* CPU: At least two 64-bit cores, Intel or AMD
* RAM: At least 8 GiB
* Storage: Enterprise-class storage system capable of at least 1 Gb/s (gigabit per second) write throughput and sufficient free space for BAM files.

### Binary Distributions

Binary distributions are available on the [GDC Transfer Tool page](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool). To install the GDC Data Transfer Tool, download the respective binary distribution and unzip the distribution's archive to a location on the target system. It is recommended that the binary be copied to a located that is in the user's path so that is it accessible from any location within the terminal or command prompt.   

### Release Notes

Release Notes are available on the [GDC Data Transfer Tool Release Notes](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Release_Notes/DTT_Release_Notes) Page.
