# Getting Started

## The GDC Data Transfer Tool: An Overview

Raw sequence data, typically stored as BAM or FASTQ files, make up the bulk of data stored at The NCI Genomic Data Commons (GDC) . The size of a single file can vary greatly. Most BAM files stored in the GDC are in the 50 MB - 40 GB size range, with some of the whole genome BAM files reaching sizes of 200-300 GB.

The GDC Data Transfer Tool provides an optimized method of transferring data to and from the GDC, and enables resumption of interrupted transfers. The tool's download and upload functionality is built upon a high-performance, UDP-based data transfer protocol known as [UDT](http://udt.sourceforge.net/), or User Datagram Protocol (UDP)-based Data Transfer. The UDT protocol is implemented in C++ for speed, but exposed via Python to provide a clean, user-friendly command line interface. The GDC Data Transfer Tool is also capable of performing TCP downloads for cases where UDP traffic may be blocked.

## Downloading the GDC Data Transfer Tool

### System Recommendations

The system recommendations for using the GDC Data Transfer Tool&nbsp;are as follows:

* OS: Linux (Ubuntu 14.x or later), OS X (10.9 Mavericks or later), or Windows (7 or later)
* CPU: At least eight 64-bit cores, Intel or AMD
* RAM: At least 8 GiB
* Storage: Enterprise-class storage system capable of at least 1 Gb/s (gigabit per second) write throughput and sufficient free space for BAM files.

### Binary Distributions

Binary distributions are available on the [GDC Transfer Tool page](https://gdc.nci.nih.gov/access-data/gdc-data-transfer-tool). To install the GDC Data Transfer Tool, download the respective binary distribution and unzip the distribution&#39;s archive to a location on the target system.

### Release Notes

Release Notes are available on the [GDC Data Transfer Tool Release Notes](../../Data_Transfer_Tool/Release_Notes/DTT_Release_Notes.md) page.
