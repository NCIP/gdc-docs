#Troubleshooting Guide

If you encounter issues when using the Data Transfer Tool for downloading files please reference the section below for helpful hints and recommendations.  

##Speed Performance During Download
The Data Transfer Tool has two performance tuning options that are presented during download operations.  The two options are:

###--n  

The "--n" option assists with assigning the number of threads to the download process.  The default is 4 and can not be lowered below three threads.

###--http-chunk-size

The "--http-chunk-size" setting can improve performance but we do not provide any hard settings due to the eclectic nature of client networks and their connections to the internet but instead encourage clients to experiment with changing the default setting of 1048576 bytes to larger size ranges.       

##Very Large Manifests
Some clients have needed to create very large manifest files to satisfy the scope of their work.  Using very large manifest files can from time to time lead to the end user experiencing  network time outs or dropped connections due to network topologies, internet connections, or load on client side systems.  Network time out or dropped network connect can manifest as a hung or unresponsive download session. To help mitigate these issues we recommended that clients breakup their manifest files into smaller chunks.  

##General Tips
To avoid running into older software bugs/conflicts we recommend you always use the latest version of the client whenever possible.  When experiencing download problems while using an access tokens try downloading a new token first before reporting it to the GDC helpdesk.  

##Logging
From time to time the GDC User Services Team will request that you run the command line application with the following flags { --debug --logfile } to assist with troubleshooting any issues that might have appeared.  These flags will run the application in debug mode and create a logfile file with the debug logs in it.  
Example Usage:
```Debug-Logfile
gdc-client download -m lung.manifest.txt -t token.file --debug --logfile logfile.txt
```

##OS Compatibility with the Data Transfer Tool
The Data transfer Tool is offered in three OS compatible versions; Mac OS, Windows, and Ubuntu Linux.  We have successfully tested the Ubuntu binary on CentOS 7.x and RHEL 7 and Scientific Linux 7 with the client but have had problems with CentOS 6.x and RHEL 6 and SL6.  To work around this problem we have asked users to build their own client from our [github](https://github.com/NCI-GDC/gdc-client) repository with the assistance of an instruction document that we provide on request via the GDC Helpdesk.        


##Network Troubleshooting
Network problems can appears as dropped network connections or even a stalled application.  The GDC Helpdesk might request more network information to assist in diagnosing the problem.  The two tests they will request the end user to run are ping and traceroute (tracert on the windows plateform) against our api servers.  Please capture the output from these tests into a text file and attach it to the reply email.   
Examples:
```Ping
>ping api.gdc.cancer.gov
 PING api.gdc.cancer.gov (192.170.230.246): 56 data bytes
 64 bytes from 192.170.230.246: icmp_seq=0 ttl=249 time=4.235 ms
 64 bytes from 192.170.230.246: icmp_seq=1 ttl=249 time=4.783 ms
```
```Traceroute
>traceroute api.gdc.cancer.gov
 traceroute to api.gdc.cancer.gov (192.170.230.246), 64 hops max, 52 byte packets
 1  h01-391-250-v1011.gw.uchicago.net (10.151.0.2)  4.595 ms  3.602 ms  3.322 ms
 2  h01-391-250-to-b65-ll129-300.p2p.uchicago.net (10.5.1.32)  13.285 ms  9.241 ms  5.156 ms
 3  b65-ll129-300-to-borderfw.p2p.uchicago.net (192.170.192.32)  3.218 ms  3.364 ms  3.396 ms
 4  borderfw-to-b65-ll129-500.p2p.uchicago.net (192.170.192.36)  3.605 ms  3.741 ms  3.833 ms
 5  b65-scidmz-01-to-b65-ll129-500.uchicago.net (128.135.247.182)  4.223 ms  5.428 ms  3.999 ms
 6  192.170.224.97 (192.170.224.97)  4.003 ms  3.970 ms  6.260 ms
 7  lnk-g30-scidist-01.scidmz.uchicago.net (192.170.224.66)  4.530 ms  4.649 ms  6.021 ms
 8  192.170.230.246 (192.170.230.246)  4.158 ms  4.273 ms  5.134 ms
```

##Common Error Codes
This is a list of the most common error codes the DTT generates and their meaning
<ul TYPE="square">
<li> Unable to connect to API – you might be running an out of data client so consider upgrading.</li>
<li> Error: Max Retries Exceeded – network connect timeouts </li>
<li>CryptographyDeprecationWarning – a warning that you should consider upgrading to a higher   
  version of python – please upgrade to 2.7.x or higher.</li>
<li>ERROR: An unexpected error has occurred during normal operation of the client -	This could  be a variety of problems and we ask you to contact our helpdesk.</li>
<li>ECONNRESET - network connection dropped.</li>
