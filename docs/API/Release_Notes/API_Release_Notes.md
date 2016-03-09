# API Release Notes






## v0.3.24.1

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: March 1, 2016

### New Features and Changes

* New entities: *Slide Data Bundle* and *Pathology Data Bundle*
* New `_dictionary/_all` endpoint displays the entire GDC dictionary for the `submission` endpoint, and project-specific dictionaries for project endpoints; e.g. `https://gdc-api.nci.nih.gov/v0/submission/_dictionary/_all` and `https://gdc-api.nci.nih.gov/v0/submission/TCGA/SKCM/_dictionary/_all`.
* Support for releases
* Project `state` property definition has been updated. The property can now have one of six values: `open` (default), `review`, `submitted`, `processing`, `closed`, and `legacy`.
* `days_to_index` property has been removed from `Case`
* `template` endpoint for submissions returns archive of all templates, e.g. `https://gdc-api.nci.nih.gov/v0/submission/TCGA/SKCM/template`
* `template` endpoint can generate JSON using `/template?format=json`
* `export` endpoint now supports POST requests
* GraphQL: lists are accepted as top level property filters
* GraphQL: `not` filter was added

### Bugs Fixed Since Last Release

* Portion is now linked to Center

### Known Issues and Workarounds

* None to report

Release details are maintained in the [GDC API change log](https://github.com/NCI-GDC/gdcapi/blob/develop/CHANGELOG.md)






## v0.3.20.1

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: January 28, 2016

### New Features and Changes

* Metadata export support for TSV, CSV and JSON formats, including:
  *  Clinical metadata
  *  Biospecimen metadata
  *  Project and case entities
  *  Annotations
* Export of metadata can be filtered by category
* Support for generating a TSV template
* Additional details are provided in JSON parsing error messages

### Bugs Fixed Since Last Release

* Biospecimen XML submissions with shipment_portion nodes now parsed correctly
* Correct HTTP code returned upon authentication failure
* API URLs work with and without trailing slash
* Miscellaneous bug fixes

### Known Issues and Workarounds

*   None to report

Release details are maintained in the [GDC API change log](https://github.com/NCI-GDC/gdcapi/blob/develop/CHANGELOG.md)






## v0.2.18

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: November 18, 2015

### New Features and Changes

*   Added support for molecular data upload.
  *   Users can submit data bundle entities to specify the molecular data and relations to be uploaded.
  *   Users can update existing data by submitting update data bundles.
  *   Manifest files can be downloaded to allow for molecular data upload via the API or Data Transfer Tool.
*   Users can delete previously submitted entities.

### Bugs Fixed Since Last Release

*   None to report

### Known Issues and Workarounds

*   None to report

Release details are maintained in the [GDC API change log](https://github.com/NCI-GDC/gdcapi/blob/develop/CHANGELOG.md)






## v0.2.18-spr1

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: November 2, 2015

### New Features and Changes

*   Entity Metadata submission added with support for JSON, CSV, TSV, and XML formats and the following data entities:
    *   Clinical metadata
    *   Biospecimen metadata
    *   Project and case entities
    *   Annotations
*   Dry run option available for metadata submission to validate metadata files prior to submission
*   Submission integrated with authorization to ensure read, create, update, and admin roles can be controlled at the project level

### Bugs Fixed Since Last Release

*   None to report

### Known Issues and Workarounds

*   None to report

Release details are maintained in the [GDC API change log](https://github.com/NCI-GDC/gdcapi/blob/develop/CHANGELOG.md)






## v0.2.15-oicr1

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: August 7, 2015


### New Features and Changes

*   None to report

### Bugs Fixed Since Last Release

*   Addressed an issue in the range aggregation feature (used by facets)

### Known Issues and Workarounds

*   Query and Download Service
    *   Universal Resource Locator (URL) length is limited.  As such, download requests with many files have to be divided into multiple requests. A feature that allows a payload to specify the file identifiers to download to overcome this limitation is planned for a future release.
    *   All TARGET, CGCI, and CCLE data is currently not available via the GDC API but will be made available as data is continously imported into the GDC
*   Authentication and Authorization Service
    *   Obtaining a token requires browser access. This is a limitation of eRA Commons and not the GDC authorization service. The Security Assertion Markup Language  (SAML) based protocol eRA Commons uses does not support [ECP](https://wiki.shibboleth.net/confluence/display/CONCEPT/ECP) (or similar)

Release details are maintained in the [GDC API change log](https://github.com/NCI-GDC/gdcapi/blob/develop/CHANGELOG.md)






## v0.2.13

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: July 23, 2015


### New Features and Changes

*   Query and Download Service
    *   Implemented server-side related files
    *   Created a quick search endpoint
*   Misc
    *   Added POST support
    *   Added range support (used by range facet)
    *   Hooked-up reports to new reports index

### Bugs Fixed Since Last Release

*   Query and Download Service
    *   Fixed error (500) on annotations without notes
    *   Correct structure for related files
    *   Fixed error (404) with annotation download
    *   Fixed issue with support for IS NOT operator

### Known Issues and Workarounds

*   Query and Download Service
    *   Universal Resource Locator (URL) length is limited.  As such, download requests with many files have to be ivided into multiple requests. A feature that allows a payload to specify the file identifiers to download to zovercome this limitation is planned for a future release.
    *   All TARGET, CGCI, and CCLE data is currently not available via the GDC API but will be made available as data is continously imported into the GDC
*   Authentication and Authorization Service
    *   Obtaining a token requires browser access. This is a limitation of eRA Commons and not the GDC authorization service. The Security Assertion Markup Language  (SAML) based protocol eRA Commons uses does not support [ECP](https://wiki.shibboleth.net/confluence/display/CONCEPT/ECP) (or similar)

Release details are maintained in the [GDC API change log](https://github.com/NCI-GDC/gdcapi/blob/develop/CHANGELOG.md)






## v0.1.10

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: March 18, 2015


### New Features and Changes

*   Query and Download Service
    *   Supports query for files, participants, annotations, and projects by metadata. Interfaces with ElasticSearch for efficient queries.
    *   Provides a RESTful interface to support the download of data, metadata, and annotations
    *   Supports the download of multiple files, manifests, or annotations in compressed (gzip) format by default
*   Authentication and Authorization Service
    *   Authenticates users using the GDC authorization service which authorizes against dbGaP
    *   Manages tokens and token-based authorization using OpenStack Keystone
    *   Implements federated authentication with eRA Commons to generate a Keystone token
    *   Downloads Comma Separated Values (CSVs) from dbGaP at regular intervals and verifies authorization using eRA Commons credentials

### Bugs Fixed Since Last Release

*   Query and Download Service
    *   Initial Release - Not Applicable
*   Authentication and Authorization Service
    *   Further developed the prototype into a production ready service

### Known Issues and Workarounds

*   Query and Download Service
    *   Universal Resource Locator (URL) length is limited.  As such, download requests with many files have to be divided into multiple requests. A feature that allows a payload to specify the file identifiers to download to overcome this limitation is planned for a future release.
    *   All TARGET, CGCI, and CCLE data is currently not available via the GDC API but will be made available as data is continously imported into the GDC
*   Authentication and Authorization Service
    *   Obtaining a token requires browser access. This is a limitation of eRA Commons and not the GDC authorization service. The Security Assertion Markup Language  (SAML) based protocol eRA Commons uses does not support [ECP](https://wiki.shibboleth.net/confluence/display/CONCEPT/ECP) (or similar)

Release details are maintained in the [GDC API change log](https://github.com/NCI-GDC/gdcapi/blob/develop/CHANGELOG.md)






## v0.1.8

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: January 22, 2015


### New Features and Changes

*   Identity Service ([Signpost](https://github.com/NCI-GDC/signpost))
    *   Provides an internal identity service used for assigning Digital IDs to all GDC data
    *   Tracks the location of data in GDC object stores
*   Authentication and Authorization Service
    *   Provides a prototype service for authenticating using eRA Commons as the identity provider and dbGaP for authorization

### Bugs Fixed Since Last Release

*   None to report

### Known Issues and Workarounds

*   Identify Service ([Signpost](https://github.com/NCI-GDC/signpost))
    *   The internal identity service does not support authentication for create, read, update, and delete operations as this is an internal service
*   Authentication and Authorization Service
*   Obtaining a token requires browser access. This is a limitation of eRA Commons, not the GDC authorization service. The SAML based protocol eRA Commons uses does not support [ECP](https://wiki.shibboleth.net/confluence/display/CONCEPT/ECP) (or similar)

Release details are maintained in the [GDC API change log](https://github.com/NCI-GDC/gdcapi/blob/develop/CHANGELOG.md).
