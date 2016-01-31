## 4.1 Overview

The GDC API supports a wide range of query string operators. Users can use a query string operation by placing the query terms and a ```?``` after the endpoint in the URL ```https://gdc-api.nci.nih.gov/projects?```. This allows users to select, filter and order the desired data.

## 4.2 GDC Query Parameters

The following query parameters can be used with all methods and resources in the GDC API. The use of any particular parameter is optional except where noted.

### 4.2.1 Facets
The _facets_ query parameter provides a list of document counts for each included facet.

A list of all valid _field_ names that can be used as facets is available in [Appendix A](/developers/gdc-application-programming-interface-api-users-guide/appendix-available-fields "Appendix A - Available Fields").

**Example:** To get a count of projects in each program, facets=program.name can be passed to the projects endpoint:

    curl  'https://gdc-api.nci.nih.gov/projects?facets=program.name&from=1&size=0&sort=program.name:asc&pretty=true'
    {
      "data": {
        "pagination": {
          "count": 0,
          "sort": "program.name:asc",
          "from": 1,
          "pages": 46,
          "total": 46,
          "page": 1,
          "size": 0
        },
        "hits": [],
        "aggregations": {
          "program.name": {
            "buckets": [
              {
                "key": "TCGA",
                "doc_count": 37
              },
              {
                "key": "TARGET",
                "doc_count": 9
              }
            ]
          }
        }
      },
      "warnings": {}
    }


### 4.2.2 Fields
This query option specifies which fields to include in the response. Using fields can help improve performance. There are about 2,000 fields currently available in the data model. A complete listing of all valid fields for each endpoint type are available in [Appendix A](/developers/gdc-application-programming-interface-api-users-guide/appendix-available-fields "Appendix A - Available Fields").

**Example:** To get back only the file names for each file, fields=file_name can be passed to the files endpoint:

    curl 'https://gdc-api.nci.nih.gov/files?fields=file_name&pretty=true'
    {
      "data": {
        "hits": [
          {
            "file_name": "Collagen_VI-R-V_GBL1112757.tif"
          },
          {
            "file_name": "DUNGS_p_TCGA_b84_115_SNP_N_GenomeWideSNP_6_C09_771624.birdseed.data.txt"
          },
          {
            "file_name": "unc.edu.ba350477-3c49-4258-8409-f39447687497.2145690.junction_quantification.txt"
          },
          {
            "file_name": "SWEDE_p_TCGAb322_23_24_25_26NSP_GenomeWideSNP_6_C06_1364960.hg18.seg.txt"
          },
          {
            "file_name": "unc.edu.a696be4d-8576-43a3-8137-30778dff9b89.2445604.junction_quantification.txt"
          },
          {
            "file_name": "MSK_252152921979_S01_CGH-v4_10_27Aug08__GCN_V3_A1__CBS_out.txt"
          },
          {
            "file_name": "9721366028_R05C01_Red.idat"
          },
          {
            "file_name": "28bdce76404c4a5a9a766c6e93f390ed.bam"
          },
          {
            "file_name": "PKC-delta_pS664-R-V_GBL9016761.tif"
          },
          {
            "file_name": "jhu-usc.edu_LAML.HumanMethylation450.2.lvl-2.TCGA-AB-3011-03A-01D-0742-05.txt"
          }
        ],
        "pagination": {
          "count": 10,
          "sort": "",
          "from": 1,
          "pages": 54777,
          "total": 547761,
          "page": 1,
          "size": 10
        }
      },
      "warnings": {}
    }

### 4.2.3 Filters

Using the query option _filters_  lets users specify criteria for the returned response. See [Filter Available Operations (op)](#avaiable_ops).

The _filters_ syntax is a JSON object that contains the query filters that are translatable to ElasticSearch JSON-based queries by the GDC middleware

Users can get a list of available values for a specific field in the filter by making a call to the appropriate API endpoint using the 'facets' parameters.

**Example**: To get a list of available values for the clinical.gender field, users can query the cases endpoint using the 'facets' parameter:

  curl https://gdc-api.nci.nih.gov/cases?facets=clinical.gender&from=1&size=0&sort=clinical.gender:asc

Filters support complex nested operations as well as simple queries on a single field. There are different types of operations available for many uses. For more examples see [7 - Examples] (/node/8207/).

It is possible to obtain multiple values from multiple fields in one single query, by listing fields ```(facets=field1,field2)```.

**Example:** Filter cases keeping only 'male'. The JSON to be passed to the filter parameter looks like:

    {"op": "=",
          "content": {
              "field": "cases.clinical.gender",
              "value": ["male"]
          }
    }

The above JSON is URL encoded and passed to the filters parameter for the API call:

  curl  'https://gdc-api.nci.nih.gov/cases?filters=%7b%22op%22%3a+%22%3d%22%2c%0d%0a++++++%22content%22%3a+%7b%0d%0a++++++++++%22field%22%3a+%22cases.clinical.gender%22%2c%0d%0a++++++++++%22value%22%3a+%5b%22male%22%5d%0d%0a++++++%7d%0d%0a%7d'

More in depth examples of various filter types supported in GDC are available in the [Appendix A](/developers/gdc-application-programming-interface-api-users-guide/appendix-available-fields "Appendix A - Available Fields")

##### 4.2.3.1 Available Filtering Operations

<a id="available_ops" name="available_ops"></a>

Operators allow users to define query conditions. These can be used to restrict facet values and then to connect these in logical statements.

Operators can relate an operation to one field (Single Field Operators, e.g. A = B) or multiple fields (e.g. [and (a,b,c,d)].

Operators (**op** in the examples in Section 6.2) can take different values depending of the context and type of data.

| Type | Possible Values |
| --- | --- |
| Single field | =, != , <, <=, =, >, >=, in, is, not, range, exclude |
| Multiple fields | and, or |

When using multiple fields, operator content requires nested data containing additional operators.

### 4.2.4 From
The GDC API uses pagination, the **from** query parameter specifies the first record to return out of the set of results. For example, if there are 20 cases returned from the case endpoint, **from** can be set to 10 and results 10-20 will be returned. The **from** parameter can be used in conjunction with the **size** parameter to return a specific subset of results. For more information see 7.1.2 Filters examples.

**Example:** To get 5 file results starting from the 100th result from a set of 500 file results.

<a id="from_example" name="from_example"></a>

    curl 'https://gdc-api.nci.nih.gov/files?fields=file_name&from=101&size=5&pretty=true'
    {
      "data": {
        "hits": [
          {
            "file_name": "unc.edu.956590e9-4962-497b-a59f-81ee0a1c0caf.1379536.junction_quantification.txt"
          },
          {
            "file_name": "MATZO_p_TCGAb40_SNP_1N_GenomeWideSNP_6_G09_667760.ismpolish.data.txt"
          },
          {
            "file_name": "GIRTH_p_TCGA_b108_137_SNP_N_GenomeWideSNP_6_D06_787864.hg18.seg.txt"
          },
          {
            "file_name": "PLENA_p_TCGAb63and64_SNP_N_GenomeWideSNP_6_B12_697382.CEL"
          },
          {
            "file_name": "TCGA-HU-8604-01A-11R-2402-13.isoform.quantification.txt"
          }
        ],
        "pagination": {
          "count": 5,
          "sort": "",
          "from": 100,
          "pages": 109553,
          "total": 547761,
          "page": 21,
          "size": 5
        }
      },
      "warnings": {}
    }

### 4.2.5 Size
The **size** query parameter specifies the number of results to return. When size is not specified the default is 10.

### 4.2.6 Sort
The **sort** query parameter specifies a single field to sort the returned results by sort order, **asc** for ascending order use **des** for descending order: ```sort=field:asc```.

**Example:** Sort returned cases by ascending order by submitter ID:
   curl  'https://gdc-api.nci.nih.gov/cases?sort=submitter_id:asc'


### 4.2.7 Pretty
Returns response with indentations and line breaks in a human-readable format.

**Example:** Sort returned cases by ascending order by submitter ID in human-readable format:
   curl  'https://gdc-api.nci.nih.gov/cases?sort=submitter_id:asc&pretty=true'

### 4.2.8 Format
Returns the response data in a format other than JSON.

The following options are available:

*   TSV
*   XML
