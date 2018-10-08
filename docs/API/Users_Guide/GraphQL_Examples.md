Introduction to GDC GraphQL #
[GraphQL](https://graphql.org/) is a query language for APIs. The [GDC REST API](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/) has structured and specifically defined query parameters as well as endpoints that have set requests and responses. The GDC GraphQL provides advanced GDC developers greater flexibility to specify the data they would like to be returned. This allows queries to be cleaner and easier to understand, especially when combining multiple queries into one request.

## Using GDC GraphQL vs GDC REST API?

If the use case does not require all of the data to be returned, GDC GraphQL may speed up requests as GraphQL queries return only the specified data, therefore this may require less work on the GDC server-side to fulfill those requests. Conversely, if you need all of the data in each request, GDC REST API may be a better fit. Either way,the data itself returned from the GDC REST API or the GraphQL query will be identical.

##  GDC GraphQL Overview ##
[GraphQL](http://graphql.org/) is not a storage model or a database query language. The graph refers to graph structures defined in the schema, where nodes define objects and edges define relationships between objects. The API traverses and returns application data based on the schema definitions, independent of how the data is stored. In other words:

## GDC GraphQL Endpoints

The GDC GraphQL has only two endpoints:
__GDC Search & Retrieval Endpoint:__ https://api.gdc.cancer.com/v0/graphql
__GDC Submission Endpoint:__ https://api.gdc.cancer.com/v0/submission/graphql

## GDC GraphQL Schema ##
All GDC GraphQL queries are validated and executed against the [[GDC GraphQL schema]( https://github.com/NCI-GDC/portal-ui/blob/92f0dfa17838746093c3c011141d08391016da91/data/schema.graphql). Since GraphQL is introspective the GDC GraphQL schema can be queried for details about itself.

Query `__schema` to list all types defined in the schema and retrieve details about each:

```Run in Explorer
    {__schema
      { types
	    {
	      name
	      kind
	      fields{
		      name
		  }
	    }
	  }
    }
```
The `__type` keyword can also be queried to retrieve details about any type such as "Explore":
```Run in Explorer

    query {
      __type(name: "Explore") {
		    name
		    kind
		    description
		    fields {
    			name
	    		}
			}
    	}
```
```Run in Explorer

    query {
      __type(name: "Case") {
    	name
    	kind
    	description
    	fields {
    	  name
    		}
    	  }
    	}
```

## Basic GraphQL queries in GDC
The two types of allowed operations in GDC GraphQL API are queries and mutations. Comparing GraphQL to REST, queries operate like `GET` requests, while mutations operate like `POST`/`PATCH`/`DELETE`.

**Note:** This guide does not cover GDC GraphQL mutation operations.

GraphQL queries return only the data that is specified. To form a query, you must specify fields within fields (also called nested *subfields*) until only scalars are returned.  Scalars are primitive values such as: `Int`, `Float`, `String`, `Boolean`, or `ID`.

## Anatomy of a typical GDC GraphQL Query ##:

(see graphql-query.png)


- **Operation type**: Describes what type of operation you’re trying to do, such as query, mutation, or subscription
- **Operation name:** Similar to a function name, gives queries meaningful names
- **Field:** Denotes the specific fields on objects that will be included with the response data
- **Arguments:** A set of key-value pairs associated with a specific field. The parameters can be literal values or variables. **NOTE**: Arguments can appear on any field, even fields nested deep in an operation.
- **Variable definitions:** As GraphQL is strong typed, it validates the variable being passed dynamically. **NOTE**:Variables are passed separately from the query document as JSON such as:
```json
    { "filters_1": {"op":"and","content":[{"op":"in","content":{"field":"projects.program.name","value":["TARGET"]}}]} }
```

## GDC Graphql Examples ##
### Nodes And Edges Example ###
A very powerful feature of GDC Graphql API is that the graph structures defined in the [GDC GraphQL schema]( https://github.com/NCI-GDC/portal-ui/blob/92f0dfa17838746093c3c011141d08391016da91/data/schema.graphql ) can be queried and traversed, where nodes define objects and edges define relationships between objects.

```Run in Explorer

    query PROJECTS_EDGES($filters_1: FiltersArgument) {
      projects {
    	hits(filters: $filters_1) {
      		total
      		edges {
    			node {
    			  primary_site
    			  disease_type
    			  project_id
    			  dbgap_accession_number
    			}
      		}
    	}
      }

    variable:
    { "filters_1": {"op": "and", "content": [{"op": "in", "content": {"field": "projects.primary_site", "value": ["Kidney"]}}]}}
```

### Query Case File Counts ###
Run in Explorer```

    query CaseFileCounts($filters: FiltersArgument) {
            viewer {
              repository {
                cases {
                  hits(first: 1, filters: $filters) {
                    edges {
                      node {
                        case_id
                        files {
                          hits(first: 0) {
                            total
                          }
                        }
                        summary {
                          experimental_strategies {
                            experimental_strategy
                            file_count
                          }
                          data_categories {
                            data_category
                            file_count
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }

    variable:
    {"filters":{"op":"and","content":[{"op":"in","content":{"field":"cases.case_id","value":["dcd5860c-7e3a-44f3-a732-fe92fe3fe300"]}}]}}
```
### Query Simple Static Mutations Based on Gene IDs ###
```Run in Explorer

    query PROJECTS_EDGES($filters_1: FiltersArgument) {
      projects {
    	hits(filters: $filters_1) {
    		total
    		edges {
    		node {
      			primary_site
      			disease_type
      			project_id
      			dbgap_accession_number
    			}
    		}
    	}
      }
    }
	variable:
    {"filters_2": {"op":"and","content":[{"op":"in","content":{"field":"consequence.transcript.gene.gene_id","value":["ENSG00000155657"]}}]}}```

