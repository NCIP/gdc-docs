# Projects

At a high level, data in the Genomic Data Commons is organized by project. Typically, a project is a specific effort to look at particular type(s) of cancer undertaken as part of a larger cancer research program. The GDC Data Portal allows users to access aggregate project-level information via the Projects View and project detail pages.

## Projects View

The Projects View provides an overview of all harmonized data available in the Genomic Data Commons, organized by project. It also provides filtering and navigation features that allow users to identify and browse projects of interest. Users can access Projects View from the GDC Data Portal front page, from the Data Portal toolbar, or directly at [https://portal.gdc.cancer.gov/projects](https://portal.gdc.cancer.gov/projects).

On the right, the Table tab provides a list of projects and select details about each project, such as the number of cases and data files. The Graph tab provides a visual representation of this information.

On the left, a panel of data facets allows users to apply filters to find projects of interest. When facet filters are applied, the tabs on the right are updated to display only the matching projects. When no filters are applied, all projects are displayed.

[![Projects View, Main Window (Table View)](images/gdc-data-portal-project-page.png)](images/gdc-data-portal-project-page.png "Click to see the full image.")

### Projects Table

The Table tab lists projects by Project ID and provides additional information about each project. If no facet filters have been applied, the table will display all available projects; otherwise it will display only those projects that match the selected criteria.

The table provides links to project detail pages in the Project ID column. Columns with file and case counts include links to open the corresponding files or cases in [Data View](Cases_and_Files.md#data-view).

### Projects Graph

The Graph tab contains an interactive view of information in the Table tab. The numerical values in Case Count, File Count, and File Size columns are represented by bars of varying length according to size. These columns are sorted independently in descending order. Mousing over an element of the graph connects it to associated elements in other columns, including project name and primary site.

[![Graph Mouseover](images/gdc-table-graph-mouse-over.png)](images/gdc-table-graph-mouse-over.png "Click to see the full image.")

Most elements in the graph are clickable, allowing the user to open the associated cases or files in [Data View](Cases_and_Files.md#data-view).

Like the projects table, the graph will reflect any applied facet filters.

### Facets Panel

Facets represent properties of the data that can be used for filtering. The facets panel on the left allows users to filter the projects presented in the Table and Graph tabs. Users can filter by the following facets:

*   __Project__: Individual project ID
*   __Primary Site__: Anatomical site of the cancer under investigation or review.
*   __Cancer Program__: Research program that the project is part of.
*   __Disease Type__: Type of cancer studied.
*   __Data Category__: Type of data available in the project.
*   __Experimental Strategy__: Experimental strategies used for molecular characterization of the cancer.

Filters can be applied by selecting values of interest in the available facets, for example "WXS" and "miRNA-Seq" in the "Experimental Strategy" facet and "Kidney" in the "Primary Site" facet. When facet filters are applied, the Table and Graph tabs are updated to display matching projects, and the banner above the tabs  summarizes the applied filters. The banner allows the user to click on filter elements to remove the associated filters, and includes a link to view the matching cases and files.

[![Panel with Applied Filters](images/panel-with-applied-filters.png)](images/panel-with-applied-filters.png "Click to see the full image.")


For information on how to use facet filters, see [Getting Started](Getting_Started.md#facet-filters).

## Project Detail Page

Each project has its own detail page that provides an overview of all cases, files and annotations available for the project. Clicking on summary numbers on the page will display the corresponding data.

[![Project Entity Page](images/gdc-project-entity-page.png)](images/gdc-project-entity-page.png "Click to see the full image.")

Three download buttons in the top right corner of the screen allow the user to download the entire project dataset along with the associated project metadata:

* __Download Manifest__: Download a manifest of all data files available in the project. The manifest can be used with the GDC Data Transfer Tool to download the files.
* __Download Clinical__: Download clinical metadata about all cases in the project.
* __Download Biospecimen__: Download metadata about all biospecimens available in the project.

### Summary Tables and Pie Charts

The Case and File Counts tables and pie charts break down the available data by experimental strategy and data category.

[![Convert Pie Chart to Table](images/gdc-pie-chart-table.png)](images/gdc-pie-chart-table.png "Click to see the full image.")

Users can switch between chart and table view using the icon in the top right corner of the chart or table.

[![Pie Charts Widget](images/gdc-pie-chart-view.png)](images/gdc-pie-chart-view.png "Click to see the full image.")

Mousing over a slice of the pie chart display detailed about the slice.

[![Mouseover a Slice](images/gdc-pie-chart-mouse-over.png)](images/gdc-pie-chart-mouse-over.png "Click to see the full image.")
