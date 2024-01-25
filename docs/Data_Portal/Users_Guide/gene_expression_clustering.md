# Gene Expression Clustering Tool

## Introduction to Gene Expression Clustering

The Gene Expression Clustering tool is a web-based tool for performing sample clustering by selecting a desired set of genes from the NCI Genomic Data Commons (GDC), and visualizing a heatmap of a z-score transformed matrix.

## Accessing the Gene Expression Clustering Heatmap

At the analysis center, click the ‘Gene Expression Clustering’ card to launch the heatmap.

[![Analysis Center Gene Expression Clustering Card](images/analysis_center_GEC_tool.png)](images/analysis_center_GEC_tool.png "Click to see the full image.")

Users can view publicly available genes as well as login with credentials to access controlled data.

## Gene Expression Clustering Tool Features

There are four main panels in the Gene Expression Clustering tool: [controls](#Controls), [heatmap](#Heatmap), [variables](#Variables), and [legend](#Legend). 

[![Gene Expression Clustering Tool Features](images/GEC_tool_features.png)](images/GEC_tool_features.png "Click to see the full image.")

### Controls

The control panel has various functionalities with which users can change or modify the appearance of the matrix. The control panel provides flexibility and a wide range of options to maximize user control.

[![Gene Expression Clustering Tool Controls](images/GEC_tool_controls.png)](images/GEC_tool_controls.png "Click to see the full image.")

__Control Panel:__

* __Clustering:__ Modify the default clustering of the heatmap (Average or Complete), alter the column and row dendrogram measurements, and change the z-score cap
* __Cases:__ Adjust the visible characters of the case labels
* __Genes:__ Modify how cases are represented for each gene (Absolute, Percent, or None), row group and label lengths, rendering style, and the existing gene set
    * __Edit Group:__ Displays a panel of currently selected genes, which can be modified by clicking on a gene to remove it from the gene set, searching for a particular gene to add, loading top variably expressed genes, or loading a pre-defined gene set provided by the MSigDB database
    * __Create Group:__ Create a new gene set by searching for a particular gene, loading top mutated genes, or loading a pre-defined gene set provided by the MSigDB database
* __Variables:__ Search and select variables to add to the matrix below the heatmap
* __Cell Layout:__ Modify the format of the cells by changing colors, cell measurements, and label formatting
* __Legend Layout:__ Alter the legend by changing the font size, measurements, and other formatting preferences
* __Download:__ Download the plot in svg format
* __Zoom:__ Adjust the zoom level by using the up and down arrows on the input box, entering a number, or using the sliding scale to view the case lables
* __Undo:__ Undo changes made to the matrix
* __Redo:__ Redo changes made to the matrix
* __Restore:__ Restore the matrix to its default settings

### Heatmap

The Gene Expression Clustering heatmap displays the active cohort's cases along the top horizontally, genes along the left column, and the z-score transformed gene expression value. 

[![Gene Expression Clustering Tool Heatmap](images/GEC_tool_heatmap.png)](images/GEC_tool_heatmap.png "Click to see the full image.")

Hovering over a cell in the heatmap displays the case submitter_id, gene name, and gene expression value.

[![Gene Expression Clustering Tool Heatmap Cell](images/GEC_tool_cell.png)](images/GEC_tool_cell.png "Click to see the full image.")

Clicking on a cell also gives users the option to launch the [Disco plot](oncomatrix.md/#disco-plot), a circos plot displaying copy number data and consequences for that case.

#### Selecting cases on the cluster

Cases on the cluster can be selected by clicking on the dendrogram. Once part of the dendrogram is selected, users can choose to zoom in to the cases, list all highlighted cases, or create a cohort of the selected cases.

[![Gene Expression Clustering Tool Heatmap Cases Dendrogram](images/GEC_tool_heatmap_cases.png)](images/GEC_tool_heatmap_cases.png "Click to see the full image.")

#### Clicking a case label

Users can click on a case in the dendrogram to showcase the Disco plot or the GDC [Case Summary Page](quick_start.md/#cohort-case-table).

[![Gene Expression Clustering Tool Heatmap Case Selection](images/GEC_tool_heatmap_case_selection.png)](images/GEC_tool_heatmap_case_selection.png "Click to see the full image.")

#### Clicking a gene label

In the column of genes on the left, users can click on a gene to rename it, launch the [ProteinPaint Lollipop plot](proteinpaint_lollipop.md), display the GDC [Gene Summary Page](mutation_frequency.md/#gene-and-mutation-summary-pages), or remove the gene.

[![Gene Expression Clustering Tool Gene Selection](images/GEC_tool_gene.png)](images/GEC_tool_gene.png "Click to see the full image.")

### Variables

Any variables added to the matrix appear below the heatmap. Users can hover over a cell to display the case submitter_id and their value for the given variable.

[![Gene Expression Clustering Tool Variables](images/GEC_tool_variables.png)](images/GEC_tool_variables.png "Click to see the full image.")

#### Clicking a Variable

Click on a variable to rename it, edit it by excluding categories, replace it with a different variable, or remove it entirely. 

[![Gene Expression Clustering Tool Variable Selection](images/GEC_tool_variable_selection.png)](images/GEC_tool_variable_selection.png "Click to see the full image.")

### Legend

In addition to the color coding system for the gene expression values, the legend displays the number of cases from the active cohort in each category for all variables that are selected to appear in the matrix. 

[![Gene Expression Clustering Tool Legend](images/GEC_tool_legend.png)](images/GEC_tool_legend.png "Click to see the full image.")

#### Interacting with legend filters

Users can click on a variable in the legend to hide a specific category, only show a specific category, or show all categories for the selected variable.

[![Gene Expression Clustering Tool Legend Selection](images/GEC_tool_legend_selection.png)](images/GEC_tool_legend_selection.png "Click to see the full image.")