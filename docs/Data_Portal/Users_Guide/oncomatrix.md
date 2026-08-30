# Oncomatrix

## Introduction

The OncoMatrix tool is a web-based tool for visualizing coding mutations such as Simple Somatic Mutations (SSM) and Copy Number Variations (CNV) from the NCI Genomic Data Commons (GDC).

## Accessing the Matrix Chart

At the Analysis Center, click on the “OncoMatrix” card to launch the app.
[![Analysis Tools with Oncomatrix Card](./images/oncomatrix/1-analysis_center.png)](./images/oncomatrix/1-analysis_center.png 'Click to see the full image.')

View publicly available genes as well as login with credentials to access controlled data.

There are three main panels in the OncoMatrix tool: [control panel](#control-panel), [matrix plot](#matrix-plot), and [legend panel](#legend-panel).

[![Oncomatrix chart overview](./images/oncomatrix/2-entire_matrix.png)](./images/oncomatrix/2-entire_matrix.png 'Click to see the full image.')

Each of the features and functionalities are described in detail in the following sections.

### Control Panel

The control panel has various functionalities with which users can change or modify the appearance of the matrix. The control panel provides flexibility and a wide range of options to maximize user control.

[![OncoMatrix Control Panel](./images/oncomatrix/18-control-panel.png)](./images/oncomatrix/18-control-panel.png 'Click to see the full image.')

__Control Panel:__

* __Cases:__ Choose how to sort the cases, specify the maximum number of cases to display, group cases according to selected GDC variables, and adjust the visible characters of the case labels
* __Genes:__ Modify how cases are represented for each gene (Absolute, Percent, or None), row group and label lengths, rendering style, how genes are sorted, the maximum number of genes displayed, and the existing gene set
    * __Edit Group:__ Displays a panel of currently selected genes, which can be modified by clicking on a gene to remove it from the gene set, searching for a particular gene to add, loading top variably expressed genes, or loading a pre-defined gene set provided by the MSigDB database
    * __Create Group:__ Create a new gene set by searching for a particular gene, loading top mutated genes, or loading a pre-defined gene set provided by the MSigDB database
* __Mutation:__ Show or hide specific mutation consequences
* __CNV:__ Show or hide specific CNVs
* __Variables:__ Search and select variables to add to the bottom of the matrix
* __Cell Layout:__ Modify the format of the cells by changing colors, cell dimensions and spacing, and label formatting
* __Legend Layout:__ Alter the legend by changing the font size, dimensions and spacing, and other formatting preferences
* __Download:__ Download the matrix in svg format
* __Zoom:__ Adjust the zoom level by using the up and down arrows on the input box, entering a number, or using the sliding scale to view the case labels.

## Matrix plot

The OncoMatrix plot displays the genes along the left panel with each column representing a case.

#### Matrix cells

Each column in the matrix represents a case. Hovering over a cell will display the corresponding case submitter_id, gene name, copy number information, and mutation class if any are provided. Clicking on a cell also gives users the option to launch the Disco Plot.

[![Disco Plot button](./images/oncomatrix/oncomatrix_disco_plot_button.png)](./images/oncomatrix/oncomatrix_disco_plot_button.png 'Click to see the full image.')

The Disco Plot is a circular plot that shows all the mutations and CNVs for a given case. The Disco Plot also displays the legend for the mutation class and the CNV.

[![Disco Plot](./images/oncomatrix/oncomatrix_disco_plot.png)](./images/oncomatrix/oncomatrix_disco_plot.png 'Click to see the full image.')

#### Automatic Zoom

To perform an automatic zoom, users can click on and hold a case column then drag the mouse from left to right to form a zoom boundary. From the pop-up window, users can choose to zoom in to the cases, list all highlighted cases, or create a cohort of the selected cases.

[![Dragging and zooming on case columns](./images/oncomatrix/oncomatrix_zoom.png)](./images/oncomatrix/oncomatrix_zoom.png 'Click to see the full image.')

The individual case columns are now visible with a demarcated boundary. Above the cases, a slider has been provided for moving from one view to another to accommodate all cases.

[![Zoom using input control and slider](./images/oncomatrix/oncomatrix_zoom2.png)](./images/oncomatrix/oncomatrix_zoom2.png 'Click to see the full image.')

#### Genes

In the panel of genes on the left, users can hover over a gene to view the number of mutated samples, a breakdown of consequence type, and copy number gain and loss counts.

[![Hover over gene](./images/oncomatrix/oncomatrix_gene_hover.png)](./images/oncomatrix/oncomatrix_gene_hover.png 'Click to see the full image.')

Clicking on a gene opens a pop-up window where users can rename it, launch the [ProteinPaint Lollipop plot](proteinpaint_lollipop.md), display the [Gene Summary Page](mutation_frequency.md#gene-and-mutation-summary-pages), and replace or remove the gene. The lollipop plot displays all cases across the GDC affected by SSMs in the selected gene.

[![Click on a gene](./images/oncomatrix/oncomatrix_gene_click.png)](./images/oncomatrix/oncomatrix_gene_click.png 'Click to see the full image.')

#### Variables

Any variables added to the matrix appear at the bottom of the plot. Users can hover over a cell in a variable row to display the case submitter_id and their value for the given variable.

[![Hover over a variable cell](./images/oncomatrix/oncomatrix_variable_cell.png)](./images/oncomatrix/oncomatrix_variable_cell.png "Click to see the full image.")

Clicking on a variable allows users to rename it, edit it by excluding categories, replace it with a different variable, or remove it entirely.

[![Click on a variable](./images/oncomatrix/oncomatrix_variable_click.png)](./images/oncomatrix/oncomatrix_variable_click.png "Click to see the full image.")

#### Drag and drop genes and variables

By default, the genes in the matrix are sorted in descending order according to which genes have the highest number of rendered cases. Users can override this by dragging and dropping gene and variable row labels to sort the rows manually.

### Legend Panel

Below the matrix, the legend displays color coding for mutation classes, CNV, as well as each variable that is selected to appear in the plot.

[![OncoMatrix Legend](./images/oncomatrix/oncomatrix_legend.png)](./images/oncomatrix/oncomatrix_legend.png 'Click to see the full image.')

Clicking on `Consequences` offers options to show only truncating mutations, show only protein-changing mutations, or hide consequences.

[![OncoMatrix Legend - Consequences button](./images/oncomatrix/oncomatrix_legend_consequences.png)](./images/oncomatrix/oncomatrix_legend_consequences.png 'Click to see the full image.')

Clicking on `CNV` allows users to hide CNV.

[![OncoMatrix Legend - CNV button](./images/oncomatrix/oncomatrix_legend_cnv.png)](./images/oncomatrix/oncomatrix_legend_cnv.png 'Click to see the full image.')

Additionally, users can click on a variable's category to hide a specific group, only show a specific group, or show all groups for the selected variable.

[![OncoMatrix Legend - variable category](./images/oncomatrix/oncomatrix_legend_var_cat.png)](./images/oncomatrix/oncomatrix_legend_var_cat.png 'Click to see the full image.')

## Features

The following features are viewable once the matrix application is loaded.
There are three main panels as outlined in the figure below i.e., the `Control panel`, `Matrix chart`, and the `Legend panel`.

[![Oncomatrix chart overview](./images/oncomatrix/2-entire_matrix.png)](./images/oncomatrix/2-entire_matrix.png 'Click to see the full image.')

Each of the features and functionalities are described in detail in the following sections.

## Matrix plot

### Hovering on sample columns

Each column in the matrix represents a sample.
Hover over sample cells/columns to display information about the sample such as case id, gene name, Copy number information and mutation/mutation class (if any provided) as shown.

[![Hovering over samples](./images/oncomatrix/3-sample_hovering.png)](./images/oncomatrix/3-sample_hovering.png 'Click to see the full image.')

### Drag to zoom

A user may click a row label and drag it while keeping the mouse button down, to sort the rows manually. Click and hold on a column of sample and drag the mouse from left to right to form a zoom boundary as shown in the image and leave the mouse.

[![Drag and zoom](./images/oncomatrix/4-drag2zoom.png)](./images/oncomatrix/4-drag2zoom.png 'Click to see the full image.')

This allows for an automatic zoom as shown. The individual sample columns are now visible with a well demarcated boundary. Above the samples, a slider (as shown in gray) has been provided for moving from one view to another to accommodate all cases.

[![Zoom using input control and slider](./images/oncomatrix/5-slide&zoom.png)](./images/oncomatrix/5-slide&zoom.png 'Click to see the full image.')

Additionally, to have a finer control on the zoom the user may follow the steps outlined in the section - Zooming

### Clicking on Sample columns

In the same zoomed in view as shown above, click on any sample column for TP53. This displays a clickable button `Disco plot` as shown.

[![Clicking on sample column label](./images/oncomatrix/6-clicksamplebtn.png)](./images/oncomatrix/6-clicksamplebtn.png 'Click to see the full image.')

Click on the disco plot button to display a circular plot that shows all the mutations for a given sample as shown.

[![Disco plot](./images/oncomatrix/7-discoplot1.png)](./images/oncomatrix/7-discoplot1.png 'Click to see the full image.')

The disco plot can also be accessed by following steps outlined in the section - Disco Plot

### Clicking on gene/variable labels

Click on `TP53` gene label to display the following options.

[![Sort gene labels](./images/oncomatrix/8-gene_sort_icons.png)](./images/oncomatrix/8-gene_sort_icons.png 'Click to see the full image.')

The first row in the options highlighted by a red box as shown in the image above allows the user to sort rows and move rows up and down (please note that rows can also be moved by dragging and dropping as outlined in section Drag and Drop Gene Label/Variable variable). Every time a sorting icon is clicked the chart will update and reload.

Click the first arrow as shown by clicking the gene label TP53. This will sort the samples against the gene at the top left corner which is TP53 in this example.

[![Click to sort against the top left corner gene](./images/oncomatrix/9-diag_arrow.png)](./images/oncomatrix/9-diag_arrow.png 'Click to see the full image.')

Next, click on the left arrow as shown. This allows for sorting samples against the gene.

[![Sort up a gene row](./images/oncomatrix/10-left_arrow.png)](./images/oncomatrix/10-left_arrow.png 'Click to see the full image.')

Now click the down arrow as shown.
The row with TP53 cases will move below ATRX.

[![Move down a gene row](./images/oncomatrix/11-down_arrow.png)](./images/oncomatrix/11-down_arrow.png 'Click to see the full image.')

Click the gene label `TP53` and click the up arrow as shown.

[![Move up a gene row](./images/oncomatrix/12-up_arrow.png)](./images/oncomatrix/12-up_arrow.png 'Click to see the full image.')

The row containing TP53 cases now moves back up in position 1 above ATRX.

Click TP53 again to showcase the edit menu.

Click on `Replace` as shown above to replace TP53 gene variable with `Primary site` as shown below. The chart updates with the first row as `Primary site` thereby replacing TP53 gene variable as shown below. User may choose to sort samples by clicking the `Primary site` label.

[![Replace a gene row](./images/oncomatrix/13-replace_gene.png)](./images/oncomatrix/13-replace_gene.png 'Click to see the full image.')

[![Replaced row](./images/oncomatrix/14-term-replaced.png)](./images/oncomatrix/14-term-replaced.png 'Click to see the full image.')

Click on the label `Primary site` and click the option `Remove` as shown to remove the row completely. 

[![Remove term](./images/oncomatrix/15-remove-primary.png)](./images/oncomatrix/15-remove-primary.png 'Click to see the full image.')

This updates the chart. User may choose to add back TP53 through the gene panel.

Click on `Replace` as shown above to replace TP53 gene variable with `Primary site` as shown below.

Moreover, a user has the option to add rows for gene expression data. Click on `TP53` gene label, click on the option `Gene expression` and search for a gene of interest as shown.

[<img src="./images/oncomatrix/13.1-replace_term_gene_exp.png" width="500"/>](./images/oncomatrix/15-remove-primary.png 'Click to see the full image.')

After a gene expression term is selected, the chart updates with the selected term shown as a separate row at the end of the matrix in a continuous distribution

[<img src="./images/oncomatrix/13.2_gene_exp_term_row.png" width="500"/>](./images/oncomatrix/15-remove-primary.png 'Click to see the full image.')

## Control panel

The control panel as shown has various functionalities with which users can change or modify the appearance of the matrix. The control panel provides flexibility and a wide range of options to maximize user control.

[![Control Panel](./images/oncomatrix/18-control-panel.png)](./images/oncomatrix/18-control-panel.png 'Click to see the full image.')

### Drag and Drop Gene Label/Variable

The genes on the matrix are sorted by default on the number of cases with the gene having the highest number of cases at the top of the matrix. A user may choose to override this by dragging a gene label and dropping it above or below any other gene in order to customize their own gene groupings. 

Select `PTEN` gene label and drag it below the gene labeled `EGFR` as shown. When dragging a gene label, hover over EGFR such that the EGFR gene label would appear blue.

[![Drag a gene row](./images/oncomatrix/16-drag-gene-row.png)](./images/oncomatrix/16-drag-gene-row.png 'Click to see the full image.')

When the EGFR gene label appears blue, then drop the PTEN gene label row. The display updates to show PTEN below EGFR as shown below.

[![Dragged gene row](./images/oncomatrix/17-dragged-gene-row.png)](./images/oncomatrix/17-dragged-gene-row.png 'Click to see the full image.')

### Cases

Within the control panel, the first button displays the number of cases that are shown as columns of the matrix. The default view is as shown.

[![Default view](./images/oncomatrix/19-default-view.png)](./images/oncomatrix/19-default-view.png 'Click to see the full image.')

Click on the `10000 Cases` button to display the following options as shown.

1. Maximum #cases
2. Case Label Character Limit
3. Group Cases by
4. Sort Case Priority

[![Cases menu](./images/oncomatrix/20-cases-controls.png)](./images/oncomatrix/20-cases-controls.png 'Click to see the full image.')

These sections are described below.

#### Maximum #cases

There is a default number of samples that are shown in the matrix chart. Users can choose to increase or decrease the number of samples. This allows the chart to re-render and display the number of columns based on the user's selection. Figure below shows increased cases to 10000. Please note that any high arbitrary number can be selected but the chart will only show the maximum cases that GDC has.

[![Cases menu](./images/oncomatrix/20.1-by_consequence.png)](./images/oncomatrix/20.1-by_consequence.png 'Click to see the full image.')

The chart will reload with new cases added.

#### Case Label Character Limit

This option allows users to increase or decrease the length of the case label. The default number is 32 characters. The chart will reload with new cases added.

#### Group cases by

This option allows users to group cases by different variables from the GDC dictionary. Click on the `+` icon shown in blue to display different variables such as demographics, diagnoses, Exposures etc. Users may also search for a variable from the search bar provided in the menu as shown by `Search Variables`.

[![Group cases by another variable](./images/oncomatrix/group-cases-by.png)](./images/oncomatrix/group-cases-by.png 'Click to see the full image.')

Search for the term `Age at diagnosis` and Click the term. The matrix reloads to show the following view.

[![Select a grouping variable](./images/oncomatrix/22-group-age-at-diag.png)](./images/oncomatrix/22-group-age-at-diag.png 'Click to see the full image.')

As shown above, labels for different age groups show up vertically and all cases get distributed with a clearcut separation according to the age bins.

Click on the `Age at diagnosis` (blue pill, as shown). This opens a short menu with action items. Click on the first item `Edit` as shown.

[![Click the variable button to edit](./images/oncomatrix/23-edit-groups.png)](./images/oncomatrix/23-edit-groups.png 'Click to see the full image.')

Drag the red lines on the density distribution to select binning or input numbers for custom binning and select `Apply`.

[![Click the variable button to edit bins](./images/oncomatrix/24-binning.png)](./images/oncomatrix/24-binning.png 'Click to see the full image.')

The matrix reloads with new bin groupings The labels for the groups are user controlled and hence can be modified according to user requirements. 

Click on the blue pill for `Age at diagnosis` again and click `Replace`. Select `Primary site` as shown.

[![Replace the grouping variable](./images/oncomatrix/25-replace-group.png)](./images/oncomatrix/25-replace-group.png 'Click to see the full image.')

The matrix reloads with the new variable distribution.

The last option on the menu is `Remove`. Click on the `917 Cases` button, followed by `Age at diagnosis` shown in blue to reveal the menu option. Click `Remove` to completely get rid of any groups.

[![Remove the grouping variable](./images/oncomatrix/26-remove-group.png)](./images/oncomatrix/26-remove-group.png 'Click to see the full image.')

This will remove all and any groupings and show the default view again.

### Adding Survival term as a grouping variable

Click on the `Age at diagnosis` (blue pill, as shown). This opens a short menu with action items. Click on the second item `Replace` as shown.

[![Click the variable button to edit](./images/oncomatrix/23-edit-groups.png)](./images/oncomatrix/23-edit-groups.png 'Click to see the full image.')

Select the term 'Overall Survival'. 

[![Selecting survival term](./images/oncomatrix/OS_term.png)](./images/oncomatrix/OS_term.png 'Click to see the full image.')

The matrix reloads with the new variable distribution. Cases are shown for categories 'Alive' vs 'Dead'.

[![Selecting survival term](./images/oncomatrix/grouped_by_survival.png)](./images/oncomatrix/grouped_by_survival.png 'Click to see the full image.')

To view a continuous distribution of the survival outcome data, replace a gene row with a survival term as shown. Click on 'Overall survival'.

[![Selecting survival term](./images/oncomatrix/replace_tp53_survival.png)](./images/oncomatrix/replace_tp53_survival.png 'Click to see the full image.')

The matrix reloads with replacement of gene row 'Tp53' with 'Overall survival'.

[![Replacing gene row](./images/oncomatrix/replace_tp53_survival2.png)](./images/oncomatrix/replace_tp53_survival2.png 'Click to see the full image.')

Once the survival term is loaded, click on the term label, and click edit to display the following options.

[![Edit survival term row](./images/oncomatrix/edit_survival_term.png)](./images/oncomatrix/edit_survival_term.png 'Click to see the full image.')

Now the user has the option to display the survival term as a continuous distribution, or to show it as a discrete distribution. The default is a discrete distribution with no z-scoring.

Click 'Apply' to load the continuous distribution for 'Time to event'

This reloads the matrix and now displays the survival outcome as a continuous distribution. 

[![Survival term as continuous distribution](./images/oncomatrix/OS_cont_distribution.png)](./images/oncomatrix/OS_cont_distribution.png 'Click to see the full image.')

#### Sort Case Priority

The default sort setting sorts the cases 'by presence' under 'Basic' sort settings. 

Click the second option `by consequence` to change the sorting. The matrix reloads with the new sorting as shown below.

[![Sort oncomatrix by consequence](./images/oncomatrix/20.1-by_consequence.png)](./images/oncomatrix/20.1-by_consequence.png 'Click to see the full image.')

To perform an advanced sorting, click 'Advanced' on the 'Sort Case Priority' menu as shown below. 

[![Advanced sorting options](./images/oncomatrix/20.2-sort_advanced.png)](./images/oncomatrix/20.2-sort_advanced.png 'Click to see the full image.')

Now user has the option to sort the cases by each selected row, gene mutation, dictionary variable or alphabetically by name. Details of each sort option are provided. 

### Sort Case Groups

Add the variable `Age at diagnosis` again using the `Group Cases by` button as shown in the previous section. By default, groups are loaded ordered by their name. Change the selection to `Case count` as shown below. 

[![Sorting Case groups](./images/oncomatrix/27-sort-case-groups.png)](./images/oncomatrix/27-sort-case-groups.png 'Click to see the full image.')

The third selection option `Hits` orders the groupings based on the number of gene variants for a particular case or case group for the genes in display. Click `Hits` under `Sort Case Groups` to change the order of groupings.

Next, hover over the group label `<=30 years (81)` as shown below.

[![Hover over lables](./images/oncomatrix/28-label-hover.png)](./images/oncomatrix/28-label-hover.png 'Click to see the full image.')

This shows the number of cases in parenthesis of the group label and the breakdown for the  number of variants and CNV for all the samples within that group for the genes in display.

## Genes

The gene panel as shown below has several options as listed below for modifying the genes visible on the plot as well as their appearance/style.

- Display Case Counts for Gene
- Genomic Alterations Rendering
- Sort Genes
- Maximum # Genes
- Gene Set

[![Genes button](./images/oncomatrix/29-genes-btn.png)](./images/oncomatrix/29-genes-btn.png 'Click to see the full image.')

### Display Case Counts for Gene

This option allows change in the number of cases that is represented in parentheses next to the gene variable label as shown below. By default, the number of cases for each gene is an `Absolute`.

Click on the button `50 Genes` to display the menu and select `Percent`

This shows the case counts as a percentage of the absolute values as shown. 

[![Counts as percent of cases](./images/oncomatrix/30-genes-percent-display.png)](./images/oncomatrix/30-genes-percent-display.png 'Click to see the full image.')

User has the option to hide the display of case counts. Click `Genes` button again and select `None` for `Display Case Counts for Gene` as shown below

[![Genes None](./images/oncomatrix/31-genes_none1.png)](./images/oncomatrix/31-genes_none1.png 'Click to see the full image.')

This hides all the case counts as shown.

[![Hide case counts](./images/oncomatrix/32-genes-none-display.png)](./images/oncomatrix/32-genes-none-display.png 'Click to see the full image.')

### Genomic Alterations Rendering

The style of rendering for the sample cells/columns is an Oncoprint style by default. Click on `Stacked` option via `50 Genes` button as shown below.

[![Rendering Style](./images/oncomatrix/33-rendering_stacked1.png)](./images/oncomatrix/33-rendering_stacked1.png 'Click to see the full image.')

The mutations and CNV are now stacked on top of each other as shown below.

[![Stacked rendering](./images/oncomatrix/34-stacked-view.png)](./images/oncomatrix/34-stacked-view.png 'Click to see the full image.')

To view the rendering in an oncoprint style, click `Oncoprint` button on the control panel. This updates the rendering as shown.

### Sort Genes

The default sorting option for genes is `By Case Count`. This means the genes are sorted by the number of cases from increasing to decreasing order. Click `50 Genes` button on the control panel, and select `By Input Data Order` under the `Sort Genes` as shown below.

[![Sort genes](./images/oncomatrix/35-sort-genes.png)](./images/oncomatrix/35-sort-genes.png 'Click to see the full image.')

The genes will now sort according to the order that is stored in the dataset and queried. However, please note that the sorting order can be overridden by the users choice as described in the section - Drag and Drop Gene Label/Variable.

### Maximum # Genes

The number of genes to display on the matrix plot can be modified by the input option as shown below. Click `50 Genes` button and change input number for `Maximum # Genes` to 70.

[![Maximum number of genes to display](./images/oncomatrix/36-max_genes.png)](./images/oncomatrix/36-max_genes.png 'Click to see the full image.')

The chart updates and loads the extra 20 genes. User can modify the set of genes by using the `Gene set` option next.

### Editing gene set

Gene groups can be edited using the `Gene set` option as shown below. Click `50 Genes` button to display this option and then click the `Edit` button in the `Gene set` as shown.

[![Edit geneset](./images/oncomatrix/37-geneset_edit1.png)](./images/oncomatrix/37-geneset_edit1.png 'Click to see the full image.')

User may choose to remove single genes one at a time by clicking over the genes.To do so, hover over TP53 as shown in the image below. A red cross mark appears with a description box. Click TP53 to delete the gene as shown below.

[![Deleting a gene](./images/oncomatrix/39-geneset_edit3.png)](./images/oncomatrix/39-geneset_edit3.png 'Click to see the full image.')

User may choose to delete all genes from view by clicking the `Clear` button as shown below. However, a gene/variable selection is mandatory for the chart to load.

[![Clear genes](./images/oncomatrix/40-clear_genes.png)](./images/oncomatrix/40-clear_genes.png 'Click to see the full image.')

For loading user specific gene sets, click the 'Load gene set' button within the geneset edit panel. This shows an input option for the user to add a comma separated gene set of their own.

[![Load gene set](./images/oncomatrix/load_gene_set.png)](./images/oncomatrix/load_gene_set.png 'Click to see the full image.')

### MSigDB genes

The MSigDB database (Human Molecular Signatures Database) has 33591 gene sets divided into 9 major collections and several subcollections. Users can choose to view the gene sets on the matrix plot.

Click on the `50 Genes` button. Then click on the `Gene set - Edit`. Here user can see a button with a dropdown for loading MSigDB genes.

[![msigdb dropdown](./images/oncomatrix/41-msigdb1.png)](./images/oncomatrix/41-msigdb1.png 'Click to see the full image.')

Click on this dropdown to display a tree for the different gene sets.

[![msigdb tree](./images/oncomatrix/42-msigdb_tree.png)](./images/oncomatrix/42-msigdb_tree.png 'Click to see the full image.')

Select `C2: curated gene sets` and select `NABA_COLLAGENS` as shown below.

[![NABA collagens gene set](./images/oncomatrix/43-naba_collagens.png)](./images/oncomatrix/43-naba_collagens.png 'Click to see the full image.')

This loads the following genes as shown below.

[![NABA collagens genes loaded](./images/oncomatrix/44-msigdb2.png)](./images/oncomatrix/44-msigdb2.png 'Click to see the full image.')

Click `Submit` and the matrix will update to reflect the selected MSigDB gene set as shown.

[![NABA collagens genes matrix](./images/oncomatrix/45-msigdb_collagen_matrix.png)](./images/oncomatrix/45-msigdb_collagen_matrix.png 'Click to see the full image.')

## Mutation

To view specific mutations click on the 'Mutations' tab in the control panel to display the following options:

[![Show mutations](./images/oncomatrix/mutations-tab.png)](./images/oncomatrix/mutations-tab.png 'Click to see the full image.')

1. Show all mutations - selecting this will display all the mutations for the cohort.
2. Show only truncating mutations - this option only shows truncated mutations and hides other mutations 
3. Show only protein-changing mutations - this option only shows mutations that alter the protein and hides other mutations. This is the default option.
4. Do not show mutations - this option hides all the mutations and only displays CNV cases
5. Show selected mutations - here a user can specifically select the consequence of mutations they want to visualize

[![Show selected mutations](./images/oncomatrix/select-mutations.png)](./images/oncomatrix/select-mutations.png 'Click to see the full image.')

# Edit mutation data on gene row

To create distinct grouping on variants within a gene row, click 'edit' on the gene row to view the following options.

[![Edit mutation data on gene row](./images/oncomatrix/edit_mutation_data_gene_row.png)](./images/oncomatrix/edit_mutation_data_gene_row.png 'Click to see the full image.')

Here user has the option to cancel any groupings on variant data or depending on the variant type (SNV/Indel or CNV), user can further create custom groups. The pre-defined groups 'Wildtype' and 'Mutated' are shown. 

## CNV

By default the CNVs are hidden. User can select the `Show all CNV` option to display the CNV data by clicking the 'CNV' tab from the control panel.

[![CNV tab](./images/oncomatrix/CNV-tab.png)](./images/oncomatrix/CNV-tab.png 'Click to see the full image.')

This re-loads the matrix with CNVs as shown.

[![CNV data](./images/oncomatrix/CNV-data.png)](./images/oncomatrix/CNV-data.png 'Click to see the full image.')

To select specific CNVs, click on the 'Show selected CNV' tab in the control panel and select the CNVs you want to visualize.

[![Select CNV](./images/oncomatrix/select-CNV.png)](./images/oncomatrix/select-CNV.png 'Click to see the full image.')

## Variables

The third button from the left called `Variables` allows user to add in additional variables in the form of rows on the matrix. Click `Variables` to display a tree of variables and select `Disease type`, `Index date` and `Primary site`. Click the button `Submit 3 variables` as shown.

[![Variables menu](./images/oncomatrix/46-variables.png)](./images/oncomatrix/46-variables.png 'Click to see the full image.')

This updates the chart to display the selected variables on the very top of the matrix as shown below. User may choose to configure these rows by following steps outlined in section Clicking on gene/variable labels.

[![Loaded variables](./images/oncomatrix/47-variable-loaded.png)](./images/oncomatrix/47-variable-loaded.png 'Click to see the full image.')

## Cell Layout

The cell layout menu enables customization of the appearance such as cell dimensions, spacing, font sizes, and borders. You may mouseover an input to see the description for that input, or try checking or editing inputs to test the effects of the control input and undo/redo as needed.

[![Cell layout](./images/oncomatrix/48-cell_layout.png)](./images/oncomatrix/48-cell_layout.png 'Click to see the full image.')

## Legend Layout

The legend layout menu enables customization of the appearance of the legend, such as dimensions, spacing, and font sizes. These customizations can help avoid or minimize the need for post-download edits when generating figures.

[![Legend layout](./images/oncomatrix/49-legend_layout.png)](./images/oncomatrix/49-legend_layout.png 'Click to see the full image.')

## Zooming

The matrix plot offers an interactive zoom panel as shown below with which a user can zoom in to view individual samples. There are two ways to use this panel. One by changing the input number and second by sliding the zoom bar to a desired zoom level as shown.

[![Zoom slider](./images/oncomatrix/50-zoom_slider.png)](./images/oncomatrix/50-zoom_slider.png 'Click to see the full image.')

Change zoom level to 10+ as shown.

[![Zoom slider plus](./images/oncomatrix/51-zoom_slider_plus.png)](./images/oncomatrix/51-zoom_slider_plus.png 'Click to see the full image.')

Scroll down to view individual samples at the bottom of the plot as shown below.

[![Zoomed in matrix view](./images/oncomatrix/52-zoomed-in.png)](./images/oncomatrix/52-zoomed-in.png 'Click to see the full image.')

The zoom action can also be implemented by following steps as outlined in section - Drag to zoom.

## Disco Plot

Click on any sample to reveal a second type of plot called as the `Disco Plot` as shown.

[![Sample label clicking to show disco plot button](./images/oncomatrix/53-disco_plot1.png)](./images/oncomatrix/53-disco_plot1.png 'Click to see the full image.')

Click on `Disco plot` as shown above in gray. This loads a new chart above the matrix plot as shown below.

[![Disco plot](./images/oncomatrix/54-disco_plot2.png)](./images/oncomatrix/54-disco_plot2.png 'Click to see the full image.')

This plot shows all the mutations and CNV associated with that sample id as shown above. The plot also displays the legend for the mutation class and the CNV.

To reset the zoom level to default, click on the `Reset` button as shown. This will reset the zoom level to a default of 1.0

[![Reset Zoom](./images/oncomatrix/55-reset_zoom.png)](./images/oncomatrix/55-reset_zoom.png 'Click to see the full image.')

## Download

The control panel shows an option to download the plot as an svg after user has specified their customizations. Select the `Download` button as shown below to save the svg.

[![Download button](./images/oncomatrix/59-download_svg1.png)](./images/oncomatrix/59-download_svg1.png 'Click to see the full image.')

If svg format is selected then the download will get saved to the default download folder as shown at the bottom of the browser window.

[![Download](./images/oncomatrix/60-download_svg2.png)](./images/oncomatrix/60-download_svg2.png 'Click to see the full image.')

## Legend

The legend for the matrix is below the plot and shows color coding for different mutation classes as well as color codes for CNV as shown here. This legend is interactive and user may choose to hide or show features such as mutation classes or copy number changes.

Click on the legend icons to hide anything.

[![Legend](./images/oncomatrix/61-legend.png)](./images/oncomatrix/61-legend.png 'Click to see the full image.')
