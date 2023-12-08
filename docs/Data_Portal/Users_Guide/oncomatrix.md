# OncoMatrix

## Introduction to OncoMatrix

The GDC OncoMatrix is a handy tool to visualize coding mutations (Simple Somatic Mutations, or SSM) and copy number variations (here onward referred to as CNV).
Each row is a gene and each cell, or a column represents a case. At any point in the tutorial, hover over a symbol or icon for two-three seconds to display more information about that icon.

## Accessing the Matrix Chart

At the Analysis Center, click on the “OncoMatrix” card to launch the app.
[![Analysis Tools with Oncomatrix Card](./images/oncomatrix/1-analysis_center.png)](./images/oncomatrix/1-analysis_center.png 'Click to see the full image.')

View publicly available genes as well as login with credentials to access controlled data.

## OncoMatrix Features

The following features are viewable once the matrix application is loaded.
There are three main panels as outlined in the figure below i.e., the `Control panel`, `Matrix chart`, and the `Legend panel`.

[![Oncomatrix chart overview](./images/oncomatrix/2-entire_matrix.png)](./images/oncomatrix/2-entire_matrix.png 'Click to see the full image.')

Each of the features and functionalities are described in detail in the following sections.

## Matrix plot

### Hovering on sample columns

Each column in the matrix represents a sample.
Hover over sample cells/columns to display information about the sample such as case id, gene name, Copy number information and mutation/mutation class (if any provided) as shown below.

[![Hovering over a sample column](./images/oncomatrix/3-sample_hovering.png)](./images/oncomatrix/3-sample_hovering.png 'Click to see the full image.')

### Drag to zoom

A user may click a row label and drag it while keeping the mouse button down, to sort the rows manually. Click and hold on a column of sample and drag the mouse from left to right to form a zoom boundary as shown in the image below and leave the mouse.

[![Dragging and zooming on sample columns](./images/oncomatrix/4-drag2zoom.png)](./images/oncomatrix/4-drag2zoom.png 'Click to see the full image.')

This allows for an automatic zoom as shown below. The individual sample columns are now visible with a well demarcated boundary. Above the samples, a slider (as shown in gray) has been provided for moving from one view to another to accommodate the 2000 cases.

[![Zoom using input control and slider](./images/oncomatrix/5-drag&zoomed.png)](./images/oncomatrix/5-drag&zoomed.png 'Click to see the full image.')

Additionally, to have a finer control on the zoom the user may follow the steps outlined in the section - Zooming

### Clicking on Sample columns

In the same zoomed in view as shown above, click on any sample column for TP53. This displays a clickable button `Disco plot` as shown below.

[![Disco plot button on clicking sample column](./images/oncomatrix/6-col_discoBtn.png)](./images/oncomatrix/6-col_discoBtn.png 'Click to see the full image.')

Click on the disco plot button to display a circular plot that shows all the mutations for a given sample as shown below.

[![Disco plot](./images/oncomatrix/7-discoplot1.png)](./images/oncomatrix/7-discoplot1.png 'Click to see the full image.')

The disco plot can also be accessed by following steps outlined in the section - Disco Plot

### Clicking on gene/variable labels

Click on `TP53` gene label to display the following options.

[![Gene sort icons](./images/oncomatrix/8-gene_sort_icons.png)](./images/oncomatrix/8-gene_sort_icons.png 'Click to see the full image.')

The first row in the options highlighted by a red box as shown in the image above allows the user to sort rows and move rows up and down (please note that rows can also be moved by dragging and dropping as outlined in section Drag and Drop Gene Label/Variable variable). Every time a sorting icon is clicked the chart will update and reload.

Click the first arrow as shown below by clicking the gene label `TP53`. This will sort the samples against the gene at the top left corner which is TP53 in this example.

[![Sorting by diagonal arrow](./images/oncomatrix/9-diag_arrow.png)](./images/oncomatrix/9-diag_arrow.png 'Click to see the full image.')

Next, click on the left arrow as shown below. This allows for sorting samples against the gene.

[![Sorting by left arrow](./images/oncomatrix/10-left_arrow.png)](./images/oncomatrix/10-left_arrow.png 'Click to see the full image.')

Now click the down arrow as shown below. This will move the row for that gene downward as shown below.

[![Sorting by down arrow](./images/oncomatrix/11-down_arrow.png)](./images/oncomatrix/11-down_arrow.png 'Click to see the full image.')

In the figure below, the row with TP53 cases has moved below MUC16 as shown below.

[![MUC16 above Tp53](./images/oncomatrix/12-MUC16_up.png)](./images/oncomatrix/12-MUC16_up.png 'Click to see the full image.')

Click the gene label `TP53` and click the up arrow as shown.

[![Sorting by up arrow](./images/oncomatrix/13-up_arrow.png)](./images/oncomatrix/13-up_arrow.png 'Click to see the full image.')

The row containing TP53 cases now moves back up in position 1 above MUC16 as shown.

[![Move Tp53 on top](./images/oncomatrix/14-TP53_up.png)](./images/oncomatrix/14-TP53_up.png 'Click to see the full image.')

Click `TP53` again to showcase the edit menu.

Click `Edit` option to display the interactive legend for mutation classes within TP53. Click on `MISSENSE` as shown below to hide those mutations and click apply. NOTE: There are plans to improve this edit menu.

This will hide the corresponding mutation class from view.

[![Gene Legend](./images/oncomatrix/15-gene_legend.png)](./images/oncomatrix/15-gene_legend.png 'Click to see the full image.')

Click on `Replace` as shown above to replace TP53 gene variable with `Primary site` as shown below.

[![Replace gene variable](./images/oncomatrix/16-replace_gene.png)](./images/oncomatrix/16-replace_gene.png 'Click to see the full image.')

The chart updates with the first row as `Primary site` thereby replacing TP53 gene variable as shown below. User may choose to sort samples by clicking the `Primary site` label.

[![Primary site](./images/oncomatrix/17-primary_site.png)](./images/oncomatrix/17-primary_site.png 'Click to see the full image.')

Click on the label `Primary site` and click the option `Remove` as shown below to remove the row completely.

[![Remove variable](./images/oncomatrix/18-remove_var.png)](./images/oncomatrix/18-remove_var.png 'Click to see the full image.')

This updates the chart and shows all 2000 samples as shown below. Please note the number of genes displayed is now `49`. User may choose to add back TP53 by following the section - Genes.

### Drag and Drop Gene Label/Variable

The genes on the matrix are sorted by default on the number of cases with the gene having the highest number of cases at the top of the matrix. A user may choose to override this by dragging a gene label and dropping it above or below any other gene to customize their own gene groupings.

Select `TP53` gene label and drag it below the gene labeled `KRAS` as shown. When dragging a gene label, hover over KRAS such that the `KRAS` gene label would appear blue.

[![Drag gene label](./images/oncomatrix/19-drag_gene_label.png)](./images/oncomatrix/19-drag_gene_label.png 'Click to see the full image.')

When the `KRAS` gene label appears blue, then drop the TP53 gene label row. The display updates to show TP53 below KRAS as shown below.

[![Dragged gene](./images/oncomatrix/20-gene_dragged.png)](./images/oncomatrix/20-gene_dragged.png 'Click to see the full image.')

## Control panel

The control panel as shown below has various functionalities with which users can change or modify the appearance of the matrix. The control panel provides flexibility and a wide range of options to maximize user control.

[![Control Panel](./images/oncomatrix/21-control_panel.png)](./images/oncomatrix/21-control_panel.png 'Click to see the full image.')

### Cases

Within the control panel, the first button displays the number of cases that are shown as columns of the matrix. The default view is two thousand samples as shown below.

[![Cases](./images/oncomatrix/22-cases.png)](./images/oncomatrix/22-cases.png 'Click to see the full image.')

Click on the `2000 Cases` button to display the following options as shown in the figure below.

- Sort Cases
- Maximum #cases
- Group cases by
- Sort Case Groups
- Case Group Label Max Length
- Case Label Max Length


[![Cases options](./images/oncomatrix/23-cases_options.png)](./images/oncomatrix/23-cases_options.png 'Click to see the full image.')

These sections are described below.

#### Sort Cases

The default sort setting sorts the cases by row with first displaying samples with both CNV and SSM followed by SSM only and lastly CNV only. 
Click the second option `CNV+SSM > SSM only` to change the sorting as shown below.

[![Sort cases by SSM+CNV > SSM only](./images/oncomatrix/24-sort_cases.png)](./images/oncomatrix/24-sort_cases.png 'Click to see the full image.')

Sorting ‘By Case name, ID or label’ sorts on the basis of the order of the name of cases or ID as shown below.

[![Sort cases by name, ID or label](./images/oncomatrix/73-sort_cases_name.png)](./images/oncomatrix/73-sort_cases_name.png 'Click to see the full image.')

#### Maximum #cases

The default number of samples that are shown in the matrix chart is 2000. Users can choose to increase or decrease the number of samples. This allows the chart to re-render and display the number of columns based on the user`s selection. Figure below shows increased cases to 10000. Please note that any high arbitrary number can be selected but the chart will only show the maximum cases that GDC has.

[![Max Cases](./images/oncomatrix/25-max_cases.png)](./images/oncomatrix/25-max_cases.png 'Click to see the full image.')

The chart will reload with new cases added. Total cases now shown are 4878.

#### Group cases by

This option allows users to group cases by different variables from the GDC dictionary. Click on the `+` icon shown in blue to display different variables such as demographics, diagnoses, Exposures etc. Users may also search for a variable from the search bar provided in the menu as shown by `Search Variables` below.

[![Group cases](./images/oncomatrix/26-group_cases.png)](./images/oncomatrix/26-group_cases.png 'Click to see the full image.')

Click `Disease type` from the options. The matrix reloads to show the following view.

[![Grouping on disease type variable](./images/oncomatrix/27-grouping_disease_type.png)](./images/oncomatrix/27-grouping_disease_type.png 'Click to see the full image.')

As shown above, labels for different disease types show up vertically and all cases get distributed with a clearcut separation according to the type of disease.

Click on the `Disease type` (blue pill, as shown below). This opens a short menu with action items. Click on the first item `Edit` as shown below.

[![Group cases edit](./images/oncomatrix/28-group_cases_edit.png)](./images/oncomatrix/28-group_cases_edit.png 'Click to see the full image.')

Select the number of groups from the dropdown. User has the option to choose 2,3 or 4 groups. Drag types into different groups as shown below. In the figure below, Group 2 contains `Myeloid Leukemias` and `Nevi and Melanoma`.

Drag the type `Acinar Cell Neoplasm` into the box that shows the title `Excluded Categories`. After all the selections have been made, click `Apply`.

[![Groups and exclude](./images/oncomatrix/29-groups&exclude.png)](./images/oncomatrix/29-groups&exclude.png 'Click to see the full image.')

The matrix reloads with new groupings where `Myleoid Leukemias` and `Nevi and Melanoma` have been grouped together as Group 2, `Acinar Cell Neoplasm` is excluded from the groupings as well as the UI and Group 1 consists of the rest of the subtypes. The labels for the groups are user controlled and hence can be modified according to user requirements.

Click on the blue pill for `Disease type` and click `Cancel Grouping`. This will clear the groups and reload the default view with the types of diseases.

[![Cancel Grouping](./images/oncomatrix/30-cancel_grouping.png)](./images/oncomatrix/30-cancel_grouping.png 'Click to see the full image.')

The next option after edit is the `Reuse`. Click on this button to display the following.

[![Grouping Reuse option 1](./images/oncomatrix/31-grouping_reuse1.png)](./images/oncomatrix/31-grouping_reuse1.png 'Click to see the full image.')

Rename the title `Setting #1` as `My default setting` and click `Save`. The chart reloads. Click the `Disease type` option and click `Reuse` again from the menu. Figure below shows that user defined setting has been saved and applied as shown by a green checkmark.

[![Grouping Reuse option 2](./images/oncomatrix/32-grouping_reuse2.png)](./images/oncomatrix/32-grouping_reuse2.png 'Click to see the full image.')

This way a user can save multiple settings or delete them as required. Click the `Reuse` button again and click the `Use` button for the default setting. Once the default setting is applied `My default setting` will show a `Delete` button. Click this button to delete a setting.

[![Grouping Delete option](./images/oncomatrix/33-grouping_delete.png)](./images/oncomatrix/33-grouping_delete.png 'Click to see the full image.')

Click on the blue pill for `Disease type` again and click `Replace`. Select `Primary site` as shown below

[![Group cases replace option](./images/oncomatrix/34-group_cases_replace.png)](./images/oncomatrix/34-group_cases_replace.png 'Click to see the full image.')

The matrix reloads with the new variable distribution.

[![Group cases by Primary site variable](./images/oncomatrix/35-group_cases_primary_site.png)](./images/oncomatrix/35-group_cases_primary_site.png 'Click to see the full image.')

The last option on the menu is `Remove`. Click on the `4912 Cases` button, followed by `Primary site` shown in blue to reveal the menu option. Click `Remove` to completely get rid of any variable as shown below.

[![Remove Grouping](./images/oncomatrix/36-grouping_remove.png)](./images/oncomatrix/36-grouping_remove.png 'Click to see the full image.')

This will remove all and any groupings and show the default view again.

#### Sort Case Groups

Add the variable `Disease type` using the `Group Cases by` button as shown in the previous section. By default, groups are loaded ordered by their name. Change the selection to `Case count` as shown below. 

[![Sort grouped cases by case counts ](./images/oncomatrix/74-sort_case_counts.png)](./images/oncomatrix/74-sort_case_counts.png 'Click to see the full image.')

Now groups are ordered by the count of cases within the groups as shown above by the number of cases in parentheses in vertically shown group labels. The third selection option `Hits` orders the groupings based on the number of gene variants for a particular case or case group for the genes in display. Click `Hits` under `Sort Case Groups` to change the order of groupings as shown below. 

[![Sort grouped cases by case counts ](./images/oncomatrix/74-sort_case_counts.png)](./images/oncomatrix/74-sort_case_counts.png 'Click to see the full image.')

Hover over the first group label `Adenomas and Adenocarcinomas (763)` as shown below.

[![Hover over group labels ](./images/oncomatrix/76-Hover_group_label.png)](./images/oncomatrix/76-Hover_group_label.png 'Click to see the full image.')

This shows the number of cases in parenthesis of the group label and the breakdown for the  number of variants and CNV for all the samples within that group for the genes in display.

## Genes

The gene panel as shown below has several options as listed below for modifying the genes visible on the plot as well as their appearance/style.

- Display Case Counts for Gene
- Row Group Label Max Length
- Row Label Max Length
- Rendering Style
- Sort Genes
- Maximum # Genes
- Gene Set

[![Geneset options](./images/oncomatrix/37-genes_options.png)](./images/oncomatrix/37-genes_options.png 'Click to see the full image.')

### Display Case Counts for Gene

This option allows change in the number of cases that is represented in parentheses next to the gene variable label as shown below. By default, the number of cases for each gene is an `Absolute`.

[![Absolute case numbers](./images/oncomatrix/38-genes_absolute.png)](./images/oncomatrix/38-genes_absolute.png 'Click to see the full image.')

Click on the button `50 Genes` to display the menu and select `Percent` as shown below.

[![Percent cases 1](./images/oncomatrix/39-genes_percent1.png)](./images/oncomatrix/39-genes_percent1.png 'Click to see the full image.')

This shows the case counts as a percentage of the absolute values as shown below. TP53 has 100% case count that is, all 2000 cases show a mutation/CNV for TP53.

[![Percent cases 2](./images/oncomatrix/40-genes_percent2.png)](./images/oncomatrix/40-genes_percent2.png 'Click to see the full image.')

User has the option to hide the display of case counts. Click `50 Genes` button again and select `None` for `Display Case Counts for Gene` as shown below

[![None option for not displaying case numbers 1](./images/oncomatrix/41-genes_none1.png)](./images/oncomatrix/41-genes_none1.png 'Click to see the full image.')

This hides all the case counts as shown below.

[![None option for not displaying case numbers 2](./images/oncomatrix/42-genes_none2.png)](./images/oncomatrix/42-genes_none2.png 'Click to see the full image.')

### Rendering Style

The style of rendering for the sample cells/columns is an Oncoprint style by default. Click on `Stacked` option via `50 Genes` button as shown below.

[![Stacked rendering 1](./images/oncomatrix/43-rendering_stacked1.png)](./images/oncomatrix/43-rendering_stacked1.png 'Click to see the full image.')

The mutations and CNV are now stacked on top of each other as shown below.

[![Stacked rendering 2](./images/oncomatrix/44-rendering_stacked2.png)](./images/oncomatrix/44-rendering_stacked2.png 'Click to see the full image.')

### Sort Genes

The default sorting option for genes is `By Case Count`. This means the genes are sorted by the number of cases from increasing to decreasing order. Click `50 Genes` button on the control panel, and select `By Input Data Order` under the `Sort Genes` as shown below.

[![Sort gene options](./images/oncomatrix/45-sort_genes_input.png)](./images/oncomatrix/45-sort_genes_input.png 'Click to see the full image.')

The genes will now sort according to the order that is stored in the dataset and queried. However, please note that the sorting order can be overridden by the users choice as described in the section - Drag and Drop Gene Label/Variable.

### Maximum # Genes

The number of genes to display on the matrix plot can be modified by the input option as shown below. Click `50 Genes` button and change input number for `Maximum # Genes` to 70.

[![Maximun genes](./images/oncomatrix/46-max_genes.png)](./images/oncomatrix/46-max_genes.png 'Click to see the full image.')

The chart updates and loads the extra 20 genes. User can modify the set of genes by using the `Gene set` option next.

### Editing gene set

Gene groups can be edited using the `Gene set` option as shown below. Click `50 Genes` button to display this option and then click the `Edit` button in the `Gene set` as shown.

[![Geneset Edit 1](./images/oncomatrix/47-geneset_edit1.png)](./images/oncomatrix/47-geneset_edit1.png 'Click to see the full image.')

Select or deselect the blue checkbox to change the display to show `Cancer Gene Census` (CGC) genes only as shown below.

[![Geneset Edit 2](./images/oncomatrix/48-geneset_edit2.png)](./images/oncomatrix/48-geneset_edit2.png 'Click to see the full image.')

More information about CGC can be found at https://cancer.sanger.ac.uk/census. Figure displays the top 50 CGC genes.

User may choose to remove single genes one at a time by clicking over the genes shown above in gray.

Hover over TP53 as shown in the image below. A red cross mark appears with a description box. Click `TP53` to delete the gene as shown below.

[![Geneset Edit 3](./images/oncomatrix/49-geneset_edit3.png)](./images/oncomatrix/49-geneset_edit3.png 'Click to see the full image.')

User may choose to delete all genes from view by clicking the `Clear` button as shown below. However, a gene/variable selection is mandatory for the chart to load.

[![Clear all genes](./images/oncomatrix/50-clear_genes.png)](./images/oncomatrix/50-clear_genes.png 'Click to see the full image.')

After deleting the genes as required and making specific selections, click the `Submit` to update the matrix chart. In the figure below, several genes have been deleted such as TP53, KRAS, PTEN and MUC16 as shown. The matrix updates to reflect the deleted genes.

[![Delete genes from UI](./images/oncomatrix/51-del_genes.png)](./images/oncomatrix/51-del_genes.png 'Click to see the full image.')

### MSigDB genes

The MSigDB database (Human Molecular Signatures Database) has 33591 gene sets divided into 9 major collections and several subcollections. Users can choose to view the gene sets on the matrix plot.

Click on the `50 Genes` button. Then click on the `Gene set - Edit`. Here user can see a button with a dropdown for loading MSigDB genes.

[![MSigDb gene tree](./images/oncomatrix/52-msigdb1.png)](./images/oncomatrix/52-msigdb1.png 'Click to see the full image.')

Click on this dropdown to display a tree for the different gene sets.

[![MSigDb gene tree nodes](./images/oncomatrix/53-msigdb_tree.png)](./images/oncomatrix/53-msigdb_tree.png 'Click to see the full image.')

Select `C2: curated gene sets` and select `NABA_COLLAGENS` as shown below.

[![Naba Collagens](./images/oncomatrix/54-naba_collagens.png)](./images/oncomatrix/54-naba_collagens.png 'Click to see the full image.')

This loads the following genes as shown below.

[![Naba collagen geneset](./images/oncomatrix/55-msigdb2.png)](./images/oncomatrix/55-msigdb2.png 'Click to see the full image.')

Click `Submit` and the matrix will update to reflect the selected MSigDB gene set as shown below.

[![Naba collagen matrix](./images/oncomatrix/56-msigdb_collagen_matrix.png)](./images/oncomatrix/56-msigdb_collagen_matrix.png 'Click to see the full image.')

## Variables

The third button from the left called ``Variables` allows user to add in additional variables in the form of rows on the matrix. Click `Variables` to display a tree of variables and select `Disease type`, `Index date` and `Primary site`. Click the button `Submit 3 variables` as shown below.

[![Variable options](./images/oncomatrix/57-variables.png)](./images/oncomatrix/57-variables.png 'Click to see the full image.')

This updates the chart to display the selected variables on the very top of the matrix as shown below. User may choose to configure these rows by following steps outlined in section Clicking on gene/variable labels.

[![Variable matrix](./images/oncomatrix/58-variables_matrix.png)](./images/oncomatrix/58-variables_matrix.png 'Click to see the full image.')

User may also choose to add numerical variables. Click the `Variable` button from the control panel to display the following menu. Search and select the variable `Days to Birth` and hit `Submit 1 term` button as shown.

[![Selecting numerical variable](./images/oncomatrix/77-numerical_variable.png)](./images/oncomatrix/77-numerical_variable.png 'Click to see the full image.')

Once the matrix updates with the numerical variable on top, click on the label `Days to Birth` and click `Edit`. This shows the following menu. 

[![Selecting numerical variable](./images/oncomatrix/78-numerical_variable_continuous.png)](./images/oncomatrix/78-numerical_variable_continuous.png 'Click to see the full image.')

From this menu, select `Continuous` and hit `Apply`. This updates the matrix as shown below.

[![Selecting numerical variable](./images/oncomatrix/79-numerical_variable_continuous2.png)](./images/oncomatrix/79-numerical_variable_continuous2.png 'Click to see the full image.')

## Cell Layout

The cell layout menu enables customization of the appearance such as cell dimensions, spacing, font sizes, and borders. You may mouseover an input to see the description for that input, or try checking or editing inputs to test the effects of the control input and undo/redo as needed.

[![Cell Layout](./images/oncomatrix/59-cell_layout.png)](./images/oncomatrix/59-cell_layout.png 'Click to see the full image.')

## Legend Layout

The legend layout menu enables customization of the appearance of the legend, such as dimensions, spacing, and font sizes. These customizations can help avoid or minimize the need for post-download edits when generating figures.

[![Legend Layout](./images/oncomatrix/60-legend_layout.png)](./images/oncomatrix/60-legend_layout.png 'Click to see the full image.')

## Zooming

The matrix plot offers an interactive zoom panel as shown below with which a user can zoom in to view individual samples. There are two ways to use this panel. One by changing the input number and second by sliding the zoom bar to a desired zoom level as shown.

[![Zoom slider](./images/oncomatrix/61-zoom_slider.png)](./images/oncomatrix/61-zoom_slider.png 'Click to see the full image.')

Change zoom level to 10+ as shown.

[![Zoom slider plus](./images/oncomatrix/62-zoom_slider_plus.png)](./images/oncomatrix/62-zoom_slider_plus.png 'Click to see the full image.')

Scroll down to view individual samples at the bottom of the plot as shown below.

[![Zoomed matrix view](./images/oncomatrix/63-zoom_matrix_view.png)](./images/oncomatrix/63-zoom_matrix_view.png 'Click to see the full image.')

The zoom action can also be implemented by following steps as outlined in section - 'Drag to zoom'.

## Disco Plot

Click on any sample to reveal a second type of plot called as the `Disco Plot` as shown.

[![Disco plot 1](./images/oncomatrix/64-disco_plot1.png)](./images/oncomatrix/64-disco_plot1.png 'Click to see the full image.')

Click on `Disco plot` as shown above in gray. This loads a new chart above the matrix plot as shown below.

[![Disco plot 2](./images/oncomatrix/65-disco_plot2.png)](./images/oncomatrix/65-disco_plot2.png 'Click to see the full image.')

This plot shows all the mutations and CNV associated with that sample id as shown above. The plot also displays the legend for the mutation class and the CNV.

To reset the zoom level to default, click on the `Reset` button as shown. This will reset the zoom level to a default of 1.0

[![Zoom reset](./images/oncomatrix/66-reset_zoom.png)](./images/oncomatrix/66-reset_zoom.png 'Click to see the full image.')

## Undo/Redo/Restore

User may also choose to undo settings by clicking the `undo` button on the control panel as shown below.

[![Undo Zoom](./images/oncomatrix/67-undo_zoom.png)](./images/oncomatrix/67-undo_zoom.png 'Click to see the full image.')

Click `undo` to go back to the previous zoom level as shown above.

User may also choose to restore the initial state by clicking the `Restore` button as shown below.

[![Restore button](./images/oncomatrix/71-restore_btn.png)](./images/oncomatrix/71-restore_btn.png 'Click to see the full image.')

This restores the chart back to default settings as shown.

[![Restored UI](./images/oncomatrix/72-restored_chart.png)](./images/oncomatrix/72-restored_chart.png 'Click to see the full image.')


## Download

The control panel shows an option to download the plot as an svg after user has specified their customizations. Select the `Download` button as shown below to save the svg.

[![Download svg 1](./images/oncomatrix/68-download_svg1.png)](./images/oncomatrix/68-download_svg1.png 'Click to see the full image.')

The download will get saved to the default download folder as shown at the bottom of the browser window.

[![Download svg 2](./images/oncomatrix/69-download_svg2.png)](./images/oncomatrix/69-download_svg2.png 'Click to see the full image.')

## Legend

The legend for the matrix is below the plot and shows color coding for different mutation classes as well as color codes for CNV as shown here. (this is not interactive yet but will be in the near future).

[![Legend UI](./images/oncomatrix/70-legend.png)](./images/oncomatrix/70-legend.png 'Click to see the full image.')
