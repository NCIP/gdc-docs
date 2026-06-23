# OncoGrid

Oncogrid displays cases and mutated genes using a grid-style graphic.  Each column represents a case and each row represents a gene. The color and shape that are displayed at each intersection represent the presence/absence, type, and impact of mutations that are found in that gene for the associated case.  Mutation frequency and CNV frequency are both displayed alongside each row and column, while clinical attributes are displayed alongside each column.

[![OncoGrid](images/OncoGrid.png)](images/OncoGrid.png "Click to see the full image.")

A toolbar at the top right of the graphic allows the user to export the data as a JSON object, PNG image, or SVG image.  Seven buttons are available in this toolbar:

* __Customize Colors:__ Users can customize the colors that represent mutation consequence types and CNV gains/losses.
* __Download:__ Users can choose to export the contents either to a static image file (PNG or SVG format) or the underlying data in JSON format.
* __Reload Grid:__ Sets all OncoGrid rows, columns, and zoom levels back to their initial positions.
* __Cluster Data:__ Clusters the rows and columns to place mutated genes with the same cases and cases with the same mutated genes together. <font color=red> Some details suggested. Clicking the button doesn't make any visible change</font>
* __Toggle Heatmap View:__ The view can be toggled between cells representing mutation consequences or the number of mutations in each gene.
* __Toggle Gridlines:__ Turn the gridlines on and off.
* __Toggle Crosshairs:__ Turns crosshairs on, so that users can zoom into specific sections of the OncoGrid.
* __Fullscreen:__ Turns Fullscreen mode on/off.
