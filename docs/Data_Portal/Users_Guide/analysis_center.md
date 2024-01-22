# GDC Analysis Center


The Analysis Center is the central hub for accessing tools to support cohort analysis. The Analysis Center can be accessed by clicking on the Analysis Center button in the GDC Data Portal header, on the "Explore Our Cancer Datasets" button on the home page, or on one of the sites in the human anatomical outline or bar graph.

[![Analysis Center View](images/FullAnalysisCenter.png)](images/FullAnalysisCenter.png "Click to see the full image.")

The Analysis Center consists of a main toolbar and a query expressions section, both of which are always displayed. The main toolbar displays the active cohort and can be used to create and manage custom cohorts. The query expression section displays the filters applied to the active cohort.

Available tools are displayed under the Query Expression section of the Analysis Center. Each analysis tool is showcased within a tool 'card', which has several items related to the analysis tool such as:

* A teal 'Play' button to launch the analysis tool on the given cohort
* A 'Demo' button that launches a demonstration of the analysis tool on an example cohort
* Clicking on the name of the analysis tool in the tool card toggles a drop down description of the analysis tool
* The number of cases from the cohort that the analysis will be performed on is at the bottom of the card

## Core Tools ##

This section contains the core GDC tools.  This includes the [Projects](Projects.md) tool, the [Cohort Builder](cohort_builder.md), and the [Repository](Repository.md). These can be selected for use with your current cohort by clicking on each from the Core Tools section.

## Analysis Tools ##

The Analysis Tools section contains the analysis tools available in the Analysis Center: [BAM Slicing Download](BAMslicing.md), [Clinical Data Analysis](clinical_data_analysis.md), [Cohort Comparison](cohort_comparison.md), [Gene Expression Clustering](LINKHERE), [Mutation Frequency](mutation_frequency.md), [OncoMatrix](oncomatrix.md), [ProteinPaint](proteinpaint_lollipop.md), [Sequence Reads](proteinpaint_bam.md), and [Set Operations](set_operations.md).

These can be launched by clicking the Play buttons in each of the tool cards.

If there is not sufficient data in the active cohort to use a particular tool, the play button will be grayed out and will not be usable until a new cohort with sufficient data is selected.

[![Analysis Center Tools](images/AnalysisCenterTools.png)](images/AnalysisCenterTool.png "Click to see the full image.")

## Tool Panel

As each tool is selected, it is loaded in the Analysis Center within a panel.

[![Analysis Center Tools Panel](images/AnalysisCenterToolPanel.png)](images/AnalysisCenterToolPanel.png "Click to see the full image.")

To close a tool and return to the default view that displays all the tool cards within the Analysis Center, click the "X" to the left of the tool's header.
