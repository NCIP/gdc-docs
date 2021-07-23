RPPA

**R**everse **P**hase **P**rotein **A**rray (RPPA), originally developed by [MD Anderson Cancer Center](https://www.mdanderson.org/), University of Texas, is a high-throughput antibody-based technique with the procedures similar to that of Western blots. Proteins are extracted from tumor tissue or cultured cells, denatured by SDS, printed on nitrocellulose-coated slides followed by antibody probe. A group (approx. several hundreds) of antibodies form a *set*, which will be used for each assay. Antibodies may be added to or removed from the set depending on feasibility/functionality

To quantify protein expression, a “standard curve” is constructed from spots on each slide (one slide probed for one antibody). These spots include serial dilutions of each sample plus QC spots of standard lysates at different concentrations.

The technique is capable of the following types of analysis:

* Classify patient tumors
* Correlate DNA, RNA and Protein
* Determine prognosis
* Predict responses to targeted therapies
* Pharmacodynamics and biologically relevant dose
* Determine appropriate handling procedures for clinical samples (based on antigen stability analysis)

The GDC RPPA dataset contains the following fields:

* ```AGID```: The antigen ID
* ```lab_id```: The unique antibody ID
* ```catalog_number```: Antibody vendor's catalog number 
* ```set_id```: The ID for a set, ie list of antibodies (eg refs [3] & [4]).
* ```peptide_target```: The ID for the target site which antigen binds to.
* ```protein expression```: Relative levels of protein expression - interpolation of each dilution curve to the “standard curve” (supercurve) of the slide (antibody).

**References:**

1. https://bioinformatics.mdanderson.org/public-software/tcpa/

2. https://www.mdanderson.org/research/research-resources/core-facilities/functional-proteomics-rppa-core/rppa-process.html

3. https://www.mdanderson.org/content/dam/mdanderson/documents/core-facilities/Functional%20Proteomics%20RPPA%20Core%20Facility/RPPA_Expanded_Ab_List_Updated.xlsx

4. https://www.mdanderson.org/content/dam/mdanderson/documents/core-facilities/Functional%20Proteomics%20RPPA%20Core%20Facility/RPPA_Standard_Ab_List_Updated.xlsx

