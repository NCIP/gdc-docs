# Data Release Notes






## Initial Data Release

* __GDC Product__: GDC Data Release
* __Release Date__: May 10, 2016

### New Datasets

* TCGA
* TARGET

### Known Issues and Workarounds

BAM files produced by our RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.

Details are provided in [Data Release Manifest](Manifests/GDC_PreLaunch_release_notes_05102016.txt)
