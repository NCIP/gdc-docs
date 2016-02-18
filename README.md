# Docs Site installation


# Repository Conventions

- All Shared content in the "Commons" directory
- One Directory per GDC product (API, Data_Portal, Data_Submission_Portal, Data_Transfer_Tool)
- Each GDC product have a Users_Guide and Release_Notes directory

## Linking

To another documentation page
```
[Authentication and Authorization](../../Commons/Authentication.md)
```

Inside another documentation page

```
[Authentication and Authorization](../../Commons/Authentication.md#internal-section)
```

# Documentation Conventions

A detailed list of all conventions is available on [GDC Website](https://gdc.nci.nih.gov/conventions-page)

# Run mkdocs

# Build PDF

Install mkdocs2pandoc, following instructions available here:
```
https://github.com/jgrassler/mkdocs-pandoc
```

Prepare a yml file dedicated to your Userguide, using Data_Portal_UG.yml as an example.

Run the following commands to:
* Convert the User Guide to Pandoc:
* Tweak the pandoc file
* Build a PDF

```
mkdocs2pandoc -f Data_Portal_UG.yml -o docs/Data_Portal/PDF/Data_portal_UG.pd
sed -i -e 's/# / /g' docs/Data_Portal/PDF/Data_portal_UG.pd
sed -i -e 's/### /## /g' docs/Data_Portal/PDF/Data_portal_UG.pd
sed -i -e 's/\/site\//\/docs\//g' docs/Data_Portal/PDF/Data_portal_UG.pd
pandoc --toc -V documentclass=report -V geometry:"top=2cm, bottom=1.5cm, left=1cm, right=1cm" -f markdown+grid_tables+table_captions -o docs/Data_Portal/PDF/Data_portal_UG.pdf docs/Data_Portal/PDF/Data_portal_UG.pd
```
