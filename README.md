![](https://gdc.cancer.gov/sites/all/themes/gdc_bootstrap/logo.png)

# GDC Open Source code

=======
GDC is Open Source, Github Repositories containing source code of GDC Applications can be found on [GDC GitHub Organization page](https://github.com/NCI-GDC/).

- GDC Data Portal: https://github.com/NCI-GDC/portal-ui
- GDC Legacy Archive: https://github.com/NCI-GDC/portal-ui-legacy
- GDC Data Transfer Tool: https://github.com/NCI-GDC/gdc-client
- GDC Data Dictionary: https://github.com/NCI-GDC/gdcdictionary
- GDC Data Model: https://github.com/NCI-GDC/gdcdatamodel
- GDC Psqlgraph: https://github.com/NCI-GDC/psqlgraph 

# Support

Please direct technical questions to [GDC Support](https://gdc.cancer.gov/support).

# GDC Documentation Site

### Technology

 - Python 2.6, 2.7, 3.3, 3.4 and 3.5.
 - [mkdocs](http://www.mkdocs.org/)
 - [BSCodeTabs for mkdocs](https://github.com/mikecules/MarkdownBSCodeTabs#for-use-in-mkdocs)

### Install & Run

(Optional) Set up virtualenv:

- [Install virtualenv](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/)
- `python -m virtualenv venv`
- `source venv/bin/activate`
- Run the installation commands below
- To leave the virtual environment: `deactivate`

Install GDC-docs:

 - `pip install -r requirements.txt`
 - `mkdocs serve` (optionally set port `--dev-addr=0.0.0.0:<PORT>`)

### Build

 - `mkdocs build --clean`

### Repository Conventions

- All Shared content in the "Commons" directory
- One Directory per GDC product (API, Data_Portal, Data_Submission_Portal, Data_Transfer_Tool)
- Each GDC product have a Users_Guide and Release_Notes directory

### Linking

To another documentation page
```
[Authentication and Authorization](../../Commons/Authentication.md)
```

Inside another documentation page

```
[Authentication and Authorization](../../Commons/Authentication.md#internal-section)
```

### Adding icons and PDFs
The convention for this, when updating mkdocs.yml is the following:
- <font-awesome-icon> <content> <url ending in .pdf>: 'index.md'
example:
- fa-file-pdf-o Download PDF /API/PDF/API_UG.pdf: 'index.md'

### Documentation Conventions

A detailed list of all conventions is available on [GDC Website](https://gdc.cancer.gov/conventions-page)

## Material Port

Note for a port to the Material theme for mkdocs.

### Resource

[Material for MkDocs documentation](https://squidfunk.github.io/mkdocs-material/getting-started/)

### Python Version

Use Python 3.8 or greater. For development, this branch is using Python 3.12.

### Installation

Create virtual environment

```bash
python -m venv venv
```

Activate virtual environment

```bash
source venv/bin/activate
```

Pre-req: Make sure `pip-tools` is installed

```bash
pip install pip-tools
```

Install MkDocs and dependencies

```bash
pip-sync requirements.txt
```

### Run Locally

For development, start a mkdocs server locally with

```bash
mkdocs serve
```

In a browser, open `http://127.0.0.1:8000/`

### Generate Site

To build the site, run

```bash
mkdocs build
```

### Building User Guide PDFs

Building PDFs requires a library named **Pango** which used for laying out and rendering text. Install it using apt. 

```bash
sudo apt install libpango1.0-dev
```

The `mkdocs-with-pdf` plugin is used to generate PDFs.  

```bash
ENABLE_PDF_EXPORT=1 mkdocs build -f API_UG.yml 
ENABLE_PDF_EXPORT=1 mkdocs build -f Data_Portal_UG.yml 
ENABLE_PDF_EXPORT=1 mkdocs build -f Data_Submission_Portal_UG.yml 
ENABLE_PDF_EXPORT=1 mkdocs build -f Data_Transfer_Tool_UG.yml 
ENABLE_PDF_EXPORT=1 mkdocs build -f Data_UG.yml 
```

Build the site again to move the generated PDFs in the correct directory.

```bash
mkdocs build
```

### Known Issues

- Dict search app does not work
- Dict viewer app does not work
- Code blocks not displaying for multiple languages
- PDF generation needs to be added
  
### Build Logs

The build logs give us many warnings that could be investigated and fixed. Categories of messages:

- The following pages exist in the docs directory, but are not included in the "nav" configuration
- "there is no such anchor on this page"
- "contains an absolute link"
- "not found among documentation files."
- "contains an unrecognized relative link"
- "does not contain an anchor"
  
