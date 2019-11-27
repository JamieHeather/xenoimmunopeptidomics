# Xeno-Immunopeptidomics

#### James M Heather *et al*, 2019, the Cobbold lab @ MGH

This repository contains data and scripts designed to test the hypothesis that growing cell lines as xenografts in immunodeficient mice is a valid replacement for cell culture for the production of material for immunopeptidomics. This is achieved by examining the overlap between different HLA ligandome datasets (both peptidomes and phosphopeptidomes) from the JY and other cell lines, grown either in mice or in flasks.

This project was done primarily in collaboration with [Agenus Inc.](https://agenusbio.com/), who performed all of the JY mass spectrometry and primary analysis as described in the paper this repo accompanies. As Agenus are [pursuing phosphopeptide antigens as targets for cancer immunotherapy](http://agenusbio.com/wp-content/uploads/2018/11/SITC_CRCPSV_Poster_PTM_version3.pdf) the IMAC-enriched phosphopeptide data has not been included in the uploaded/shared version of this repo. Similarly a small number of phosphopeptides (not currently covered by appropriate IP measures) have been redacted from the shared **non**-IMAC-enriched data which is being shared. However these minor redactions make no discernable difference to the outcome of the analyses undertaken here, and we believe that the lack of the phosphopeptide enriched data does not detract from the overall findings of the work.    

Confirmation of the utility of mice to recapitulate the phospho-immunopeptidome in other cell lines was performed in collaboration with the [Hunt lab at UVA](https://med.virginia.edu/faculty/faculty-listing/dfh/).

## Other JY data sources

JY immunopeptidome datasets published by other groups have been collated here for comparison to our cultured and mouse results.

Note that there were minor edits to downloaded files prior to processing:
* All xls files were resaved as xslx (for easier parsing with the `xlrd` module), and any additional non-data/non-header rows were moved to other worksheets. 
    * Extraneous columns/worksheets have also been deleted, to reduce file size 
* All csv files downloaded from SysteMHC were first modified using the following one liner to remove quote marks:

````
    for f in *csv; do echo $f; sed -i 's/\"//g' $f; done
````

* [Bassani-Sternberg *et al.*, 2015 Mol. Cell Proteomics](http://dx.doi.org/10.1074/mcp.M114.042812)
* [Hassan *et al.*, 2013 Mol. Cell Proteomics](http://dx.doi.org/10.1074/mcp.M112.024810)
    * Note that this raw data was actually presented in one spreadsheet, with one list of peptides and intensity values for both cell types
    * For convenience of analysis, the file was duplicated and the list of peptides trimmed down such that each version only has the peptides found on that line
* [Bourdetsky *et al.*, 2014 PNAS](http://dx.doi.org/10.1073/pnas.1321902111)
* [Caron *et al.*, 2015 eLife](http://dx.doi.org/10.7554/eLife.07661)

### Non-JY control datasets

Downloaded from [SysteMHC](https://systemhcatlas.org/datasets), accessed in January 2019.
* MM16 - SYSMHC00023
    * File labelled '1_A'
* THP1 - SYSMHC00002
* C1866 - SYSMHC00011
    * File labelled 'neg_all'

# Running the scripts

All scripts are designed to run on Python 2.7, tested on 2.7.12/2.7.15 on both Mac and Linux OSs.

The scripts can be run individually, by navigating to the 'Scripts' directory and running the appropriate ```python [script-name].py``` command, after ensuring the appropriate dependencies are installed. The plots produced should appear in a dated folder in the 'Plots' directory, along with certain other data files in the 'Data' directory. Some files need to be run in order, as they make use of data files produced by previous scripts (see below).

Alternatively you can run ```python run-analysis.py``` in the same directory, which should run each individual script in turn. The code for this script contains a suggested run order for the scripts to ensure proper intermediate file production. 

Note that an internet connection is required for some scripts, as certain larger data files are downloaded from the internet.

Also note that the 'netmhc-check.py' script will need to have NetMHC installed and run on the data prior to running, in order to generate the input data, and the 'phospho-analysis.py' will not run as its raw data has not been shared.

# Dependencies

Versions provided where available.

* scipy.__version__ '1.1.0'
* numpy.__version__ '1.15.3'
* matplotlib.__version__ '2.2.3'
* pandas.__version__ '0.23.4'
* seaborn.__version__ '0.9.0'
* upsetplot.__version__ '0.1'
* mygene.__version__ '3.0.0'
* mhcflurry.__version__ '1.2.2'
* matplotlib_venn 
* acora

All maintained using **anaconda2** (4.5.11) installing via conda where available, and pip (18.1) for everything else.
