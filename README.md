flowcytometry_in_R

Script for standard analysis of FACS data, later to be included in more generic functions. Stared on 16-October-2018.

User Input Section: Directories specify locations of raw data, gating template, and saving locations of plots created. Idea is for this section to be the only place the user has to specify their variables, from which point onwards the script will run automatically.

Note:

reading of .fcs files uses regex which may not fit all file names, user has to specify these
'markers' is a .txt matrix with raw data .fcs file colnames in first column and desired names in second column
gating template needs to be a .csv file in this format: https://bioconductor.org/packages/devel/bioc/vignettes/openCyto/inst/doc/HowToWriteCSVTemplate.html
See the following link for further details: https://www.bioconductor.org/packages/release/bioc/vignettes/openCyto/inst/doc/openCytoVignette.html

See session info tab for citation of all packages used in the flow cytometry and statistical analysis.
