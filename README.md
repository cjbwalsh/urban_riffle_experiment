# Constructed rock riffles increase habitat heterogeneity but not biodiversity in streams constrained by urban impacts: code and data.

*Christopher J Walsh, J Angus Webb, Daniel C Gwinn, and Peter F Breen*

This repository holds all of the code and smaller data used to produce the paper (currently in review).  The paper comprises for quarto documents:  

-  *urban_riff_exp_ms.qmd*; The manuscript.  

-  *urban_riff_exp_S1.qmd*; appendix S1 (Study site details)

-  *urban_riff_exp_S2.qmd*; appendix S2 (Notes on biological data preparation)

-  *urban_riff_exp_S3.qmd*; appendix S3 (Methods and code for fitting and assessing the models)


Each of these should run and render from a cloned version of this repository (tested in RStudio), allowing you to render all the text, figures and tables of the documents. The manuscript, S1 and S3 documents with code that checks for a data directory, and if necessary, creates the directory and downloads four large data files from the (Open Science Framework repository)[https://osf.io/cms84/] into it.

The figures and tables in the manuscript are constructed from small tables derived from the Stan models described in S3. To run the Stan models used in the paper, the unevaluated code chunks (eval: false) in appendix S3 need to be run to create the large model object files.  

