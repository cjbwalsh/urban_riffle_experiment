# Constructed rock riffles increase habitat heterogeneity but not biodiversity in streams constrained by urban impacts: code and data.

*Christopher J Walsh, J Angus Webb, Daniel C Gwinn, and Peter F Breen*

This repository holds all of the code and smaller data used to produce the paper (currently in review).  The paper comprises for quarto documents:  

-  urban_riff_exp_ms.qmd; the main manuscript.  
-  urban_riff_exp_S1.qmd; appendix S1
-  urban_riff_exp_S2.qmd; appendix S2
-  urban_riff_exp_S3.qmd; appendix S3

Each of these should run and render from a cloned version of this repository. Each begins with code that checks for the existence of a data directory, and if absent, creates one and downloads three large data files from the linked [Open Science Framework Repository](https://osf.io/cms84/).

To run the Stan models used in the paper, the unevaluated code chunks (eval: false) in appendix S3 need to be run to create the large model object files.  


