## **Overview**<br>
This repository contains code associated with "Functional analysis of O-GlcNAcylation by networking of OGT interactors and substrates" by Griffin ME, Thompson JW, et al.

It includes our novel algorithms for parsimonious enumeration of PTMs in mass spectrometry data (see Extended Data Figure 2. Main Text, and Methods):<br>
1. SiteAndRegions.py<br>
2. SitesAndRegionsMultiExperiment.py (includes the additional option to separate conditions/replicates)<br>

It also includes Jupyter Notebooks:<br>
1. Adaptor Protein Ranking.ipynb (used to rank potential hub proteins)<br>
2. Network Generation.ipynb (used to make OGT substrate and interactor networks) <br>
3. Sites and Regions Analysis.ipynb (used to quantify and localize O-GlcNAc sites with the algorithm above)<br>

Finally, sample input and output data is also provided:<br>
1. 293T_Glycomics_Full_PSMs.txt (Proteome Discoverer PSMs file from our O-GlcNAcomics experiments in HEK 293T cells)<br>
2. 293T_Preprocessed.txt (sample input derived from '293T_Glycomics_Full_PSMs.txt' for 'SiteAndRegions.py' as described in 'Sites and Regions Analysis.ipynb')<br>
3. 293T_Preprocessed_2.txt (sample input derived from '293T_Glycomics_Full_PSMs.txt' for 'SitesAndRegionsMultiExperiment.py' as described in 'Sites and Regions Analysis.ipynb')<br>

## **System Requirements**<br>
All code has been verified to run in Python v3.7.3 with all Anaconda (v4.8.3) packages installed. It has been tested in both Microsoft Windows 10 (v22H2) and Linux (Ubuntu v20.04).

## **Run Time**<br>
Install time for Python and Anaconda varies based on network speed. Run time for all code above is typically <1 minute depending on the size of the input data. The entire example notebooks above can be typically run using the input data from this manuscript in <5 minutes.

## **License**<br>
This project is covered under the MIT License
