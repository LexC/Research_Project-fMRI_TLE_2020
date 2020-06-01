# Network Study of Patients with Temporal Lobe Epilepsy: A fMRI Analyses

## About

This Repository contains self-made scripts (majorly in Python) to perform graph theory analyses, and data visualization, on fMRI data of patients with Temporal Lobe Epilepsy (TLE) and healthy volunteers.


## Quick Introduction

On the [Biomedical Optics Laboratory  (LOB)](http://light4brain.com/), of the [University of Campinas](http://www.internationaloffice.unicamp.br/hotsite/), we collected dozens of Structurals and Functionals Magnetic Resonance Imaging (MRI) scans patients with left TLE (LTLE), right TLE (RTLE) or non-lateral TLE (NTLE) and healthy volunteers, in the resting-state (RS). After that, the images were reoriented, pre-processed and processed with the software [SPM12](https://www.fil.ion.ucl.ac.uk/spm/) and the toolbox [UF<sup>2</sup>C](https://www.lniunicamp.com/uf2c) (with Adds-on made by the LOB's team). Each brain image was parcellated in regions of interest (ROI) based on the [Shen Atlas](https://doi.org/10.1016/j.neuroimage.2013.05.081) and created a Pearson Correlation Matrices.

The scripts in this repository act on the Correlation Matrices, to perform individual and group analyses, employing seed-based technique and graphs theory. 


## A teaser of some Results

### Seed-Based thechinique



![Image of DMNs](https://github.com/LexC/Research_Project-fMRI_TLE_2020/blob/master/Images/DMN.gif?raw=true)
FIgure 1: The Default Mode Netwotk (DMN) of the Volunteers, LTLE and RTLE group analyses. On the gradient, the Pearson Correlation values of each ROI (on this image it's only shown ROIs with correlation values greater then 0.2), and on *blue* the ROI-seed used to the analyses.

<p>&nbsp;</p>

<p>&nbsp;</p>

<p>&nbsp;</p>


#### Disclamer: No personal information was disclosed
