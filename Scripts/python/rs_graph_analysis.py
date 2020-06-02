# -*- coding: utf-8 -*-
"""
This script will:
    * Unite the correlation matrices of the UF2C processed data and it's main information.
    * Organize, categorize and group the data.
    * Calculate Group Average Correlation Matrices
    * Make Ajacency Matrices
"""

# %% Libraries

import numpy as np
import fmri_tle_2020_tools as pf # pf for "Project Functions"


# %% Variables and loads

####### Correlation Matrices and grouping
FMRI_FOLDERS_DIR = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\Coleta_NIRS_fMRI_2015-2017\\Processed_data\\fMRI'

SUBJ_INFO_FILE = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Programing\\Variables_and_Data_Info\\Info_Subjects\\subjectsInformation.csv'

injuryClassifications = ['R', 'L', 'N', 'X']

thresholds = list(np.arange(0.05, 0.95, 0.05))

####### Seed-Based Variables

# This set of rois is in the MNI152 template.
NIIBASEFILE = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Programing\\Variables_and_Data_Info\\fMRI_Processing_Var\\Shen2013.nii'

SEEDB_RESULTS_FOLDER = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Seed_based_analysis'



# %% Main

# # Cronstructing correlation matrices variable
matricesfilesdir = pf.get_matrices_files_dir(FMRI_FOLDERS_DIR)
subjInfo = pf.subjects_information(SUBJ_INFO_FILE)
CorrMats = pf.load_subjects_correlation_matrices(matricesfilesdir, subjInfo)


# # Separeting Groups and creating
CorrMats_Groups = pf.grouping_corrmaps(CorrMats, injuryClassifications)


# # Calculating Group Mean and Adjacency Matrices of the Correlation Matrices
CorrMats_Groups = pf.group_parameters(CorrMats_Groups, thresholds)

