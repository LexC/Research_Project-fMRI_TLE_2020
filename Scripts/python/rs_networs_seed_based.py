# -*- coding: utf-8 -*-
"""
This script will:
    * Unite the correlation matrices of the UF2C processed data and it's main information.
    * Organize, categorize and group the data.
    * Calculate Group Average Correlation Matrices
    * Make Ajacency Matrices
    * Create and save function networks
"""

# %% Libraries

import numpy as np
import fmri_tle_2020_tools as pf # pf for "Project Functions"


# %% Variables and loads

####### Correlation Matrices and grouping
FMRI_FOLDERS_DIR = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\Coleta_NIRS_fMRI_2015-2017\\Processed_data\\fMRI'

SUBJ_INFO_FILE = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Programing\\Variables_and_Data_Info\\Info_Subjects\\subjectsInformation.csv'

thresholds = list(np.arange(0.05, 0.95, 0.05))

####### Seed-Based Variables

# This set of rois is in the MNI152 template.
NIIBASEFILE = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Programing\\Variables_and_Data_Info\\fMRI_Processing_Var\\Shen2013.nii'

SEEDB_RESULTS_FOLDER = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Seed_based_analysis\\TEST'

networks_names = ['DMN', 'SMOTOR', 'VISUAL', 'AUDIT', 'DAN_L', 'DAN_R', 'ECN_L', 'ECN_R']
seeds_mniloc = [[-7, -57, 18], [-25, -12, 64], [-3, -85, -3], [63, -12, 8], [-22, -53, 56], [27, -54, 55], [-42, 34, 20], [44, 36, 20]]


# %% Main

# # Cronstructing correlation matrices variable
matricesfilesdir = pf.get_matrices_files_dir(FMRI_FOLDERS_DIR)
subjInfo = pf.get_csv_information(SUBJ_INFO_FILE)
"""Reads and save subjects informations

subj_infofile = str(folder diretory)

The information file must be a csv file with 4 colluns:
Protocol, Type, Number, Injury Side
    Protocal=1 or 2
    Type = Patient or Volunteer
    Number = the subject #ID
    Injury Side = R,L,B,N,U and X for respectively Right, Left, Both, Neither, Undefined and Not a TLE Patient
"""

CorrMats = pf.CoorMatrices(matricesfilesdir, subjInfo)


# # Separeting Groups and Calculating Group Mean and Adjacency Matrices of the Correlation Matrices
CorrMats = CorrMats.get_group_by_injury(thresholds)


# # Seed-based technique

# # # Finding ROI's seeds
rois_seeds = pf.roi_location(seeds_mniloc, NIIBASEFILE, 'mni')

# # # Calculating and saving networks of Group Averages
# pf.make_seed_based_networks(CorrMats.group_by_injury, rois_seeds, 0.3, networks_names, NIIBASEFILE, SEEDB_RESULTS_FOLDER)
