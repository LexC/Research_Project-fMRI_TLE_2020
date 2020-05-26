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
import fmri_tle_2020_tools as pf # pf for Project Functions


# %% Variables and Inputs


FMRI_FOLDERS_DIR = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\Coleta_NIRS_fMRI_2015-2017\\Processed_data\\fMRI'

SUBJ_INFO_FILE = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Programing\\Variables_and_Data_Info\\Info_Subjects\\subjectsInformation.csv'

injuryClassifications = ['R', 'L', 'N', 'X']

thresholds = list(np.arange(0.05, 0.95, 0.05))

# %% Main

# # Cronstruct correlation matrices variable
matricesfilesdir = pf.get_matrices_files_dir(FMRI_FOLDERS_DIR)
subjInfo = pf.subjects_information(SUBJ_INFO_FILE)
CorrMats = pf.load_subjects_correlation_matrices(matricesfilesdir, subjInfo)


# # Separeting Groups and creating
CorrMats_Groups = pf.grouping_corrmaps(CorrMats, injuryClassifications)


# # Calculating Group Mean of the Correlation Matrices and Adjacency Matrices
CorrMats_Groups = pf.group_mean(CorrMats_Groups)
CorrMats_Groups = pf.group_mean_adjmats(CorrMats_Groups, thresholds)


# %% TESTS 
"""
Essa é a ideia de como fazer pra gerar as redes em .nii, agora o que falta é definir as redes/ sementes. 
"""

NIIBASEFILE = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Programing\\Variables_and_Data_Info\\fMRI_Processing_Var\\Shen2013.nii'

NIISAVENAME = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Programing\\Variables_and_Data_Info\\fMRI_Processing_Var\\TEST.nii'

SEED = 120

network=CorrMats_Groups[0].adjMats[thresholds.index(0.3)][SEED,:]*100


pf.make_nii_network(NIIBASEFILE, network, NIISAVENAME)
