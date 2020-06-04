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


####### Calculations
# thresholds = list(np.arange(0.05, 1, 0.05))
thresholds = list(np.arange(0.05, 1, 0.05))
# thresholds = [0.3,0.5,0.9]

####### Graph Theory

MNI_COOR_DIR = "D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Programing\\Variables_and_Data_Info\\fMRI_Processing_Var\\Shen2013_mnicoor.csv"

GRAPHS_SAVES_DIR1 = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Graphs\\Step1'

GRAPHS_SAVES_DIR2 = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Graphs\\Step2'


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

coormats = pf.CoorMatrices(matricesfilesdir)
coormats = coormats.get_raw_data(subjInfo)


# # Separeting Groups and Calculating Group Mean and Adjacency Matrices of the Correlation Matrices
coormats = coormats.get_group_by_subject(thresholds)

# # Generate Brain Graphs images as a function of the thresholds
# pf.graphs_for_all(coormats, MNI_COOR_DIR, GRAPHS_SAVES_DIR1, thresholds,69)

# # Trnasforming Adjacency Matrices into Graphs
G = pf.LexGraphs(coormats.group_by_subject)

# # Graphs Parameters
G = G.get_average_degrees()

# # Graphs Results
pf.gfigure_parameter_vs_threshold(G.average_degrees, G.thresholds, coormats.group_by_subject, GRAPHS_SAVES_DIR2, 'Average_K')


# %%

