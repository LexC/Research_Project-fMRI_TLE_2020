# -*- coding: utf-8 -*-
"""
This script will unite the correlation matrices of the UF2C processed data and
it's main information. Beyond that some minimal calculations were made to
help on the organization and categorization.


The folders containing the files must fall the satandad structure:
    [...]\\Protocolo_xx\\Ys\\Processed\\Y_Y_x\\Y_Y_xx_xx\\Correlation_map
where x are single numbers and Y are strings. Exemple:
    [...]\\Protocolo_02\\Patients\\Processed\\PRS_ASL_12\\PRS_ASL_12_02\\
    Correlation_map

The information file must be a csv file with 4 colluns:
    Protocol, Type, Number, Injury Side
        Protocal=1 or 2
        Type = Patient or Volunteer
        Number = the subject #ID
        Injury Side = R,L,B,N,U and X for respectively Right, Left, Both,
                    Neither, Undefined and Not a TLE Patient
"""
# %% Libraries


import os
import scipy.io as sio
import csv
import numpy as np


# %% Variables and Inputs


fmriFoldersdir ='D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\Coleta_NIRS_fMRI_2015-2017\\Processed_data\\fMRI'
subjInfoFile = 'D:\\GD_UNICAMP\\IC_NeuroFisica\\Projetos\\fMRI_TLE_2020\\Variables_and_Data_Info\\\Info_Subjects\\subjectsInformation.csv'


# %% Classes


class corrMap:

    def __init__(self, matricesfilesdir, subjInfo):
        mat = sio.loadmat(matricesfilesdir)
        self.map = mat['map']

        matricesfilesdir = matricesfilesdir.split('\\')

        self.protocol = str(int(matricesfilesdir[-7].split('_')[1]))
        self.subjType = matricesfilesdir[-6][0:-1]
        self.subjNumber = str(int(matricesfilesdir[-3].split('_')[-2]))
        self.subjRun = str(int(matricesfilesdir[-3].split('_')[-1]))

        for i in subjInfo:
            if i[0:3] == [self.protocol, self.subjType, self.subjNumber]:
                self.injury = i[3]


# %% Functions


# Return a list with the files location of all relevant matrices
def getMatricesFilesDir(fmriFoldersdir):
    filesdir = []
    for root, subfolder, files in os.walk(fmriFoldersdir):
        for filename in files:
            if 'Protocolo_01\\Volunteers' not in root and filename == 'Matrix-VAR.mat':
                filesdir += [os.path.join(root, filename)]
    return filesdir


# Reads and save the subjects information
def subjectsInformation(subjInfoFile):
    subjInfo = []
    with open(subjInfoFile, 'r') as csv_subjinfo:
        for line in csv.reader(csv_subjinfo):
            subjInfo += [line]
    return subjInfo


# Make a variable with all the Correlation matrices and subject information
def loadSubjectsCorrelationMatrices(matricesfilesdir, subjInfo):
    corrmatrices = []
    for x in matricesfilesdir:
        corrmatrices += [corrMap(x, subjInfo)]
    return corrmatrices


# Fisher's z-Transformation
def zTransform(r):
    z = np.log((1 + r) / (1 - r)) * (1/2)
    return z


# Inverse Fisher's z Transformation
def inverse_zTransform(z):
#    r = (np.exp(2 * z) - 1)/(np.exp(2 * z) + 1)
    r = np.tanh(z)
    return r


# Group all the Correlation Matrices of a group in one 3D array. If injury is
# is not defined, it will group all matrices.
def GroupData_3Darray(CorrMats, injury):
    group = []
    if "injury" in locals():
        for x in CorrMats:
            if x.injury == injury:
                group += [x.map]
    else:
        for x in CorrMats:
            group += [x.map]
    return np.array(group)


# Return the mean correlation matrix
def fisherMean(array, dim):
    array = zTransform(array)
    array = np.mean(array, axis=dim)
    array = inverse_zTransform(array)
    return array


# %% Main


matricesfilesdir = getMatricesFilesDir(fmriFoldersdir)
subjInfo = subjectsInformation(subjInfoFile)
CorrMats = loadSubjectsCorrelationMatrices(matricesfilesdir, subjInfo)

CorrMats_RPatients = GroupData_3Darray(CorrMats, 'R')
CorrMats_LPatients = GroupData_3Darray(CorrMats, 'L')
#CorrMats_BPatients = GroupData_3Darray(CorrMats, 'B')
CorrMats_NPatients = GroupData_3Darray(CorrMats, 'N')
CorrMats_Volunteers = GroupData_3Darray(CorrMats, 'X')



























