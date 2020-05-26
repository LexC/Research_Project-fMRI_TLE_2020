# -*- coding: utf-8 -*-
"""
asdasdasd
"""
# %% Libraries

import os
import csv
import numpy as np
import scipy.io as sio
import nibabel as nib


# %% Variables Classes


class CorrMap: # pylint: disable=too-few-public-methods
    """ DESCRIÇÃO
    
    The folders containing the files must fall the satandad structure:
    [...]\\Protocolo_xx\\Ys\\Processed\\Y_Y_x\\Y_Y_xx_xx\\Correlation_map
where x are single numbers and Y are strings. Exemple:
    [...]\\Protocolo_02\\Patients\\Processed\\PRS_ASL_12\\PRS_ASL_12_02\\
    Correlation_map
    """
    
    def __init__(self, matricesfilesdir, subj_info):

        mat = sio.loadmat(matricesfilesdir)
        self.map = mat['map']

        matricesfilesdir = matricesfilesdir.split('\\')

        self.protocol = str(int(matricesfilesdir[-7].split('_')[1]))
        self.subj_type = matricesfilesdir[-6][0:-1]
        self.subj_number = str(int(matricesfilesdir[-3].split('_')[-2]))
        self.subj_run = str(int(matricesfilesdir[-3].split('_')[-1]))

        for i in subj_info:
            if i[0:3] == [self.protocol, self.subj_type, self.subj_number]:
                self.injury = i[3]


class GroupCorrMap: # pylint: disable=too-few-public-methods
    """ DESCRIÇÃO """

    def __init__(self, CorrMat, injuryClassification):

        self.group = CorrMat
        self.injury_side = injuryClassification


# %% Functions

####### Math


def ztransform(r):
    ''' Fisher's z-Transformation '''
    z = np.log((1 + r) / (1 - r)) * (1/2)
    return z


def inverse_ztransform(z):
    """ Inverse Fisher's z Transformation """
#    r = (np.exp(2 * z) - 1)/(np.exp(2 * z) + 1)
    r = np.tanh(z)
    return r


def fisher_mean(array, dim):
    """ Return the mean correlation matrix """
    array = ztransform(array)
    array = np.mean(array, axis=dim)
    array = inverse_ztransform(array)
    return array



####### Getting Data


def get_matrices_files_dir(fmri_foldersdir):
    """Return a list with the files location of all relevant matrices"""
    filesdir = []
    for root, subfolder, files in os.walk(fmri_foldersdir):
        for filename in files:
            if 'Protocolo_01\\Volunteers' not in root and filename == 'Matrix-VAR.mat':
                filesdir += [os.path.join(root, filename)]
    return filesdir


def subjects_information(subj_infofile):
    """Reads and save the subjects information
    
    The information file must be a csv file with 4 colluns:
    Protocol, Type, Number, Injury Side
        Protocal=1 or 2
        Type = Patient or Volunteer
        Number = the subject #ID
        Injury Side = R,L,B,N,U and X for respectively Right, Left, Both,
                    Neither, Undefined and Not a TLE Patient
    """
    subj_info = []
    with open(subj_infofile, 'r') as csv_subjinfo:
        for line in csv.reader(csv_subjinfo):
            subj_info += [line]
    return subj_info


def load_subjects_correlation_matrices(matricesfilesdir, subj_info):
    """ Make a variable with all the Correlation matrices and subject information """
    corrmatrices = []
    for x in matricesfilesdir:
        corrmatrices += [CorrMap(x, subj_info)]
    return corrmatrices



####### Grouping data and analyses: step 1


def group_data_3darray(corrmats, injury):
    """ Group all the Correlation Matrices of a group in one 3D array. If injury is not defined, it will group all matrices. """
    group = []
    if "injury" in locals():
        for x in corrmats:
            if x.injury == injury:
                group += [x.map]
    else:
        for x in corrmats:
            group += [x.map]
    return np.array(group)


def grouping_corrmaps(corrmats, injury_classifications):
    """ DESCRIÇÃO """
    corrmats_groups = []

    for  x in injury_classifications:
        corrmats_groups += [GroupCorrMap(group_data_3darray(corrmats, x), x)]
    return corrmats_groups


def group_mean(corrmats_groups):
    """ DESCRIÇÃO """
    for i, x in enumerate(corrmats_groups):
        corrmats_groups[i].group_mean = fisher_mean(x.group, 0)
    return corrmats_groups


def adjacency_matrices(matrix, threshold):
    """ DESCRIÇÃO """
    aux = np.array(matrix)
    aux[matrix >= threshold] = 1
    aux[matrix < threshold] = 0
    return aux


def group_mean_adjmats(corrmats_groups, thresholds):
    """ DESCRIÇÃO """
    for i, x in enumerate(corrmats_groups):
        corrmats_groups[i].adjMats = []
        for r in thresholds:
            corrmats_groups[i].adjMats += [adjacency_matrices(x.group_mean, r)]
    return corrmats_groups



####### Nifti functions (.nii)


def change_rois_values(data, network):
    """ DESCRIÇÃO """
    for i in range(1, int(data.max())+1):
        data[data == i] = network[i-1]
    return data


def make_nii_network(niifiledir, network, savename):
    """
    Change the ROI's IDvalue to the network's list values and save a .nii file
    with the new information.

    niifiledir = the location of a nifti file, with ROI'S enumerated from 1 to n.
    network = a list | len(network)=n
    savename = a string with the save name of the new nifti file
    """

    img = nib.load(niifiledir)

    new_data = change_rois_values(img.get_fdata(), network)

    img = nib.Nifti1Image(new_data, img.affine, img.header)
    nib.save(img, savename)

    return img
