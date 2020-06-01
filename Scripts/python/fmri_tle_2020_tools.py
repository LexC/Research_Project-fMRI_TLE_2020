# -*- coding: utf-8 -*-
"""

"""
# %% Libraries

import os
import csv
import numpy as np
import scipy.io as sio
import nibabel as nib


# %% Classes

class CorrMap: # pylint: disable=too-few-public-methods
    """ A Class to load and arquive the fmri correlation matrices and relevant infromations.

    The folders containing the files must fall the satandad structure:
    [...]\\Protocolo_xx\\Ys\\Processed\\Y_Y_x\\Y_Y_xx_xx\\Correlation_map
where x are single numbers and Y are strings. Exemple:
    [...]\\Protocolo_02\\Patients\\Processed\\PRS_ASL_12\\PRS_ASL_12_02\\
    Correlation_map
    """

    def __init__(self, matricesfilesdir, subj_info):
        """
        matricesfilesdir = str(matrix diretory)
        subj_info = list(subjects informations)

        Output:
            self.map = np.array( Correlation Matrices)
            self.protocol = int()
            self.subj_type = str(Normaly Patient or Volunteer)
            self.subj_number = int()
            self.subj_run = int()
            self.injury = str()
        """
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
    """ A class to group the correlatin Matrices"""

    def __init__(self, CorrMat, injuryClassification):


        self.group = group_data_3darray(CorrMat, injuryClassification)
        self.injury_side = injuryClassification

    def group_mean (self):

        self.group_mean = fisher_mean(self.group, 0)

    def group_mean_adjmats (self,thresholds):

        self.adjmats = []
        for r in thresholds:
            self.adjmats += [adjacency_matrices(self.group_mean, r)]

# %% Functions

####### Math ------------------------------------


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


def transform_coor_space(cor, transform_matrix,transform_direction):
    """
    Transform coordenates space

    coor = list[0:2] -> coordenates
    transform_matrix = nympy.ndarry(), .shape =(4, 4) -> Transformation matrix
    transform_direction = 0 or 1 - > 0 for direct transformation, 1 for inverse transformation.

    output: list[0:2]
 """
    cor = np.array(cor)

    corfix = np.array(np.append(cor, 1))


    #mni2coor
    if transform_direction == 0:
        new_cor = corfix.dot(np.transpose(np.linalg.inv(transform_matrix)))
    #cor2mni
    else:
        new_cor = np.transpose(transform_matrix.dot(np.transpode(corfix)))

    return np.array(np.round(new_cor[0:3]), dtype=int)



####### Getting Data ------------------------------------


def get_matrices_files_dir(fmri_foldersdir):
    """Return a list with the files location of all relevant matrices

    fmri_foldersdir = str(folder diretory)
    """
    filesdir = []
    for root, subfolder, files in os.walk(fmri_foldersdir):
        for filename in files:
            if 'Protocolo_01\\Volunteers' not in root and filename == 'Matrix-VAR.mat':
                filesdir += [os.path.join(root, filename)]
    return filesdir


def subjects_information(subj_infofile):
    """Reads and save subjects informations

    subj_infofile = str(folder diretory)

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
    """ Make a variable with all the Correlation matrices and subject information

    matricesfilesdir = list(matrices diretories), from get_matrices_files_dir
    subj_info = list(subjects informations), from subjects_information
    """
    corrmatrices = []
    for x in matricesfilesdir:
        corrmatrices += [CorrMap(x, subj_info)]
    return corrmatrices



####### Grouping data and analyses ------------------------------------


def group_data_3darray(corrmats, injury):
    """
    Group all the Correlation Matrices of a group in one 3D array. If injury is not defined, it will group all matrices.

    corrmats = list(class CorrMap), from load_subjects_correlation_matrices
    injury = str()

    Outuput: np.array(), shape =(3,:,:)
    """
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
    """
    corrmats = list(class CorrMap), from load_subjects_correlation_matrices
    injury_classifications = list(str())

    Output: list (class GroupCorrMap)
    """
    corrmats_groups = []

    for  x in injury_classifications:
        # corrmats_groups += [GroupCorrMap(group_data_3darray(corrmats, x), x)]
        corrmats_groups += [GroupCorrMap(corrmats, x)]


    return corrmats_groups


def adjacency_matrices(matrix, threshold):
    """ Returns an adjacency matrix based on a threshold"""
    adjmtx = np.array(matrix)
    adjmtx[matrix >= threshold] = 1
    adjmtx[matrix < threshold] = 0
    return adjmtx



def group_parameters(corrmats_groups, thresholds):
    """
    Return the mean correlation matrix of each group from the list and adjacency matrices based on a range of threshold

    corrmats_groups = list(class GroupCorrMap)
    thresholds = list(int())

    Output: list(class GroupCorrMap)
        added structures:
            .group_mean
            .adjmats
    """
    for x in corrmats_groups:
        x.group_mean()
        x.group_mean_adjmats(thresholds)

    return corrmats_groups



####### Nifti functions (.nii) ------------------------------------



def change_rois_values(data, network):
    """ Alters the the roi's valeus (roi's IDs) to the values of the network, Assuming IDs are continuous integers and ranging from 0."""
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


def roi_location(coor, niifiledir, coortype):
    """
    Return the ROI's number value (ID). The nifti file must be in agreement with the coordinates.

    coor = list[n][0:2] -> coordenates
    niifiledir = str() -> diretory of the nii file with ROIs informations
    coortype = 'mni' or 'cor'

    Output: list(int())
    """

    img = nib.load(niifiledir)
    mtx = img.get_fdata()
    coor=np.array(coor)

    if coortype == 'mni':
        #coor = mni2cor(coor,niifiledir)
        for i,x in enumerate(coor):
            coor[i] = transform_coor_space(x, img.affine,0)

    if len(coor.shape) > 1:
        rois = []
        for x in coor:
            rois += [int(mtx[x[0], x[1], x[2]])]
    else:
        rois = [int(mtx[coor[0], coor[1], coor[2]])]

    return rois


def make_seed_based_networks (CorrMats_Groups, seed, threshold, thresholds_list, netnames, NIIFILEDIR, SAVEFOLDER):
    """ Save nifti files of networks for the Group analyses """

    threshold_index = thresholds_list.index(threshold)

    for i,s in enumerate(seed):

        # Save a nii file with the seed
        seed_network = [0 for i in CorrMats_Groups[0].adjmats[threshold_index][s]]
        seed_network[s-1] = 1
        make_nii_network(NIIFILEDIR, seed_network, SAVEFOLDER + '\\' + netnames[i] + '_seed.nii')

        for y in CorrMats_Groups:
            network = y.adjmats[threshold_index][s-1]

            if y.injury_side == 'X':
                SAVENAME = SAVEFOLDER + '\\' + netnames[i] + '_Group1' + y.injury_side + '.nii'
            elif y.injury_side == 'L':
                SAVENAME = SAVEFOLDER + '\\' + netnames[i] + '_Group2' + y.injury_side + '.nii'
            elif y.injury_side == 'R':
                SAVENAME = SAVEFOLDER + '\\' + netnames[i] + '_Group3' + y.injury_side + '.nii'
            elif y.injury_side == 'N':
                SAVENAME = SAVEFOLDER + '\\' + netnames[i] + '_Group4' + y.injury_side + '.nii'


            make_nii_network(NIIFILEDIR, network, SAVENAME)


