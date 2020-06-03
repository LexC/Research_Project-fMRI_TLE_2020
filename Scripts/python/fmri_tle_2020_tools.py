# -*- coding: utf-8 -*-
"""
DESCRIÇÃO
"""
# %% Libraries

import os
import csv
import numpy as np
import scipy.io as sio
import nibabel as nib

from nilearn import plotting


# %% Classes

class CoorMatrices:
    """
    DESCRIÇÃO
    """
    def __init__(self, matricesfilesdir):
        """
        DESCRIÇÃO
        """
        self.filesdir = matricesfilesdir
        self.raw_data = None
        self.group_by_injury = None
        self.group_by_subject = None


    def get_raw_data(self, subj_info):
        """ Make a variable with all the Correlation matrices and subject information

        matricesfilesdir = list(matrices diretories), from get_matrices_files_dir
        subj_info = list(subjects informations), from get_csv_information
        """
        corrmatrices = []
        for x in self.filesdir:
            corrmatrices += [CoorMatrices.DataAndInfo(x, subj_info)]

        self.raw_data = corrmatrices

        return self

    class DataAndInfo:
        """
        DESCRIÇÃO
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


####### Grouping data by injury  ------------------------------------


    def get_group_by_injury(self, injury_classifications, thresholds=list(np.arange(0.05, 0.95, 0.05))):
        """ DESCRIÇÃO """
        corrmats_groups = []

        for x in injury_classifications:
            corrmats_groups += [CoorMatrices.Injury(self.raw_data, x, thresholds)]

        self.group_by_injury = corrmats_groups

        return self



    class Injury:
        """ DESCRIÇÃO """
        def __init__(self, data, injuryClassification, thresholds):

            self.group = group_data_3darray_injury(data, injuryClassification)
            self.injury_side = injuryClassification

            self.group_mean = fisher_mean(self.group, 0)

            self.adjmats = []
            for r in thresholds:
                self.adjmats += [adjacency_matrices(self.group_mean, r)]


####### Grouping data by volunteer ------------------------------------


    def get_group_by_subject(self, thresholds=list(np.arange(0.05, 0.95, 0.05))):

        i=0
        vollist=[]
        while i<len(self.raw_data):

            aux = [i]

            for j in range(i+1,len(self.raw_data)):

                aux1=self.raw_data[i]
                aux1=[aux1.protocol, aux1.subj_type, aux1.subj_number]

                aux2=self.raw_data[j]
                aux2=[aux2.protocol, aux2.subj_type, aux2.subj_number]

                if aux1 == aux2:
                    aux += [j]

                else:
                    vollist+=[aux]
                    i=j-1
                    break
            i+=1

        subjs=[]
        for x in vollist:
            subjs += [CoorMatrices.Subject(self.raw_data, x, thresholds)]

        self.group_by_subject=subjs
        return self

    class Subject:
        """ DESCRIÇÃO """
        def __init__(self,raw_data,vollist, thresholds):

            self.protocol = raw_data[vollist[0]].protocol
            self.subj_type = raw_data[vollist[0]].subj_type
            self.subj_number = raw_data[vollist[0]].subj_number
            self.injury = raw_data[vollist[0]].injury

            self.runs=group_data_3darray_bynumber(raw_data, vollist)

            self.runs_mean = fisher_mean(self.runs, 0)

            self.adjmats = []
            for r in thresholds:
                self.adjmats += [adjacency_matrices(self.runs_mean, r)]



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


def transform_coor_space(cor, transform_matrix, transform_direction):
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


def get_csv_information(subj_infofile):
    """
    Reads and retuns a list of each line contained in the csv file.

    subj_infofile = str(CSV diretory)

    Output: list(list())
    """
    subj_info = []
    with open(subj_infofile, 'r') as csv_subjinfo:
        for line in csv.reader(csv_subjinfo):
            subj_info += [line]
    return subj_info



####### Grouping data and analyses ------------------------------------


def group_data_3darray_injury(corrmats, injury):
    """
    Group all the Correlation Matrices of a group in one 3D array. If injury is not defined, it will group all matrices.

    corrmats = list(class CorrMap), from load_subjects_correlation_matrices
    injury = str()

    Outuput: np.array(), shape =(3,:,:)
    """
    group = []
    for x in corrmats:
        if x.injury == injury:
            group += [x.map]
    return np.array(group)


def adjacency_matrices(matrix, threshold):
    """ Returns an adjacency matrix based on a threshold"""
    adjmtx = np.array(matrix)
    adjmtx[matrix >= threshold] = 1
    adjmtx[matrix < threshold] = 0
    return adjmtx


def group_data_3darray_bynumber(corrmats, subjnumb):
    """

    """
    group = []
    for i in subjnumb:
        group += [corrmats[i].map]
    return np.array(group)


####### Seed-Based / Nifti functions (.nii) --------------------------------



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
    coor = np.array(coor)

    if coortype == 'mni':
        #coor = mni2cor(coor,niifiledir)
        for i, x in enumerate(coor):
            coor[i] = transform_coor_space(x, img.affine, 0)

    if len(coor.shape) > 1:
        rois = []
        for x in coor:
            rois += [int(mtx[x[0], x[1], x[2]])]
    else:
        rois = [int(mtx[coor[0], coor[1], coor[2]])]

    return rois


def make_seed_based_networks(coormats, seed, threshold, thresholds_list, netnames, NIIFILEDIR, SAVEFOLDER):
    """ Save nifti files of networks for the Group analyses """

    threshold_index = thresholds_list.index(threshold)

    for i, s in enumerate(seed):

        # Save a nii file with the seed
        seed_network = [0 for i in coormats[0].adjmats[threshold_index][s]]
        seed_network[s-1] = 1
        make_nii_network(NIIFILEDIR, seed_network, SAVEFOLDER + '\\' + netnames[i] + '_seed.nii')

        for y in coormats:
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


#%% ####### Grapht Theory ------------------------------------



def graphs_for_all(coormats,MNICOORDIR,SAVEDIR,thresholds,subj_start=0):

    path = os.getcwd()
    os.chdir(SAVEDIR)

    mnicoor = np.array(get_csv_information(MNICOORDIR),dtype=int)

    for i, subj in enumerate(coormats.group_by_subject, ):

        if i >= subj_start:
            if int(subj.subj_number) < 10:
                subj_folder = 'prot'+ str(subj.protocol) + '_' + subj.subj_type + subj.injury + '_0' + str(subj.subj_number)
            else:
                subj_folder = 'prot'+ subj.protocol + '_' + subj.subj_type + subj.injury + '_' + subj.subj_number


            if not os.path.exists(subj_folder):
                os.makedirs(subj_folder)

            print('-'*40)
            print("Ploting ",str(i)," subj "+subj_folder)
            save_graphs_models(subj,subj_folder,mnicoor,thresholds)

    os.chdir(path)

    return None



def save_graphs_models(subj,subj_folder,mnicoor,thresholds):

    for i,r in enumerate(thresholds):
        correlation_matrix=subj.adjmats[i]

        aux=int(round(thresholds[i]*100))

        if aux<10:
            aux = '0'+ str(aux)
        else:
            aux = str(aux)

        print('r='+aux+',',end = ' ')
        SAVEIMAGENAME = subj_folder + "\\" + subj_folder+'_r'+aux+'.png'

        plotting.plot_connectome(correlation_matrix, mnicoor, colorbar=False,
                                  node_color='black',
                                  edge_vmin=0, edge_vmax=1,
                                  node_size=10,display_mode="lzry",
                                  title='r=.'+aux,
                                  output_file=SAVEIMAGENAME)

    return None
