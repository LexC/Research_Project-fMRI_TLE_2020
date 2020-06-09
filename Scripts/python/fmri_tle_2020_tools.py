# -*- coding: utf-8 -*-
"""
DESCRIÇÃO
"""
# %% Libraries

import os
import csv
import numpy as np
import pandas as pd
import nibabel as nib
import networkx as nx
import scipy.io as sio
import matplotlib.pyplot as plt

from nilearn import plotting


# %% Classes

class CoorMatrices:
    """
    DESCRIÇÃO
    """
    def __init__(self, matricesfilesdir,subj_info):
        """
        DESCRIÇÃO
        """
        # self.filesdir = matricesfilesdir
        # self.data = None
        self.group_by_injury = None
        self.group_by_subject = None

        subj_data = {'ID':[],'Protocol':[],'Type':[],'Injury':[],
             'Number':[],'Run':[],'CorrMatrix':[],'FileDir':[]}

        subj_data['FileDir']=matricesfilesdir

        for i,file in enumerate(subj_data['FileDir']):
            # corrmatrices += [CoorMatrices.DataAndInfo(x, subj_info)]

            mat = sio.loadmat(file)
            subj_data['CorrMatrix'] += [mat['map']]


            file = file.split('\\')

            subj_data['Protocol'] += [str(int(file[-7].split('_')[1]))]
            subj_data['Type'] += [file[-6][0:-1]]
            subj_data['Number'] += [str(int(file[-3].split('_')[-2]))]
            subj_data['Run'] += [str(int(file[-3].split('_')[-1]))]

            subj_data['ID'] += [subj_data['Protocol'][i] + '_' + subj_data['Type'][i][0] + subj_data['Number'][i]+ '_' +subj_data['Run'][i]]

            for j in subj_info:
                if j[0:3] == [subj_data['Protocol'][i], subj_data['Type'][i],
                              subj_data['Number'][i]]:
                    subj_data['Injury'] += j[3]

        subj_data=pd.DataFrame(subj_data)
        subj_data=subj_data.set_index('ID')

        self.data = subj_data



####### Grouping data by injury  ------------------------------------

    def get_group_by_injury(self, thresholds=list(np.arange(0.05, 0.95, 0.05))):
        """ DESCRIÇÃO """

        injuryset=list(set(self.data['Injury']))

        group_by_injury = {'Injury':[],'Type':[],'3DMatrix':[],
             'GroupMean':[],'Thresholds':[],'AdjMats':[]}

        for i,injury in enumerate(injuryset):

            group_by_injury['Type'] += [self.data[self.data['Injury']==injury]['Type'][0]]
            group_by_injury['3DMatrix'] += [pandascolumn_to_nparray(self.data[(self.data['Injury']==injury)]['CorrMatrix'])]

            group_by_injury['Injury'] += injury

            group_by_injury['GroupMean'] += [fisher_mean(group_by_injury['3DMatrix'][i], 0)]

            group_by_injury['Thresholds'] += [thresholds]

            adjmats = []
            for r in thresholds:
                adjmats += [adjacency_matrices(group_by_injury['GroupMean'][i], r)]

            group_by_injury['AdjMats'] += [adjmats]

        group_by_injury=pd.DataFrame(group_by_injury)
        group_by_injury=group_by_injury.set_index('Injury')

        self.group_by_injury = group_by_injury

        return self



####### Grouping data by subject ------------------------------------

    def get_group_by_subject(self, thresholds=list(np.arange(0.05, 0.95, 0.05))):


        group_by_subject = {'Protocol':[],'Type':[],'Injury':[],'Number':[],
                           '3DMatrix':[],'RunsMean':[],'Thresholds':[],'AdjMats':[]}

        for prot in set(self.data['Protocol']):
            for typ in set(self.data['Type']):
                for numb in set(self.data['Number']):

                    subj=self.data[(self.data['Protocol'] == prot) & (self.data['Type'] == typ) & (self.data['Number'] == numb)]
                    if subj.empty == False:

                        group_by_subject['Protocol'] += [prot]
                        group_by_subject['Type'] += [typ]
                        group_by_subject['Number'] += [numb]
                        group_by_subject['Injury'] += [subj['Injury'][0]]


                        group_by_subject['3DMatrix'] += [pandascolumn_to_nparray(subj['CorrMatrix'])]

                        group_by_subject['RunsMean'] += [fisher_mean(group_by_subject['3DMatrix'][-1],0)]

                        group_by_subject['Thresholds'] += [thresholds]

                        adjmats = []
                        for r in thresholds:
                            adjmats += [adjacency_matrices(group_by_subject['RunsMean'][-1], r)]

                        group_by_subject['AdjMats'] += [adjmats]

        group_by_subject=pd.DataFrame(group_by_subject)
        self.group_by_subject = group_by_subject

        return self






class LexGraphs:


    def __init__(self, subj_list):


        graphs = {'Graphs_per_r':[]}

        for index, row in subj_list.iterrows():

            graphs_per_threshold = []
            for adjmats in row['AdjMats']:

                aux=tr_zero(adjmats)
                graphs_per_threshold += [nx.from_numpy_array(aux)]
            graphs['Graphs_per_r'] += [graphs_per_threshold]

        graphs = pd.DataFrame(graphs)
        graphs[['Protocol','Type','Injury','Number','Thresholds']] = subj_list[['Protocol','Type','Injury','Number','Thresholds']]

        self.graphs = graphs[['Protocol','Type','Injury','Number','Thresholds',
                              'Graphs_per_r']]


    def get_degrees(self):

        degrees = {'Degrees':[]}

        for index, row in self.graphs.iterrows():

            deg_per_threshold = []
            for graph_r in row['Graphs_per_r']:
                deg_per_threshold += [g_degree(graph_r)]

            degrees['Degrees'] += [deg_per_threshold]

        degrees = pd.DataFrame(degrees)
        self.graphs['Degrees'] = degrees

        return self


    def get_average_degrees(self):

        try:
            self.graphs['Degrees']
        except:
            self = self.get_degrees()

        average_degrees = {'Average_Degrees':[]}

        for index, row in self.graphs.iterrows():

            avdeg_per_r = []
            for deg_r in row['Degrees']:
                avdeg_per_r += [[np.mean(deg_r),np.std(deg_r)]]

            average_degrees['Average_Degrees'] += [avdeg_per_r]

        average_degrees = pd.DataFrame(average_degrees)
        self.graphs['Average_Degrees'] = average_degrees

        return self



# %% Functions

####### Tools ------------------------------------


def change_char(string,old,new):

    new_str=str()
    for i in string:
        if i in old:
            new_str+=new
        else:
            new_str+=i

    return new_str


def pdcolum_normalization(graphs,param):

    graphs=pd.DataFrame(graphs)
    for index, row in graphs.iterrows():
        aux=[]
        for values in row[param]:
            aux += [normalization(values,len(graphs['Graphs_per_r'][index][0].nodes))]
        graphs[param][index]=aux

    return graphs


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


def tr_zero(adjmat):

    for i in range(len(adjmat)):
        adjmat[i,i]=0

    return adjmat


def normalization(array,numb):
    aux = np.array(array)/numb
    return aux.tolist()


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


def pandascolumn_to_nparray(corrmats):
    """
    Group all the Correlation Matrices of a group in one 3D array. If injury is not defined, it will group all matrices.

    corrmats = list(class CorrMap), from load_subjects_correlation_matrices
    injury = str()

    Outuput: np.array(), shape =(3,:,:)
    """
    group = []
    for mtx in corrmats:
            group += [mtx]
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


def make_seed_based_networks(coormats, seed, threshold, netnames, NIIFILEDIR, SAVEFOLDER):
    """ Save nifti files of networks for the Group analyses """

    thresholds_list =coormats['Thresholds'][0]
    threshold_index = thresholds_list.index(threshold)

    for i, s in enumerate(seed):

        # Save a nii file with the seed
        seed_network = [0 for i in range(coormats['GroupMean'][0].shape[1])]
        seed_network[s-1] = 1
        make_nii_network(NIIFILEDIR, seed_network, SAVEFOLDER + '\\' + netnames[i] + '_seed.nii')

        for j,y in enumerate(coormats['AdjMats']):
            network = y[threshold_index][s-1]

            if coormats.index[j] == 'X': # index = injury
                SAVENAME = SAVEFOLDER + '\\' + netnames[i] + '_Group1' + coormats.index[j] + '.nii'
            elif coormats.index[j] == 'L':
                SAVENAME = SAVEFOLDER + '\\' + netnames[i] + '_Group2' + coormats.index[j] + '.nii'
            elif coormats.index[j] == 'R':
                SAVENAME = SAVEFOLDER + '\\' + netnames[i] + '_Group3' + coormats.index[j] + '.nii'
            elif coormats.index[j] == 'N':
                SAVENAME = SAVEFOLDER + '\\' + netnames[i] + '_Group4' + coormats.index[j] + '.nii'


            make_nii_network(NIIFILEDIR, network, SAVENAME)


#%% ####### Grapht Theory ------------------------------------


####### Calculations

def g_degree(G):

    return [deg[1] for deg in G.degree()]


def g_average_degree(G):

    degrees = g_degree(G)

    return [np.mean(degrees),np.std(degrees)]


def g_param_pergroups(param,subjs):

    groups_ids = list(set([subj.injury for subj in subjs]))

    groups_type = [[]]
    groups= [[]]
    for i in range(len(groups_ids)-1):
        groups_type.append([])
        groups.append([])

    for i,subj in enumerate(subjs):

        groups[groups_ids.index(subj.injury)] += [param[i]]
        groups_type[groups_ids.index(subj.injury)] = subj.subj_type

    return list(zip(groups_type,groups_ids,groups))


####### Figures


def save_graphs_models(subj,subj_folder,mnicoor,thresholds):

    for i,r in enumerate(thresholds):
        correlation_matrix=subj[i]

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



def graphs_for_all(coormats,MNICOORDIR,SAVEDIR,subj_start=0):


    path = os.getcwd()
    os.chdir(SAVEDIR)

    mnicoor = np.array(get_csv_information(MNICOORDIR),dtype=int)


    df=coormats[subj_start:]


    for ind in df.index:

        if int(df['Number'][ind]) < 10:
            subj_folder = 'prot'+ str(df['Protocol'][ind]) + '_' + df['Type'][ind] + df['Injury'][ind] + '_0' + str(df['Number'][ind])
        else:
            subj_folder = 'prot'+ df['Protocol'][ind] + '_' + df['Type'][ind] + df['Injury'][ind] + '_' + df['Number'][ind]


        if not os.path.exists(subj_folder):
            os.makedirs(subj_folder)

        print('-'*40)
        print("Ploting ",str(ind)," subj "+subj_folder)
        save_graphs_models(df['AdjMats'][ind],subj_folder,mnicoor, df['Thresholds'][ind])

    os.chdir(path)

    return None



def gfig_parameter_vs_threshold(graphs, param, save_folder, parameter_name,norm = 0):
#(param, thresholds, subjs, save_folder, parameter_name,normalization = 0):
    """

    """

    if norm == 1:

        graphs = pdcolum_normalization(graphs, param)
        lowlim=-0.05
        parameter_name = 'Normalized ' + parameter_name
    else:
        lowlim=-5



    uplim=np.max([np.sum(np.array(row),1) for row in graphs[param]])

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']


    labels_size=30
    for i,y in enumerate(graphs[param]):


        y_value = [j[0] for j in y]
        y_std = [j[1] for j in y]

        thresholds=graphs['Thresholds'][i]

        fig = plt.figure(figsize=(15, 12))
        plt.errorbar(thresholds, y_value, y_std,
                     marker='D', markerfacecolor = 'w', markersize=8,
                     capsize = 5, #ecolor = 'k',
                     ls='dotted', color = colors[0]#,dash_capstyle	='round'
                     )


        plt.xlim(0, 1)
        plt.ylim(lowlim, uplim)
        plt.grid()

        plt.xticks(fontsize=labels_size*0.65)
        plt.yticks(fontsize=labels_size*0.65)

        plt.xlabel('Thershold', fontsize=labels_size)
        plt.ylabel(parameter_name, fontsize=labels_size)


        if graphs['Type'][i] == 'Patient':

            plt.title(graphs['Type'][i]+' '+graphs['Number'][i]+', ' +graphs['Injury'][i]+' Injury',fontsize=labels_size)

            parameter_name=change_char(parameter_name,' ','_')
            if int(graphs['Number'][i]) < 10:
                SAVENAME= save_folder+'\\'+parameter_name+'_'+graphs['Type'][i]+'_'+graphs['Injury'][i]+'_0'+graphs['Number'][i]+'.png'
            else:
                SAVENAME= save_folder+'\\'+parameter_name+'_'+graphs['Type'][i]+'_'+graphs['Injury'][i]+'_'+graphs['Number'][i]+'.png'

        else:

            plt.title(graphs['Type'][i]+' '+graphs['Number'][i],fontsize=labels_size)

            parameter_name=change_char(parameter_name,' ','_')
            if int(graphs['Number'][i]) < 10:
                SAVENAME= save_folder+'\\'+parameter_name+'_'+graphs['Type'][i]+'_0'+graphs['Number'][i]+'.png'
            else:
                SAVENAME= save_folder+'\\'+parameter_name+'_'+graphs['Type'][i]+'_'+graphs['Number'][i]+'.png'



        plt.show()
        # fig.savefig(SAVENAME,format='png')
        plt.close()
        parameter_name=change_char(parameter_name,'_',' ')


    print('-'*40)
    print('All '+parameter_name+' Single Figures Saved !')
    print('-'*40)




def gfig_ofall_parameter_vs_threshold(graphs, param, save_folder, parameter_name,norm = 0):
        #param, thresholds, subjs, save_folder, parameter_name,normalization = 0):

    if norm == 1:
        graphs = pdcolum_normalization(graphs, param)
        lowlim=-0.05
        parameter_name = 'Normalized '+parameter_name
    else:
        lowlim=-5

    thresholds=graphs['Thresholds'][0]


    group_by_injury={'injuryset':list(set(graphs['Injury'])),'arrays':[]}
    for i,injury in enumerate(group_by_injury['injuryset']):
        group_by_injury['arrays'] += [graphs[(graphs['Injury']==injury)][param].tolist()]
    group_by_injury=pd.DataFrame(group_by_injury)

    uplim=np.max([np.sum(np.array(row),1) for row in graphs[param]])


    for index, row in group_by_injury.iterrows():

        labels_size=30
        fig = plt.figure(figsize=(15, 12))

        for i,y in enumerate(row['arrays']):

            y_value = [j[0] for j in y]
            y_std = [j[1] for j in y]


            plt.errorbar(thresholds, y_value, y_std,
                         marker='D', markerfacecolor = 'w', markersize=8,
                         capsize = 5, #ecolor = 'k',
                         ls='dotted'#, color = colors[0]#,dash_capstyle	='round'
                         )

        plt.ylim(lowlim, uplim)
        plt.xlim(0, 1)
        plt.grid()

        plt.xticks(fontsize=labels_size*0.65)
        plt.yticks(fontsize=labels_size*0.65)

        plt.xlabel('Thershold', fontsize=labels_size)
        plt.ylabel(parameter_name, fontsize=labels_size)

        if row['injuryset'] == 'X':
            plt.title('All Volunteers',fontsize=labels_size)
            SAVENAME=save_folder+'\\z_'+parameter_name+'_Volunteerss.png'
        else:
            plt.title('All Patients '+row['injuryset'],fontsize=labels_size)
            SAVENAME=save_folder+'\\z_'+parameter_name+'_Patients'+row['injuryset']+'.png'

        plt.show()
        # fig.savefig(SAVENAME,format='png')
        plt.close()



def gfig_averageofall_parameter_vs_threshold(graphs, param, save_folder, parameter_name,norm = 0):


    if norm == 1:
        graphs = pdcolum_normalization(graphs, param)
        parameter_name = 'Normalized '+parameter_name


    thresholds=graphs['Thresholds'][0]


    group_by_injury={'injuryset':list(set(graphs['Injury'])),'arrays':[]}
    for i,injury in enumerate(group_by_injury['injuryset']):

        group_by_injury['arrays'] += [np.array(graphs[(graphs['Injury']==injury)][param].tolist())]
    group_by_injury=pd.DataFrame(group_by_injury)


    fig = plt.figure(figsize=(15, 12))

    for index, row in group_by_injury.iterrows():

        aux=row['arrays']

        y=[np.mean(row['arrays'][:,i,:]) for i in range(row['arrays'].shape[1])]
        y_std=[np.std(row['arrays'][:,i,:]) for i in range(row['arrays'].shape[1])]

        labels_size=30

        if row['injuryset'] == 'X':
            LABEL = "Volunteers"
        else:
            LABEL = 'Patients' + ' '+ row['injuryset']

        plt.errorbar(thresholds, y, y_std,
                     marker='D', markerfacecolor = 'w', markersize=8,
                     capsize = 5, #ecolor = 'k',
                     ls='dotted',#, color = colors[0]#,dash_capstyle	='round'
                     label = LABEL)


    plt.xlim(0, 1)
    # plt.ylim(-3, 120)
    plt.grid()
    plt.legend(fontsize=labels_size*0.75)

    plt.title('Average of the Groups',fontsize=labels_size)
    plt.xticks(fontsize=labels_size*0.65)
    plt.yticks(fontsize=labels_size*0.65)

    plt.xlabel('Thershold', fontsize=labels_size)
    plt.ylabel(parameter_name, fontsize=labels_size)


    plt.show()

    parameter_name=change_char(parameter_name,' ','_')
    SAVENAME=save_folder+'\\z_'+parameter_name+'_zGroups_mean.png'
    # fig.savefig(SAVENAME,format='png')
    plt.close()



